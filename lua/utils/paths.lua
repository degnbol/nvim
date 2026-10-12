local util = require "utils/init"

local M = {}

-- A maximal run of path characters. Narrower than 'isfname', which also
-- includes `,`.
M.PATH_RUN = "[%w/%._%-+~=@#%%$]+"

-- Longest first. A leading 0 is excluded: line and column 0 do not exist in the
-- 1-based suffix convention, and nvim_win_set_cursor errors on row 0.
local SUFFIXED = {
    M.PATH_RUN .. ":[1-9]%d*:[1-9]%d*",
    M.PATH_RUN .. ":[1-9]%d*",
}

-- A longer run is an identifier, not a line number — and past 2^53 it stops
-- being an integral double, which vim rejects (E805: Using a Float as a
-- Number) wherever the value reaches a builtin.
local MAX_DIGITS = 9

-- Git project root. vim.fs.root walks the tree in pure Lua — no subprocess.
-- `source` is a bufnr or path. Matches .git as file or directory, so inside a
-- submodule it returns the submodule root (same as the git toplevel command).
function M.git_root(source)
    return vim.fs.root(source, ".git")
end

-- Look up a shell/Makefile-style `NAME=value` assignment in the buffer.
-- Matches line-anchored assignments with no spaces around `=` (zsh/bash/make
-- convention), strips surrounding single or double quotes. Returns nil if not
-- found. Scans from the bottom so later assignments win.
function M.buffer_var(name, bufnr)
    local lines = vim.api.nvim_buf_get_lines(bufnr, 0, -1, false)
    for i = #lines, 1, -1 do
        local val = lines[i]:match("^%s*" .. name .. "=(.+)")
        if val then
            val = val:match("^%s*(.-)%s*$")
            return val:match("^['\"](.-)['\"]$") or val
        end
    end
end

--- Expand `$(git root)` and `$(git rev-parse --show-toplevel)` to the git root
--- (else $ROOT), then each $NAME from the environment, else from a NAME=value
--- assignment in the buffer. Unresolvable parts stay literal. Arbitrary command
--- substitution is intentionally NOT run: executing shell from a file would be
--- a footgun.
---@param token string
---@param bufnr integer buffer for the git root and NAME=value lookups
---@return string
function M.expand_dollars(token, bufnr)
    token = token:gsub("%$%b()", function(m)
        local cmd = m:sub(3, -2):match("^%s*(.-)%s*$")
        if cmd == "git root" or cmd:match("^git rev%-parse%s+%-%-show%-toplevel$") then
            return M.git_root(bufnr) or vim.env.ROOT or m
        end
        return m
    end)
    return (token:gsub("%$([%w_]+)", function(name)
        return vim.env[name] or M.buffer_var(name, bufnr) or ("$" .. name)
    end))
end

-- The $(...) or $VAR token under the cursor, or nil if none sits there.
local function dollar_token()
    local line = vim.api.nvim_get_current_line()
    local col = vim.api.nvim_win_get_cursor(0)[2] + 1
    -- Tail class differs from PATH_RUN on purpose: no $, or $A/$B would be one
    -- token.
    for _, pat in ipairs({
        "%$%b()[%w/%._%-+~=@#%%]*",
        "%$[%w_]+[%w/%._%-+~=@#%%]*",
    }) do
        local match = util.match_covering(line, col, pat)
        if match then return match end
    end
end

---Whether `path` is a regular file, following symlinks.
---@param path string
---@return boolean
function M.is_file(path)
    local stat = vim.uv.fs_stat(path)
    return stat ~= nil and stat.type == "file"
end

---Whether `path` is a directory, following symlinks.
---@param path string
---@return boolean
function M.is_dir(path)
    local stat = vim.uv.fs_stat(path)
    return stat ~= nil and stat.type == "directory"
end

-- The directory a buffer's relative paths resolve from: its file's directory,
-- or the cwd for a buffer with no file — nameless, or named by a URL
-- (`scheme://…`), whose dirname is no directory.
---@param bufnr integer 0 = current buffer
---@return string
function M.buf_dir(bufnr)
    local name = vim.api.nvim_buf_get_name(bufnr)
    if name == "" or name:match("^%a[%w+.-]*://") then return vim.uv.cwd() end
    return vim.fs.dirname(name)
end

-- `path` with a leading `~` expanded to $HOME. Unlike expand(), it never runs a
-- `backtick` shell command, which a buffer NAME=value can carry.
---@param path string
---@return string
local function expand_tilde(path)
    if path == "~" or vim.startswith(path, "~/") then return vim.env.HOME .. path:sub(2) end
    return path
end

--- The regular file at relative path `name` in the buffer's directory
--- (M.buf_dir), as an absolute path; nil when there is none or `name` is
--- absolute.
---@param name string
---@param bufnr integer 0 = current buffer
---@return string|nil path
function M.file_in_buf_dir(name, bufnr)
    if vim.startswith(name, "/") then return nil end
    local joined = vim.fs.joinpath(M.buf_dir(bufnr), name)
    -- fnamemodify, not vim.fs.normalize: `..` must apply after following a
    -- symlink, as the kernel does, not remove the segment before it.
    if M.is_file(joined) then return vim.fn.fnamemodify(joined, ":p") end
    return nil
end

-- Resolve a path token to an existing file: absolute, then joined with the
-- buffer's directory, then findfile() — which honours 'path'/'suffixesadd'.
---@param token string ~, $(git root) and $NAME are expanded (M.expand_dollars)
---@param bufnr integer
---@return string|nil path
local function resolve_file(token, bufnr)
    local expanded = expand_tilde(M.expand_dollars(token, bufnr))
    if vim.startswith(expanded, "/") and M.is_file(expanded) then
        return expanded
    end
    local sibling = M.file_in_buf_dir(expanded, bufnr)
    if sibling then return sibling end
    -- findfile hands back a URL untouched, without a filesystem lookup, so the
    -- same existence test the branches above use has to gate this one too.
    local found = vim.fn.findfile(expanded) --[[@as string]]
    found = found ~= "" and vim.fn.fnamemodify(found, ":p") or ""
    if M.is_file(found) then return found end
    return nil
end

-- Split the `path:lnum[:col]` token covering the cursor. Only this form is
-- parsed — gF's looser separators (`f (30)`, `f @ 20`, `f line 10`) invite
-- false positives in prose. <cfile> cannot serve here: ':' is not in 'isfname'
-- on unix, so the suffix is cut off before it is ever seen (and with the cursor
-- on the digits, <cfile> is the bare number).
---@return string|nil path the token with its suffix removed
---@return integer|nil lnum 1-based
---@return integer|nil col 1-based
local function parse_suffixed_token()
    local line = vim.api.nvim_get_current_line()
    local pos = vim.api.nvim_win_get_cursor(0)[2] + 1
    for _, pat in ipairs(SUFFIXED) do
        local match = util.match_covering(line, pos, pat)
        if match then
            local path, lnum, col = match:match("^(.*):(%d+):(%d+)$")
            if not path then
                path, lnum = match:match("^(.*):(%d+)$")
            end
            -- An over-long run is an identifier (abc:12345678901234567890),
            -- not a line number. Falling through lets the shorter pattern
            -- still salvage a valid :lnum.
            if #lnum <= MAX_DIGITS and (not col or #col <= MAX_DIGITS) then
                return path, tonumber(lnum), col and tonumber(col)
            end
        end
    end
    return nil
end

-- Whether `name` exists where gf would look: an absolute path as is, a relative
-- one through findfile()/finddir(), which search 'path' (fs_stat would use the
-- cwd).
---@param name string
---@return boolean
local function exists_for_gf(name)
    if vim.startswith(name, "/") then return vim.uv.fs_stat(name) ~= nil end
    return vim.fn.findfile(name) ~= "" or vim.fn.finddir(name) ~= ""
end

--- The first of these candidates that differs from `fname` and exists, as an
--- absolute path or on 'path': `fname` with its `$(git root)` / $NAME forms
--- expanded, then that with a trailing `#anchor` also removed.
---@param fname string
---@param bufnr integer buffer for the git root and NAME=value lookups
---@return string|nil name nil when no candidate qualifies
function M.transform_name(fname, bufnr)
    local expanded = M.expand_dollars(fname, bufnr)
    for _, name in ipairs({ expanded, (expanded:gsub("#%S*$", "")) }) do
        if name ~= fname and exists_for_gf(name) then return name end
    end
end

--- An 'includeexpr' that chains: the current buffer's b:ft_includeexpr (an
--- 'includeexpr' value) when its result exists, else `fname` with $ forms
--- expanded or a trailing `#anchor` removed when that exists, else `fname`.
--- Errors from b:ft_includeexpr propagate.
---@param fname string v:fname, which b:ft_includeexpr reads directly
---@return string
function M.includeexpr(fname)
    local bufnr = vim.api.nvim_get_current_buf()
    local ft_expr = vim.b[bufnr].ft_includeexpr
    if ft_expr then
        local name = tostring(vim.fn.eval(ft_expr))
        if exists_for_gf(name) then return name end
    end
    return M.transform_name(fname, bufnr) or fname
end

M.INCLUDEEXPR = "v:lua.require'utils.paths'.includeexpr(v:fname)"

--- Location from the token under the cursor in the current window, first hit
--- wins:
---   1. The $(...) / $VAR token, with `$(git root)`, $NAME and ~ expanded.
---   2. The path half of the `path:lnum[:col]` token.
---@return string|nil path existing path (file, or dir via case 1)
---@return integer|nil lnum 1-based, nil unless the token carried a suffix
---@return integer|nil col 1-based
function M.cursor_location()
    -- Parsed up front so case 1 gets the suffix too: dollar_token's own
    -- character class stops at the ':', so it never sees one.
    local suffixed, lnum, col = parse_suffixed_token()

    local token = dollar_token()
    if token then
        local path = expand_tilde(M.expand_dollars(token, 0))
        if vim.uv.fs_stat(path) then return path, lnum, col end
    end

    if suffixed then
        -- Only when it resolves: PATH_RUN is narrower than 'isfname', so the
        -- token can be a tail of the real path (`/tmp/a,b.py` -> `b.py`), and
        -- <cfile> is then the better seed for whatever resolves next.
        -- Residual: a tail that happens to exist relative to the buffer's
        -- directory wins wrongly.
        local path = resolve_file(suffixed, 0)
        if path then return path, lnum, col end
    end
    return nil
end

--- Resolve the location under the cursor in the current window:
--- M.cursor_location, else <cfile> with ~/$ENV expanded (no suffix: 'isfname'
--- has already cut it), as is, else with buffer $NAME forms expanded or a
--- trailing `#anchor` removed.
---@return string|nil path existing path under cursor (file, or dir from a $ token)
---@return integer|nil lnum 1-based, nil unless the token carried a suffix
---@return integer|nil col 1-based
function M.resolve_location_under_cursor()
    local path, lnum, col = M.cursor_location()
    if path then return path, lnum, col end

    local cfile = vim.fn.expand("<cfile>") --[[@as string]]
    if cfile == "" then return nil end
    -- Cheap shape gate: bare identifiers (the common LSP-hover case) never
    -- resolve to a file, so skip the filesystem walk entirely.
    if not cfile:match("[/~]") and not cfile:match("%.%w+$") then return nil end
    path = resolve_file(cfile, 0)
    if path then return path end
    local name = M.transform_name(cfile, 0)
    return name and resolve_file(name, 0)
end

-- Resolve a path token to its single absolute candidate, as pure string math —
-- no filesystem access. vim.fs.normalize expands ~ / $ENV and resolves . / ..;
-- a relative token is anchored to M.buf_dir. Unlike
-- resolve_location_under_cursor there is no 'path' search, no
-- git-root or buffer-var expansion, and no cursor dependency, so the result is
-- a pure function of (token, base dir) — a valid cache key for an existence
-- check.
---@param token string a path-shaped token
---@param bufnr integer 0 = current buffer
---@return string|nil abspath normalized absolute path (existence not checked)
function M.resolve_path(token, bufnr)
    if token == "" then return nil end
    local path = vim.fs.normalize(token)
    if vim.startswith(path, "/") then return path end
    return vim.fs.normalize(vim.fs.joinpath(M.buf_dir(bufnr), path))
end

return M
