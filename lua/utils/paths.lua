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

-- Find the $(...) or $VAR token under the cursor and expand the git-root
-- subshells ($(git root) / $(git rev-parse ...)), $ENV, and buffer-local
-- NAME=value assignments. Returns the expanded token (ready for vim.fn.expand)
-- or nil if no such token sits under the cursor. Arbitrary command
-- substitution is intentionally NOT run — only the git-root forms, resolved via
-- vim.fs.root — since executing shell from a file would be a footgun.
local function dollar_token(bufnr)
    local line = vim.api.nvim_get_current_line()
    local col = vim.api.nvim_win_get_cursor(0)[2] + 1
    local match
    -- Tail class differs from PATH_RUN on purpose: no $, or $A/$B would be one
    -- token.
    for _, pat in ipairs({
        "%$%b()[%w/%._%-+~=@#%%]*",
        "%$[%w_]+[%w/%._%-+~=@#%%]*",
    }) do
        match = util.match_covering(line, col, pat)
        if match then break end
    end
    if not match then return nil end
    local token = match:gsub("%$%b()", function(m)
        local cmd = m:sub(3, -2):match("^%s*(.-)%s*$")
        if cmd == "git root" or cmd:match("^git rev%-parse%s+%-%-show%-toplevel$") then
            return M.git_root(bufnr) or vim.env.ROOT or m
        end
        return m
    end)
    -- Env var first; only fall back to buffer assignment if unset.
    token = token:gsub("%$([%w_]+)", function(name)
        if vim.env[name] then return "$" .. name end
        return M.buffer_var(name, bufnr) or ("$" .. name)
    end)
    return token
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

-- Resolve a path token to an existing file: absolute, then joined with the
-- buffer's directory, then findfile() — which honours 'path'/'suffixesadd'.
---@param token string ~ and $ENV are expanded
---@param bufnr integer
---@return string|nil path
local function resolve_file(token, bufnr)
    local expanded = vim.fn.expand(token) --[[@as string]]
    if vim.startswith(expanded, "/") and M.is_file(expanded) then
        return expanded
    end
    local joined = vim.fs.joinpath(vim.fs.dirname(vim.api.nvim_buf_get_name(bufnr)), expanded)
    -- fnamemodify, not vim.fs.normalize: `..` must apply after following a
    -- symlink, as the kernel does, not remove the segment before it.
    if M.is_file(joined) then
        return vim.fn.fnamemodify(joined, ":p")
    end
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

-- Resolve the location under the cursor: an existing path, plus the line and
-- column of a `:lnum` / `:lnum:col` suffix where the token carries one.
-- Resolution order, first hit wins:
--   1. git-root / $VAR / buffer NAME=value expansion (dollar_token above).
--   2. The path half of the `path:lnum[:col]` token covering the cursor.
--   3. <cfile> with ~/$ENV expanded (no suffix: 'isfname' has already cut it).
-- A trailing #anchor is stripped before the readability test.
---@param bufnr integer|nil 0 or nil = current buffer
---@return string|nil path existing path under cursor (file, or dir via case 1)
---@return integer|nil lnum 1-based, nil unless the token carried a suffix
---@return integer|nil col 1-based
function M.resolve_location_under_cursor(bufnr)
    if not bufnr or bufnr == 0 then
        bufnr = vim.api.nvim_get_current_buf()
    end

    -- Parsed up front so case 1 gets the suffix too: dollar_token's own
    -- character class stops at the ':', so it never sees one.
    local suffixed, lnum, col = parse_suffixed_token()

    local token = dollar_token(bufnr)
    if token then
        local path = vim.fn.expand(token)
        if vim.uv.fs_stat(path) then return path, lnum, col end
    end

    if suffixed then
        -- Only when it resolves: PATH_RUN is narrower than 'isfname', so the
        -- token can be a tail of the real path (`/tmp/a,b.py` -> `b.py`) and
        -- <cfile> below is then the better seed. Residual: a tail that happens
        -- to exist relative to the buffer's directory wins wrongly.
        local path = resolve_file(suffixed, bufnr)
        if path then return path, lnum, col end
    end

    local cfile = vim.fn.expand("<cfile>") --[[@as string]]
    if cfile == "" then return nil end
    -- Cheap shape gate: bare identifiers (the common LSP-hover case) never
    -- resolve to a file, so skip the filesystem walk entirely.
    if not cfile:match("[/~]") and not cfile:match("%.%w+$") then return nil end
    -- Strip a trailing #anchor so README.md#install still tests readable.
    return resolve_file(cfile:gsub("#%S*$", ""), bufnr)
end

-- Resolve a path token to its single absolute candidate, as pure string math —
-- no filesystem access. vim.fs.normalize expands ~ / $ENV and resolves . / ..;
-- a relative token is anchored to the buffer's own directory when it has a
-- name, else the cwd (the nameless-scratch-buffer case). Unlike
-- resolve_location_under_cursor (the gf resolver) there is no 'path' search, no
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
    local name = vim.api.nvim_buf_get_name(bufnr)
    local base = name ~= "" and vim.fs.dirname(name) or vim.uv.cwd()
    return vim.fs.normalize(vim.fs.joinpath(base, path))
end

return M
