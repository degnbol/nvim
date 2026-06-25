local M = {}

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
    for _, pat in ipairs({
        "%$%b()[%w/%._%-+~=@#%%]*",
        "%$[%w_]+[%w/%._%-+~=@#%%]*",
    }) do
        local start = 1
        while true do
            local s, e = line:find(pat, start)
            if not s then break end
            if col >= s and col <= e then
                match = line:sub(s, e)
                break
            end
            start = e + 1
        end
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

-- Resolve the path token under the cursor to an existing path, or nil.
-- Resolution order, first hit wins:
--   1. git-root / $VAR / buffer NAME=value expansion (dollar_token above).
--   2. Plain <cfile> with ~/$ENV expanded, if absolute & a readable file.
--   3. <cfile> joined with the buffer's directory.
--   4. findfile(<cfile>) — honours 'path' / 'suffixesadd'.
-- A trailing :lnum or #anchor is stripped before the readability test.
---@param bufnr integer|nil 0 or nil = current buffer
---@return string|nil path existing path under cursor (file, or dir via case 1)
function M.resolve_path_under_cursor(bufnr)
    if not bufnr or bufnr == 0 then
        bufnr = vim.api.nvim_get_current_buf()
    end

    local token = dollar_token(bufnr)
    if token then
        local path = vim.fn.expand(token)
        if vim.uv.fs_stat(path) then return path end
    end

    local cfile = vim.fn.expand("<cfile>")
    if cfile == "" then return nil end
    -- Cheap shape gate: bare identifiers (the common LSP-hover case) never
    -- resolve to a file, so skip the filesystem walk entirely.
    if not cfile:match("[/~]") and not cfile:match("%.%w+$") then return nil end
    -- Strip a trailing :lnum or #anchor so foo.lua:42 / README.md#install still
    -- test readable. ponytail: jump-to-line/anchor deferred (YAGNI).
    cfile = cfile:gsub("#%S*$", ""):gsub(":%d+$", "")

    local function is_file(p)
        local stat = vim.uv.fs_stat(p)
        return stat ~= nil and stat.type == "file"
    end

    local expanded = vim.fn.expand(cfile) --[[@as string]] -- ~ and $ENV
    if vim.startswith(expanded, "/") and is_file(expanded) then
        return expanded
    end
    local bufdir = vim.fn.fnamemodify(vim.api.nvim_buf_get_name(bufnr), ":h")
    local joined = bufdir .. "/" .. expanded
    if is_file(joined) then
        return vim.fn.fnamemodify(joined, ":p")
    end
    local found = vim.fn.findfile(expanded) -- honours 'path'/'suffixesadd'
    if found ~= "" then return vim.fn.fnamemodify(found, ":p") end
    return nil
end

return M
