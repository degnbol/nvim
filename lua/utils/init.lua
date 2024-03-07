#!/usr/bin/env lua
M = {}

-- strip whitespace at both ends of text
function M.strip(text)
    return text:match("^[\t%s]*(.-)[\t%s]*$")
end

function M.end_visual()
    -- magic from
    -- https://github.com/neovim/neovim/issues/19770
    vim.api.nvim_feedkeys('\027', 'xt', false)
end

---@return r1 int 0-index
---@return c1 int 0-index
---@return r2 int 0-index
---@return c2 int 0-index
function M.last_changeyank_range()
    local r1, c1 = unpack(vim.api.nvim_buf_get_mark(0, "["))
    local r2, c2 = unpack(vim.api.nvim_buf_get_mark(0, "]"))
    return r1-1, c1, r2-1, c2
end
---@return r1 int 0-index
---@return c1 int 0-index
---@return r2 int 0-index
---@return c2 int 0-index
function M.last_visual_range()
    local r1, c1 = unpack(vim.api.nvim_buf_get_mark(0, "["))
    local r2, c2 = unpack(vim.api.nvim_buf_get_mark(0, "]"))
    return r1-1, c1, r2-1, c2
end

-- Return (1,0)-indexed start and end position of visual selection:
-- r1, c1, r2, c2
function M.get_visual_range()
    -- gives wrong coordinates if we don't end visual first.
    -- If you want to have visual mode unaffected then call gv() below.
    M.end_visual()
    local r1, c1 = unpack(vim.api.nvim_buf_get_mark(0, "<"))
    local r2, c2 = unpack(vim.api.nvim_buf_get_mark(0, ">"))
    -- don't allow selection beyond line
    c2 = math.min(c2, #vim.fn.getline(r2))
    -- handle edge-case where final char is unicode or other multibyte char
    local char = vim.api.nvim_buf_get_text(0, r2-1, c2, r2-1, c2+1, {})[1]
    -- if multibyte then char is only half the symbol and won't match a broad pattern like:
    if not char:match('[%w%p%s]') then c2=c2+1 end
    return r1, c1, r2, c2
end

--- vim.cmd.normal 'gv' doesn't seem to work. This function tries to do a simple gv.
--- Limitations: It can currently only place cursor at the correct end if 
--- cursor hasn't been moved since visual mode was ended. Otherwise it defaults 
--- to placing cursor and end of selection.
function M.gv()
    local r1, c1 = unpack(vim.api.nvim_buf_get_mark(0, '<'))
    local r2, c2 = unpack(vim.api.nvim_buf_get_mark(0, '>'))
    local r, c = unpack(vim.api.nvim_win_get_cursor(0))
    if r == r1 and c == c1 then
        r1 = r2
        c1 = c2
        r2 = r
        c2 = c
    end
    vim.api.nvim_win_set_cursor(0, {r1, c1})
    vim.cmd.normal 'v'
    vim.api.nvim_win_set_cursor(0, {r2, c2})
end

-- return whether r1, c1 is before r2, c2
function M.before(r1, c1, r2, c2)
    return r1 < r2 or (r1 == r2 and c1 < c2)
end

---@return r, c both 0-indexed
function M.get_cursor()
    local r, c = unpack(vim.api.nvim_win_get_cursor(0))
    return r-1, c
end

---@r int 0-indexed
function M.get_line(r)
    return vim.api.nvim_buf_get_lines(0, r, r+1, true)[1]
end

---@r int 0-indexed
---@c1 int 0-indexed
---@c2 int 0-indexed
---@return get text from a single line
function M.get_text(r, c1, c2)
    return vim.api.nvim_buf_get_text(0, r, c1, r, c2, {})[1]
end

function M.get_mode()
    return vim.api.nvim_get_mode().mode
end

-- open with the OS-specific shell command
local opener
if vim.fn.has("macunix") == 1 then
    opener = "open"
elseif vim.fn.has("linux") == 1 then
    opener = "xdg-open"
elseif vim.fn.has("win64") == 1 or vim.fn.has("win32") == 1 then
    opener = "start"
end
function M.open(urlOrPath)
    local openCommand = string.format("%s '%s' >/dev/null 2>&1", opener, urlOrPath)
    os.execute(openCommand)
end

return M
