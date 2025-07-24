local M = {}

---Convenience for plain matching.
---@param a string
---@param b string
---@return integer? start
---@return integer? stop
function string.contains(a, b)
    return a:find(b, 1, true)
end

function M.readtext(path)
    local file = io.open(path, "rb") -- r read mode and b binary mode
    if not file then return nil end
    local content = file:read "*a"   -- *a or *all reads the whole file
    file:close()
    return content
end

-- strip whitespace at both ends of text
function M.strip(text)
    return text:match("^[\t%s]*(.-)[\t%s]*$")
end

---Repeat calls to a given function as many times as the vim count value (default once).
---Optionally pass arguments.
---@param fn function
---@param ... unknown
---@return function
function M.fncount(fn, ...)
    local args = ...
    return function()
        for _ = 1, M.count() do
            fn(args)
        end
    end
end

function M.end_visual()
    -- magic from
    -- https://github.com/neovim/neovim/issues/19770
    -- alternatively call vim.cmd.stopinsert()
    vim.api.nvim_feedkeys('\027', 'xt', false)
end

---@return integer r1 0-index
---@return integer c1 0-index
---@return integer r2 0-index
---@return integer c2 0-index
function M.last_changeyank_range()
    local r1, c1 = unpack(vim.api.nvim_buf_get_mark(0, "["))
    local r2, c2 = unpack(vim.api.nvim_buf_get_mark(0, "]"))
    return r1 - 1, c1, r2 - 1, c2
end

---@return integer r1 0-index
---@return integer c1 0-index
---@return integer r2 0-index
---@return integer c2 0-index
function M.last_visual_range()
    local r1, c1 = unpack(vim.api.nvim_buf_get_mark(0, "<"))
    local r2, c2 = unpack(vim.api.nvim_buf_get_mark(0, ">"))
    return r1 - 1, c1, r2 - 1, c2
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
    local char = vim.api.nvim_buf_get_text(0, r2 - 1, c2, r2 - 1, c2 + 1, {})[1]
    -- if multibyte then char is only half the symbol and won't match a broad pattern like:
    if not char:match('[%w%p%s]') then c2 = c2 + 1 end
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
    vim.api.nvim_win_set_cursor(0, { r1, c1 })
    vim.cmd.normal 'v'
    vim.api.nvim_win_set_cursor(0, { r2, c2 })
end

-- return whether r1, c1 is before r2, c2
function M.before(r1, c1, r2, c2)
    return r1 < r2 or (r1 == r2 and c1 < c2)
end

---@return integer r 0-indexed
---@return integer c 0-indexed
function M.get_cursor()
    local r, c = unpack(vim.api.nvim_win_get_cursor(0))
    return r - 1, c
end

---@param r integer 0-indexed
---@param c integer 0-indexed
function M.set_cursor(r, c)
    vim.api.nvim_win_set_cursor(0, { r + 1, c })
end

---Get all lines of the current buffer.
---@return table array of strings.
function M.get_all_lines()
    return vim.api.nvim_buf_get_lines(0, 0, -1, true)
end

---Is current buffer empty?
---@return boolean empty
function M.is_empty()
    -- only use max 2 lines for efficiency
    local lines = vim.api.nvim_buf_get_lines(0, 0, 2, false)
    return #lines == 1 and lines[1] == ""
end

---Is current buffer named?
---@return boolean named
function M.is_named()
    -- Note that "" evaluates to true in lua
    return vim.api.nvim_buf_get_name(0) ~= ""
end

---@param r integer 0-indexed
function M.get_line(r)
    return vim.api.nvim_buf_get_lines(0, r, r + 1, true)[1]
end

---@param r  integer 0-indexed
---@param c1 integer 0-indexed
---@param c2 integer 0-indexed
---@return string text Get text from a single line
function M.get_text(r, c1, c2)
    return vim.api.nvim_buf_get_text(0, r, c1, r, c2, {})[1]
end

---Get char at r, c (and its beggining byte)
---@param r integer 0-based
---@param c integer 0-based
---@return string char
---@return integer c1 start of byte or multibyte char returned
function M.get_char(r, c)
    if c == 0 then return "", 0 end
    local char = vim.api.nvim_buf_get_text(0, r, c - 1, r, c, {})[1]
    -- if multibyte then char is only half the symbol and won't match a broad pattern like:
    if char:match('[%w%p%s]') then
        return char, c - 1
    else
        return vim.api.nvim_buf_get_text(0, r, c - 2, r, c, {})[1], c - 2
    end
end

---Get char right before cursor, i.e. most recently typed.
---@return string char
---@return integer c1 starting byte of char
function M.get_current_char()
    local r, c = M.get_cursor()
    return M.get_char(r, c)
end

---Write char right after cursor.
---@param char string
function M.put_char(char)
    vim.api.nvim_put({ char }, "c", false, true)
end

---Shortcut for the most typical call to get nvim internal represenations of
---keycodes given a string, e.g. "<C-V>".
---@param str string
function M.nvim_code(str)
    return vim.api.nvim_replace_termcodes(str, true, false, true)
end

---Press keys
---@param keys string e.g. "<C-^>"
---@param opts table remap: set to false to make sure key is pressed like in stock vim.
function M.press(keys, opts)
    local mode
    if opts and (opts.noremap or opts.remap == false) then
        mode = 'n'
    else
        mode = 'm'
    end
    vim.api.nvim_feedkeys(M.nvim_code(keys), mode, false)
end

---Get current mode.
---@return string mode "n", "v", "V", ^V, etc.
function M.get_mode()
    return vim.api.nvim_get_mode().mode
end

---Set mode
---@param mode string "i", "n", "v", or "V", or ^V
function M.set_mode(mode)
    if vim.startswith(mode, "i") then
        vim.cmd.startinsert()
    elseif vim.startswith(mode, "n") then
        -- also stops visual mode
        vim.cmd.stopinsert()
    else
        vim.cmd.normal(mode)
    end
end

M.ctrl_v = M.nvim_code("<C-v>")
---Check if current mode is visual blockwise.
---@return boolean
function M.is_visual_blockwise()
    return M.get_mode() == M.ctrl_v
end

---Notify of stderr or stdout from a vim.system call obj.
---Intended for stderr. Uses stdout if stderr is empty, which may be the case if program doesn't utilise stderr at all.
---@param obj any
function M.schedule_notify(obj)
    local text = obj.stderr:gsub("\n$", "")
    if text == "" then
        text = obj.stdout:gsub("\n$", "")
    end
    vim.schedule(function() -- notify when we are ready
        vim.notify(text)    -- vim.notify instead of print to see multiple lines
    end)
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
