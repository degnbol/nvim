local M = {}

---Read a whole file.
---@param path string
---@return string|nil content nil if the file can't be opened
---@return string|nil err why the file can't be opened
function M.readtext(path)
    local file, err = io.open(path, "rb") -- r read mode and b binary mode
    if not file then return nil, err end
    local content = file:read "*a"   -- *a or *all reads the whole file
    file:close()
    return content
end

---Scan `s` left to right for non-overlapping matches of `pattern` and return
---the one covering `pos`, or nil if none does. A match starting inside an
---earlier match is never reached — as in `string.gmatch`.
---@param s string
---@param pos integer 1-indexed byte position
---@param pattern string lua pattern
---@return string|nil match
function M.match_covering(s, pos, pattern)
    local start = 1
    while true do
        local s1, e1 = s:find(pattern, start)
        if not s1 then return nil end
        if pos >= s1 and pos <= e1 then return s:sub(s1, e1) end
        start = math.max(s1, e1) + 1 -- max: an empty match ends before it starts
    end
end

---Repeat calls to a given function as many times as the vim count value (default once).
---Optionally pass arguments.
---@param fn function
---@param ... unknown
---@return function
function M.fncount(fn, ...)
    local args = { ... }
    return function()
        for _ = 1, vim.v.count1 do
            fn(unpack(args))
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

---Start and end of the last visual selection. Ends visual mode, call `gv()`
---afterwards to reselect. A `V` selection spans whole lines. For a `<C-v>`
---selection, the start is that of the first line's part of the block and the
---end is that of the last line's part, which is short of the block's right
---edge when the last line is short or the block extends with `$`.
---@return integer r1 1-indexed line of the first char
---@return integer c1 0-indexed first byte of the first char
---@return integer r2 1-indexed line of the last char
---@return integer c2 0-indexed last byte of the last char, -1 on an empty line
function M.get_visual_range()
    -- gives wrong coordinates if we don't end visual first.
    M.end_visual()
    local mode = vim.fn.visualmode()
    -- visualmode() is "" before the first visual selection, which getregionpos rejects.
    local segments = vim.fn.getregionpos(vim.fn.getpos("'<"), vim.fn.getpos("'>"), { type = mode == "" and "v" or mode })
    -- Columns are 1-indexed bytes, and 0 on an empty line.
    local first, last = segments[1][1], segments[#segments][2]
    return first[2], math.max(0, first[3] - 1), last[2], last[3] - 1
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

---Set cursor column without changing row.
---@param c integer 0-indexed
function M.set_col(c)
    local r, _ = unpack(vim.api.nvim_win_get_cursor(0))
    vim.api.nvim_win_set_cursor(0, { r, c })
end

---Set the view of a window. Wraps `winrestview()`; fields not listed
---below pass through unchanged.
---
---When `topline` is set but `lnum` is not, the cursor is placed at the
---closest line within scrolloff of the new viewport — vim would
---otherwise drift `topline` forward to satisfy its margin constraints,
---silently overriding the requested topline. The cursor row stays as
---close as possible to its current value.
---
---Omitted fields keep their current values: cursor column (`col`),
---horizontal offset (`leftcol`/`skipcol`), diff filler (`topfill`),
---etc. Set them explicitly to override.
---
---All line/column fields use winrestview's convention: lines are
---1-indexed, columns are 0-indexed byte offsets.
---
---Field descriptions:
---  winid: window (0 or nil = current).
---  topline: 1-indexed buffer line to place at viewport top.
---  lnum: cursor line (1-indexed). When set, suppresses auto-clamping.
---  col: cursor byte column (0-indexed).
---  coladd: 'virtualedit' offset past end of line.
---  curswant: preferred column for vertical motion.
---  leftcol: leftmost visible column when 'wrap' is off.
---  skipcol: columns skipped on a wrapped long line.
---  topfill: diff filler line count.
---  silent: skip autocmds (WinScrolled, CursorMoved, ...).
---@param opts {
---  winid: integer?,
---  topline: integer?,
---  lnum: integer?,
---  col: integer?,
---  coladd: integer?,
---  curswant: integer?,
---  leftcol: integer?,
---  skipcol: integer?,
---  topfill: integer?,
---  silent: boolean?,
---}
function M.set_view(opts)
    opts = opts or {}
    local winid = opts.winid
    if not winid or winid == 0 then
        winid = vim.api.nvim_get_current_win()
    end
    if not vim.api.nvim_win_is_valid(winid) then
        return
    end

    local view = {
        topline = opts.topline,
        lnum = opts.lnum,
        col = opts.col,
        coladd = opts.coladd,
        curswant = opts.curswant,
        leftcol = opts.leftcol,
        skipcol = opts.skipcol,
        topfill = opts.topfill,
    }

    if opts.topline and not opts.lnum then
        local scrolloff = math.max(
            0,
            vim.api.nvim_get_option_value("scrolloff", { win = winid })
        )
        local height = vim.api.nvim_win_get_height(winid)
        local last = vim.api.nvim_buf_line_count(
            vim.api.nvim_win_get_buf(winid)
        )
        local current = vim.api.nvim_win_get_cursor(winid)[1]
        local lo = math.min(opts.topline + scrolloff, last)
        local hi = math.min(opts.topline + height - 1 - scrolloff, last)
        -- Pathological: window height < 2 * scrolloff + 1 (no valid range).
        if hi < lo then
            hi = lo
        end
        view.lnum = math.max(lo, math.min(current, hi))
    end

    local function apply()
        vim.fn.winrestview(view)
    end

    if opts.silent then
        local saved = vim.o.eventignore
        vim.o.eventignore = "all"
        local ok, err = pcall(vim.api.nvim_win_call, winid, apply)
        vim.o.eventignore = saved
        if not ok then
            error(err)
        end
    else
        vim.api.nvim_win_call(winid, apply)
    end
end

---Add current cursor position to the jumplist.
---Useful for setting cursor and having the change function as a jump.
function M.jumplist_add()
    vim.cmd.normal { "m`", bang = true }
end

---Show a file in the current window, loading it if it is not loaded yet.
---`bufadd` maps the path onto the buffer already holding it where there is one,
---resolving symlinks on the way (/tmp/x finds a buffer named /private/tmp/x),
---so the identity test below is exact and a file is never opened twice. It does
---not expand `~` — hence the normalize, which `:edit` would have done itself.
---@param filepath string
function M.edit(filepath)
    local bufnr = vim.fn.bufadd(vim.fs.normalize(filepath))
    vim.bo[bufnr].buflisted = true -- bufadd() leaves a new buffer unlisted
    -- Re-showing the current buffer is not free: it drops the cursor column.
    if bufnr ~= vim.api.nvim_get_current_buf() then
        vim.api.nvim_win_set_buf(0, bufnr)
    end
end

---Jump to position in specific file, leaving the departure in the jumplist.
---@param filepath string
---@param row integer 0-indexed. Clamped to the buffer's lines, as `:edit +N` does.
---@param col integer 0-indexed. Negative is clamped to 0.
function M.jump(filepath, row, col)
    M.jumplist_add()
    M.edit(filepath)
    M.set_cursor(math.max(0, math.min(row, vim.api.nvim_buf_line_count(0) - 1)), math.max(0, col))
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

---Whether the `width` display cells starting at byte `col` of buffer row `row`
---cross a screen-line wrap in `win`. Works for off-screen lines.
---@param win integer
---@param row integer 0-indexed
---@param col integer 0-indexed byte column
---@param width integer display cells
---@return boolean
function M.crosses_wrap(win, row, col, width)
    local s = vim.fn.virtcol({ row + 1, col + 1 }, true, win)[1] - 1
    return vim.api.nvim_win_text_height(win, {
        start_row = row, end_row = row, start_vcol = s, end_vcol = s + width,
    }).all > 1
end

---The char left of column `c`, i.e. the one covering byte `c - 1`.
---@param r integer 0-indexed line
---@param c integer 0-indexed byte column
---@return string char "" at column 0
---@return integer c1 0-indexed first byte of `char`
function M.get_char(r, c)
    if c == 0 then return "", 0 end
    local c1 = c - 1 + vim.str_utf_start(M.get_line(r), c)
    return M.get_text(r, c1, c), c1
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

---Press keys
---@param keys string e.g. "<C-^>"
---@param opts table|nil remap: set to false to make sure key is pressed like in stock vim.
function M.press(keys, opts)
    local mode
    if opts and (opts.noremap or opts.remap == false) then
        mode = 'n'
    else
        mode = 'm'
    end
    vim.api.nvim_feedkeys(vim.keycode(keys), mode, false)
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

local ctrl_v = vim.keycode("<C-v>")
---Check if current mode is visual blockwise.
---@return boolean
function M.is_visual_blockwise()
    return M.get_mode() == ctrl_v
end

---Notify of stderr or stdout from a vim.system call obj, on the main loop.
---Intended for stderr. Uses stdout if stderr is empty, which may be the case if program doesn't utilise stderr at all.
---@param obj vim.SystemCompleted
---@param level integer|nil `vim.log.levels`, default ERROR
function M.schedule_notify(obj, level)
    local text = (obj.stderr or ""):gsub("\n$", "")
    if text == "" then
        text = (obj.stdout or ""):gsub("\n$", "")
    end
    vim.schedule(function() -- notify when we are ready
        vim.notify(text, level or vim.log.levels.ERROR) -- vim.notify instead of print to see multiple lines
    end)
end

---Notify an error for a command that exited non-zero, naming the command, its
---exit code and its stderr. Does nothing on exit 0. Scheduled, so it is safe in
---a fast context such as a `vim.system` exit callback.
---@param cmd string[] the command that ran
---@param obj vim.SystemCompleted its result
function M.notify_failure(cmd, obj)
    if obj.code == 0 then return end
    local msg = ("%s failed (exit %d):\n%s"):format(table.concat(cmd, " "), obj.code, vim.trim(obj.stderr or ""))
    vim.schedule(function() vim.notify(msg, vim.log.levels.ERROR) end)
end

---Byte column of the start of the keyword under the cursor.
---Uses the `\<` start-of-word regex atom rather than text-matching <cword>,
---so it is correct when the word occurs more than once on the line.
---@return integer col 1-indexed byte column (cursor's own column if not on a word)
function M.cword_start_col()
    return vim.fn.searchpos([[\<]], "bcn")[2]
end

function M.is_mac()
    return vim.uv.os_uname().sysname == "Darwin"
end

---Build a shell command running a script from the directory it lives in.
---Escaped for `:!`, so paths containing spaces, `;`, `(`, `!` or `%` are safe.
---@param file string Absolute path to the script
---@param interpreter string|nil Command to run it with; nil executes the script
---  itself, relying on its shebang
---@return string Shell command
function M.script_command(file, interpreter)
    local name = vim.fs.basename(file)
    local run = interpreter and interpreter .. ' ' .. vim.fn.shellescape(name, true)
        or vim.fn.shellescape('./' .. name, true)
    return 'cd ' .. vim.fn.shellescape(vim.fs.dirname(file), true) .. ' && ' .. run
end

return M
