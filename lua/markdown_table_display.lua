--- Line up markdown pipe tables on screen under conceal, with inline virtual
--- spaces from `utils.markdown_table.pads`. The file is not changed.
---
--- Pads depend on the window's 'conceallevel' and on its unconcealed cursor
--- row, so each window has its own namespace, scoped to it with
--- `nvim__ns_set`. The marks are persistent: ephemeral inline marks set from a
--- decoration provider do not render.
---
--- Disable in a buffer with `vim.b.markdown_table_display = false`.
local M = {}

local md_table = require "utils.markdown_table"
local api = vim.api

local OPTS = { max_width = md_table.MAX_WIDTH, align_cells = true }

--- Letter of 'concealcursor' per first character of `nvim_get_mode().mode`.
local MODE_LETTER = {
    n = "n",
    v = "v",
    V = "v",
    ["\22"] = "v",
    s = "v",
    S = "v",
    ["\19"] = "v",
    i = "i",
    R = "i",
    c = "c",
}

local win_ns = {} --- @type table<integer, integer> window → namespace
local free_ns = {} --- @type integer[] namespaces of closed windows
local win_buf = {} --- @type table<integer, integer> window → buffer last rendered in it
local win_cursor = {} --- @type table<integer, integer|false> window → unconcealed row last rendered
local pending = {} --- @type table<integer, true> buffers with a scheduled render
local attached = {} --- @type table<integer, true> buffers with a change listener
--- Per buffer, `{ layout, conceals }` per table, keyed by its row texts. A
--- change render moves the entries it uses from `previous` to `current`
--- and drops the rest.
--- @type table<integer, { current: table<string, table>, previous: table<string, table> }>
local cache = {}

--- Whether the display renders in `buf`.
--- @param buf integer
--- @return boolean
local function enabled(buf)
    return vim.tbl_contains(md_table.FILETYPES, vim.bo[buf].filetype)
        and vim.b[buf].markdown_table_display ~= false
end

--- The namespace of `win`, taken from the free ones or created.
--- @param win integer
--- @return integer
local function namespace(win)
    local ns = win_ns[win]
    if not ns then
        ns = table.remove(free_ns) or api.nvim_create_namespace("")
        api.nvim__ns_set(ns, { wins = { win } })
        win_ns[win] = ns
    end
    return ns
end

--- Buffer row drawn without conceal in `win`: the cursor row of the current
--- window, unless 'concealcursor' has the letter of the current mode.
--- @param win integer
--- @return integer|nil row 0-indexed
local function unconcealed_row(win)
    if win ~= api.nvim_get_current_win() then
        return nil
    end
    local letter = MODE_LETTER[api.nvim_get_mode().mode:sub(1, 1)]
    if not letter or vim.wo[win].concealcursor:find(letter, 1, true) then
        return nil
    end
    return api.nvim_win_get_cursor(win)[1] - 1
end

--- Layout and conceals of `tbl`, from the cache when its text is unchanged.
--- Not cached without a treesitter highlighter, which may start later.
--- @param buf integer
--- @param tbl utils.markdown_table.Table
--- @return { layout: utils.markdown_table.Layout, conceals: table<integer, agentic.utils.Conceal.Range[]> }
local function measured(buf, tbl)
    cache[buf] = cache[buf] or { current = {}, previous = {} }
    local entries = cache[buf]
    local lines = {}
    for i, row in ipairs(tbl.rows) do
        lines[i] = row.line
    end
    local key = table.concat(lines, "\n")
    local entry = entries.current[key] or entries.previous[key]
    if not entry then
        local conceals = md_table.conceals(buf, tbl)
        entry = { conceals = conceals, layout = md_table.layout(tbl, conceals, OPTS) }
    end
    if vim.treesitter.highlighter.active[buf] then
        entries.current[key] = entry
    end
    return entry
end

--- Last row of a `pipe_table` node.
--- @param node TSNode
--- @return integer row 0-indexed
local function last_row(node)
    local _, _, er, ec = node:range()
    return ec == 0 and er - 1 or er
end

--- Set the pads of the table `node` in `win` on `rows`, or on all its rows.
--- @param win integer
--- @param buf integer
--- @param node TSNode
--- @param rows table<integer, true>|nil 0-indexed
local function set_pads(win, buf, node, rows)
    local tbl = md_table.parse(buf, node)
    local entry = measured(buf, tbl)
    local ns = namespace(win)
    local level = vim.wo[win].conceallevel
    for _, pad in ipairs(md_table.pads(tbl, entry.layout, entry.conceals, level, unconcealed_row(win))) do
        if not rows or rows[pad.row] then
            api.nvim_buf_set_extmark(buf, ns, pad.row, pad.col, {
                virt_text = { { string.rep(" ", pad.width) } },
                virt_text_pos = "inline",
                hl_mode = "combine",
                undo_restore = false,
            })
        end
    end
end

--- Clear the pads of `buf` in every window namespace.
--- @param buf integer
local function clear_buf(buf)
    for _, ns in pairs(win_ns) do
        api.nvim_buf_clear_namespace(buf, ns, 0, -1)
    end
end

--- Render the pads of `win`: of the tables in its visible rows, or only on
--- `rows` of the tables there.
--- @param win integer
--- @param rows integer[]|nil 0-indexed
local function render(win, rows)
    -- The autocmd window is never drawn and gets no WinClosed. A float's
    -- width is usually fitted to its raw text, which pads would overflow.
    if vim.fn.win_gettype(win) == "autocmd" or api.nvim_win_get_config(win).relative ~= "" then
        return
    end
    local buf = api.nvim_win_get_buf(win)
    local old_buf = win_buf[win]
    if old_buf ~= buf then
        if old_buf and api.nvim_buf_is_valid(old_buf) then
            api.nvim_buf_clear_namespace(old_buf, win_ns[win], 0, -1)
        end
        rows = nil
    end
    if not enabled(buf) then
        if old_buf == buf then
            clear_buf(buf)
        end
        win_buf[win] = nil
        return
    end
    win_buf[win] = buf
    win_cursor[win] = unconcealed_row(win) or false
    local ns = namespace(win)

    if rows then
        local row_set, nodes = {}, {}
        for _, row in ipairs(rows) do
            row_set[row] = true
            api.nvim_buf_clear_namespace(buf, ns, row, row + 1)
            for _, node in ipairs(md_table.tables(buf, row, row)) do
                nodes[node:id()] = node
            end
        end
        for _, node in pairs(nodes) do
            set_pads(win, buf, node, row_set)
        end
        return
    end

    -- 'w$' is stale before the window's first redraw; a window height of
    -- rows below the top covers it.
    local top = vim.fn.line("w0", win) - 1
    local bottom = math.min(
        math.max(vim.fn.line("w$", win) - 1, top + api.nvim_win_get_height(win)),
        api.nvim_buf_line_count(buf) - 1
    )
    api.nvim_buf_clear_namespace(buf, ns, top, bottom + 1)
    for _, node in ipairs(md_table.tables(buf, top, bottom)) do
        api.nvim_buf_clear_namespace(buf, ns, node:start(), last_row(node) + 1)
        set_pads(win, buf, node)
    end
end

--- Render every window showing `buf`, and drop cache entries of table texts
--- none of them shows.
--- @param buf integer
local function render_buf(buf)
    if cache[buf] then
        cache[buf] = { current = {}, previous = cache[buf].current }
    end
    for _, win in ipairs(vim.fn.win_findbuf(buf)) do
        render(win)
    end
end

--- Render `buf` once, after the current batch of changes.
--- @param buf integer
local function schedule_render(buf)
    if pending[buf] then
        return
    end
    pending[buf] = true
    vim.schedule(function()
        pending[buf] = nil
        if api.nvim_buf_is_valid(buf) then
            render_buf(buf)
        end
    end)
end

--- Re-render the old and new unconcealed rows of `win` when they differ.
--- @param win integer
local function update_cursor(win)
    local buf = api.nvim_win_get_buf(win)
    if win_buf[win] ~= buf or pending[buf] then
        return
    end
    local old, new = win_cursor[win], unconcealed_row(win) or false
    if old == new then
        return
    end
    local rows = {}
    for _, row in ipairs({ old, new }) do
        if row then
            rows[#rows + 1] = row
        end
    end
    render(win, rows)
end

--- Render `buf` after each change to it. Once per buffer.
--- @param buf integer
local function attach(buf)
    if attached[buf] then
        return
    end
    attached[buf] = true
    api.nvim_buf_attach(buf, false, {
        on_lines = function()
            schedule_render(buf)
        end,
        on_reload = function()
            schedule_render(buf)
        end,
        on_detach = function()
            attached[buf] = nil
        end,
    })
end

--- Create the autocmds that keep the pads of every window current.
function M.setup()
    local group = api.nvim_create_augroup("markdown_table_display", { clear = true })
    api.nvim_create_autocmd({ "BufWinEnter", "FileType" }, {
        group = group,
        desc = "markdown_table_display: render the windows of a buffer",
        callback = function(ev)
            if enabled(ev.buf) then
                attach(ev.buf)
            end
            for _, win in ipairs(vim.fn.win_findbuf(ev.buf)) do
                render(win)
            end
        end,
    })
    api.nvim_create_autocmd("WinNew", {
        group = group,
        desc = "markdown_table_display: render a new window",
        callback = function()
            render(api.nvim_get_current_win())
        end,
    })
    api.nvim_create_autocmd("OptionSet", {
        group = group,
        pattern = { "conceallevel", "concealcursor" },
        desc = "markdown_table_display: render the window whose conceal changed",
        callback = function()
            -- Setting a window's option from another window makes it current
            -- only while OptionSet runs.
            local win = api.nvim_get_current_win()
            vim.schedule(function()
                if api.nvim_win_is_valid(win) then
                    render(win)
                end
            end)
        end,
    })
    api.nvim_create_autocmd("WinScrolled", {
        group = group,
        desc = "markdown_table_display: render scrolled windows",
        callback = function()
            for win in pairs(vim.v.event) do
                local id = tonumber(win)
                if id and api.nvim_win_is_valid(id) then
                    render(id)
                end
            end
        end,
    })
    api.nvim_create_autocmd({ "CursorMoved", "CursorMovedI", "ModeChanged" }, {
        group = group,
        desc = "markdown_table_display: re-render the rows the cursor conceal moved between",
        callback = function()
            update_cursor(api.nvim_get_current_win())
        end,
    })
    api.nvim_create_autocmd("WinEnter", {
        group = group,
        desc = "markdown_table_display: re-render the cursor rows of the windows that gained or lost focus",
        callback = function()
            for win in pairs(win_buf) do
                if api.nvim_win_is_valid(win) then
                    update_cursor(win)
                end
            end
        end,
    })
    api.nvim_create_autocmd("WinClosed", {
        group = group,
        desc = "markdown_table_display: free the namespace of a closed window",
        callback = function(ev)
            local win = tonumber(ev.match) --[[@as integer]]
            local ns, buf = win_ns[win], win_buf[win]
            if not ns then
                return
            end
            if buf and api.nvim_buf_is_valid(buf) then
                api.nvim_buf_clear_namespace(buf, ns, 0, -1)
            end
            win_ns[win], win_buf[win], win_cursor[win] = nil, nil, nil
            free_ns[#free_ns + 1] = ns
        end,
    })
    api.nvim_create_autocmd("BufWipeout", {
        group = group,
        desc = "markdown_table_display: drop the table cache of a buffer",
        callback = function(ev)
            cache[ev.buf] = nil
        end,
    })
end

return M
