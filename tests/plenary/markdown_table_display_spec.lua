---@diagnostic disable: undefined-global
-- Screen tests for lua/markdown_table_display.lua: a child nvim read with
-- `screenstring`.

local CONFIG = vim.uv.cwd()

local LINES = {
    "text",
    "| correlation with K50 | trypsin    | chymotrypsin |",
    "| -------------------- | ---------- | ------------ |",
    "| K/R count            | **-0.215** | -0.064       |",
    "| F/Y/W/L count        | -0.168     | **-0.414**   |",
}
--- Screen columns of the pipes of every table row, with the pads drawn.
local PIPES = { 1, 24, 39, 56 }
local BOLD_LINE = 4

local chan, child

--- Run `code` in the child.
--- @param code string Lua chunk
--- @param ... any arguments, as `...` in the chunk
--- @return any
local function lua(code, ...)
    return vim.rpcrequest(chan, "nvim_exec_lua", code, { ... })
end

--- Let the child process pending events, then redraw it.
local function flush()
    vim.wait(50)
    lua("vim.cmd.redraw()")
end

--- Screen text of the table rows of `win`.
--- @param win integer child window
--- @return string[]
local function table_rows(win)
    return lua([[
        local win = ...
        local top = vim.api.nvim_win_get_position(win)[1]
        local rows = {}
        for line = 2, 5 do
            local chars = {}
            for col = 1, vim.o.columns do
                chars[col] = vim.fn.screenstring(top + line, col)
            end
            rows[#rows + 1] = table.concat(chars)
        end
        return rows
    ]], win)
end

--- Columns of the `|` characters of a screen row.
--- @param row string
--- @return integer[]
local function pipes(row)
    local cols = {}
    for col in row:gmatch("()|") do
        cols[#cols + 1] = col
    end
    return cols
end

--- Assert every table row of `win` has its pipes at `PIPES`.
--- @param win integer
--- @return string[] rows
local function assert_aligned(win)
    local rows = table_rows(win)
    for i, row in ipairs(rows) do
        assert.are.same(PIPES, pipes(row), ("row %d: %s"):format(i, row))
    end
    return rows
end

describe("markdown_table_display", function()
    before_each(function()
        local sock = vim.fn.tempname()
        -- A headless nvim redraws into its screen grid at `:redraw`.
        child = vim.system({ "nvim", "--clean", "--listen", sock, "--headless", "--cmd", "set lines=24 columns=80" })
        assert.is_true(vim.wait(5000, function() return vim.uv.fs_stat(sock) ~= nil end, 10))
        chan = vim.fn.sockconnect("pipe", sock, { rpc = true })
        lua([[
            local config = ...
            vim.opt.rtp:prepend(config)
            vim.opt.rtp:append(config .. "/after")
            require("dev_plugins").add(config)
            vim.cmd.runtime("plugin/treesitter.lua")
            require("markdown_table_display").setup()
        ]], CONFIG)
        lua([[
            vim.api.nvim_buf_set_lines(0, 0, -1, false, ...)
            vim.bo.filetype = "markdown"
            vim.treesitter.start(0, "markdown")
            vim.wo.concealcursor = ""
        ]], LINES)
        flush()
    end)

    after_each(function()
        vim.fn.chanclose(chan)
        child:kill(9)
        child:wait()
    end)

    it("lines up the table at conceallevel 0 to 3", function()
        local win = lua("return vim.api.nvim_get_current_win()")
        for level = 0, 3 do
            lua("vim.wo.conceallevel = ...", level)
            flush()
            local rows = assert_aligned(win)
            assert.are.equal(level == 0, rows[3]:find("**", 1, true) ~= nil, "level " .. level)
        end
    end)

    it("keeps the columns while the cursor moves over a bold row", function()
        local win = lua("return vim.api.nvim_get_current_win()")
        lua("vim.wo.conceallevel = 2")
        flush()
        vim.rpcrequest(chan, "nvim_input", BOLD_LINE .. "G")
        flush()
        assert.is_truthy(assert_aligned(win)[3]:find("**-0.215**", 1, true))
        vim.rpcrequest(chan, "nvim_input", "j")
        flush()
        assert.is_truthy(assert_aligned(win)[4]:find("**-0.414**", 1, true))
        vim.rpcrequest(chan, "nvim_input", "gg")
        flush()
        assert.is_falsy(assert_aligned(win)[3]:find("**", 1, true))
    end)

    it("lines up two windows with different conceallevels", function()
        local top, bottom = unpack(lua([[
            vim.cmd.split()
            local top = vim.api.nvim_get_current_win()
            vim.wo[top].conceallevel = 0
            vim.cmd.wincmd("j")
            local bottom = vim.api.nvim_get_current_win()
            vim.wo[bottom].conceallevel = 3
            return { top, bottom }
        ]]))
        flush()
        assert.is_truthy(assert_aligned(top)[3]:find("**", 1, true))
        assert.is_falsy(assert_aligned(bottom)[3]:find("**", 1, true))
    end)

    it("leaves floats unpadded", function()
        local float_buf = lua([[
            local buf = vim.api.nvim_create_buf(false, true)
            vim.api.nvim_buf_set_lines(buf, 0, -1, false, ...)
            local win = vim.api.nvim_open_win(buf, true, { relative = "editor", row = 1, col = 1, width = 60, height = 6 })
            vim.bo[buf].filetype = "markdown"
            vim.treesitter.start(buf, "markdown")
            vim.wo[win].conceallevel = 2
            return buf
        ]], LINES)
        flush()
        assert.are.equal(0, #lua(
            'return vim.api.nvim_buf_get_extmarks(..., -1, 0, -1, { type = "virt_text" })', float_buf))
    end)

    it("moves the unconcealed cursor row with the current window", function()
        local top, bottom = unpack(lua([[
            vim.wo.conceallevel = 2
            vim.cmd.split()
            local wins = { vim.api.nvim_get_current_win() }
            vim.api.nvim_win_set_cursor(0, { ..., 0 })
            vim.cmd.wincmd("j")
            wins[2] = vim.api.nvim_get_current_win()
            vim.api.nvim_win_set_cursor(0, { ..., 0 })
            return wins
        ]], BOLD_LINE))
        flush()
        assert.is_falsy(assert_aligned(top)[3]:find("**", 1, true))
        assert.is_truthy(assert_aligned(bottom)[3]:find("**", 1, true))
        lua('vim.cmd.wincmd("k")')
        flush()
        assert.is_truthy(assert_aligned(top)[3]:find("**", 1, true))
        assert.is_falsy(assert_aligned(bottom)[3]:find("**", 1, true))
    end)
end)
