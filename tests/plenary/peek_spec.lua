---@diagnostic disable: undefined-global
-- Tests for the K peek float in plugin/keymaps.lua: a `path:lnum` under the
-- cursor opens the file scrolled to that line, with the buffer still holding
-- the file's own lines. open_floating_preview rewrites what it is handed --
-- trimempty for a plain syntax, _normalize_markdown for markdown -- so the
-- line-number alignment is what these cases pin.

describe("K peek", function()
    -- The runner starts nvim with --noplugin (plenary test_harness.lua:90), so
    -- the K mapping is not installed until this file is sourced.
    vim.cmd.runtime("plugin/keymaps.lua")
    -- The markdown float sets 'filetype', and ftplugin/markdown.lua requires
    -- plugins that the harness does not install.
    vim.cmd("filetype plugin off")

    local dir

    --- Write `lines` to a fixture named `name`, then leave the window on a
    --- buffer holding only `token` (PATH substituted by the fixture path) with
    --- the cursor on its first byte.
    local function fixture(name, lines, token)
        local path = dir .. "/" .. name
        assert.are.equal(0, vim.fn.writefile(lines, path))
        vim.cmd.enew()
        vim.api.nvim_buf_set_lines(0, 0, -1, false, { (token:gsub("PATH", function() return path end)) })
        vim.api.nvim_win_set_cursor(0, { 1, 0 })
        return path
    end

    --- The float's buffer and its topline, or nil when K opened no float.
    local function peek()
        vim.cmd("normal K")
        for _, win in ipairs(vim.api.nvim_list_wins()) do
            if vim.api.nvim_win_get_config(win).relative ~= "" then
                local topline = vim.api.nvim_win_call(win, function()
                    return vim.fn.winsaveview().topline
                end)
                return vim.api.nvim_win_get_buf(win), topline
            end
        end
    end

    --- 300 numbered lines, the first two blank -- trimempty drops leading
    --- blanks, so every line below is off by two unless they are restored.
    local function numbered()
        local lines = { "", "" }
        for i = 3, 300 do
            lines[i] = "line " .. i
        end
        return lines
    end

    before_each(function()
        dir = vim.fn.tempname()
        vim.fn.mkdir(dir, "p")
        vim.cmd("syntax enable") -- open_floating_preview only stylizes with it
    end)

    after_each(function()
        for _, win in ipairs(vim.api.nvim_list_wins()) do
            if vim.api.nvim_win_get_config(win).relative ~= "" then
                vim.api.nvim_win_close(win, true)
            end
        end
        vim.fs.rm(dir, { recursive = true, force = true })
    end)

    it("scrolls to the suffix line, keeping the file's own line numbers", function()
        fixture("a.lua", numbered(), "PATH:170")
        local buf, topline = peek()
        assert.are.equal(170, topline)
        assert.are.equal("line 170", vim.api.nvim_buf_get_lines(buf, 169, 170, false)[1])
    end)

    it("keeps the alignment through the markdown rewrite", function()
        -- markdown takes _normalize_markdown instead, which also collapses
        -- blank runs and rewrites thematic breaks.
        local lines = numbered()
        lines[4] = "" -- a blank run, collapsed by the markdown path
        lines[5] = "---"
        fixture("a.md", lines, "PATH:170")
        local buf, topline = peek()
        assert.are.equal(170, topline)
        assert.are.equal("line 170", vim.api.nvim_buf_get_lines(buf, 169, 170, false)[1])
    end)

    it("opens at the top without a suffix", function()
        fixture("b.lua", numbered(), "PATH")
        local _, topline = peek()
        assert.are.equal(1, topline)
    end)

    it("opens at the top when the suffix is past the end of the file", function()
        fixture("c.lua", numbered(), "PATH:9999")
        local _, topline = peek()
        assert.are.equal(1, topline)
    end)

    it("reports an empty file instead of opening a float", function()
        fixture("d.lua", {}, "PATH")
        local notified
        local notify = vim.notify
        ---@diagnostic disable-next-line: duplicate-set-field
        vim.notify = function(msg) notified = msg end
        local buf = peek()
        vim.notify = notify
        assert.is_nil(buf)
        assert.is_truthy(notified and notified:match("empty"))
    end)
end)
