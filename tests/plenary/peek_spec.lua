---@diagnostic disable: undefined-global
-- Tests for the K peek float in plugin/keymaps.lua: a `path:lnum` under the
-- cursor opens the file at that line, which is the float's first line.
-- open_floating_preview rewrites what it is handed — trimempty for a plain
-- syntax, _normalize_markdown for markdown — and either would shift that line
-- away from the top, so the float's leading lines are what these cases pin.

describe("K peek", function()
    -- The runner starts nvim with --noplugin (plenary test_harness.lua:90), so
    -- the K mapping is not installed until this file is sourced.
    vim.cmd.runtime("plugin/keymaps.lua")
    -- The markdown float sets 'filetype', and ftplugin/markdown.lua requires
    -- plugins that the harness does not install.
    vim.cmd("filetype plugin off")

    local dir

    --- Write `lines` to a fixture named `name` — gzipped when `name` ends in
    --- .gz — then leave the window on a buffer holding only `token` (PATH
    --- substituted by the fixture path) with the cursor on its first byte.
    local function fixture(name, lines, token)
        local path = dir .. "/" .. name
        local plain = (path:gsub("%.gz$", ""))
        assert.are.equal(0, vim.fn.writefile(lines, plain))
        if plain ~= path then
            assert.are.equal(0, vim.system({ "gzip", "-f", plain }):wait().code)
        end
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

    --- 300 numbered lines, with `blank` and the line after it emptied —
    --- trimempty drops a leading blank run and _normalize_markdown collapses
    --- it, so a float opening on those two shows the wrong lines unless the
    --- lines as read are restored.
    local function numbered(blank)
        local lines = {}
        for i = 1, 300 do
            lines[i] = "line " .. i
        end
        lines[blank], lines[blank + 1] = "", ""
        return lines
    end

    before_each(function()
        dir = vim.fn.tempname()
        vim.fn.mkdir(dir, "p")
        vim.cmd("syntax enable") -- open_floating_preview only stylizes with it
    end)

    local max_bytes = require("utils.file").MAX_BYTES

    after_each(function()
        require("utils.file").MAX_BYTES = max_bytes
        for _, win in ipairs(vim.api.nvim_list_wins()) do
            if vim.api.nvim_win_get_config(win).relative ~= "" then
                vim.api.nvim_win_close(win, true)
            end
        end
        vim.fs.rm(dir, { recursive = true, force = true })
    end)

    local function head(buf)
        return vim.api.nvim_buf_get_lines(buf, 0, 3, false)
    end

    --- The messages `fn` notifies, in order.
    local function notifications(fn)
        local msgs = {}
        local notify = vim.notify
        ---@diagnostic disable-next-line: duplicate-set-field
        vim.notify = function(msg) msgs[#msgs + 1] = msg end
        local ok, err = pcall(fn)
        vim.notify = notify
        assert(ok, err)
        return msgs
    end

    --- 300 lines past the 64 kB the read cap is weighed against.
    local function padded()
        local lines = {}
        for i = 1, 300 do
            lines[i] = "line " .. i .. string.rep("-", 400)
        end
        return lines
    end

    it("opens the float on the suffix line", function()
        fixture("a.lua", numbered(170), "PATH:170")
        local buf, topline = peek()
        assert.are.equal(1, topline)
        assert.are.same({ "", "", "line 172" }, head(buf))
    end)

    it("opens on the suffix line through the markdown rewrite", function()
        -- markdown takes _normalize_markdown instead, which also collapses
        -- blank runs and rewrites thematic breaks.
        local lines = numbered(170)
        lines[172] = "---"
        fixture("a.md", lines, "PATH:170")
        local buf = peek()
        assert.are.same({ "", "", "---" }, head(buf))
    end)

    it("opens at the top without a suffix", function()
        fixture("b.lua", numbered(1), "PATH")
        local buf, topline = peek()
        assert.are.equal(1, topline)
        assert.are.same({ "", "", "line 3" }, head(buf))
    end)

    it("opens at the top when the suffix is past the end of the file", function()
        fixture("c.lua", numbered(1), "PATH:9999")
        local buf, topline
        local msgs = notifications(function() buf, topline = peek() end)
        assert.are.equal(1, topline)
        assert.are.same({ "", "", "line 3" }, head(buf))
        assert.are.same({ "c.lua: line 9999 is past the end of the file or the read cap" }, msgs)
    end)

    it("opens at the top when the read cap stops short of the suffix", function()
        fixture("f.lua", padded(), "PATH:290")
        require("utils.file").MAX_BYTES = 1 -- the first chunk read overshoots this
        local buf
        local msgs = notifications(function() buf = peek() end)
        assert.is_truthy(vim.api.nvim_buf_get_lines(buf, 0, 1, false)[1]:match("^line 1%-+$"))
        assert.are.same({ "f.lua: line 290 is past the end of the file or the read cap" }, msgs)
    end)

    it("decompresses a gzipped file", function()
        fixture("d.tsv.gz", numbered(170), "PATH:170")
        local buf = peek()
        assert.are.same({ "", "", "line 172" }, head(buf))
        assert.are.equal("tsv", vim.bo[buf].syntax) -- the .gz is not the content
    end)

    it("reports an empty file instead of opening a float", function()
        fixture("d.lua", {}, "PATH")
        local buf
        local msgs = notifications(function() buf = peek() end)
        assert.is_nil(buf)
        assert.are.same({ "d.lua is empty" }, msgs)
    end)

    it("reports an empty file once, suffix or not", function()
        fixture("d.lua", {}, "PATH:5")
        assert.are.same({ "d.lua is empty" }, notifications(peek))
    end)

    it("enters the float on a second K, and leaves it on a third", function()
        fixture("g.lua", numbered(170), "PATH:170")
        local source = vim.api.nvim_get_current_win()
        vim.cmd("normal K")
        assert.are.equal(source, vim.api.nvim_get_current_win())
        vim.cmd("normal K")
        local float = vim.api.nvim_get_current_win()
        assert.are_not.equal(source, float)
        assert.are_not.equal("", vim.api.nvim_win_get_config(float).relative)
        assert.are.same({ "", "", "line 172" }, head(vim.api.nvim_win_get_buf(float)))
        vim.cmd("normal K")
        assert.are.equal(source, vim.api.nvim_get_current_win())
    end)
end)
