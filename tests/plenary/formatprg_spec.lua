---@diagnostic disable: undefined-global

describe("utils.formatprg.formatexpr", function()
    local notify = vim.notify
    local notes

    before_each(function()
        notes = {}
        vim.notify = function(msg, level) notes[#notes + 1] = { msg = msg, level = level } end
        vim.cmd("enew!")
        vim.bo.formatexpr = "v:lua.require'utils.formatprg'.formatexpr()"
        vim.api.nvim_buf_set_lines(0, 0, -1, true, { "a", "b", "c" })
    end)

    after_each(function() vim.notify = notify end)

    local function lines() return vim.api.nvim_buf_get_lines(0, 0, -1, true) end

    it("filters only the operated lines", function()
        vim.bo.formatprg = "tr a-z A-Z"
        vim.cmd("normal! 2Ggqq")
        assert.are.same({ "a", "B", "c" }, lines())
        assert.are.same({}, notes)
    end)

    it("notifies stderr at WARN instead of inserting it", function()
        vim.bo.formatprg = "tr a-z A-Z; echo careful >&2"
        vim.cmd("normal! gggqG")
        assert.are.same({ "A", "B", "C" }, lines())
        assert.are.same({ { msg = "careful", level = vim.log.levels.WARN } }, notes)
    end)

    it("leaves the lines on a failing first pipeline stage", function()
        vim.bo.formatprg = "echo bad >&2; false | cat"
        vim.cmd("normal! gggqG")
        vim.wait(50, function() return #notes > 0 end)
        assert.are.same({ "a", "b", "c" }, lines())
        assert.are.same(vim.log.levels.ERROR, notes[1].level)
        assert.truthy(notes[1].msg:find("bad", 1, true))
    end)
end)
