---@diagnostic disable: undefined-global
local util = require "utils/init"

describe("utils.notify_failure", function()
    local notify = vim.notify
    local notes

    before_each(function()
        notes = {}
        vim.notify = function(msg, level) notes[#notes + 1] = { msg = msg, level = level } end
    end)

    after_each(function() vim.notify = notify end)

    it("does nothing on exit 0", function()
        util.notify_failure({ "true" }, { code = 0, signal = 0, stderr = "", stdout = "" })
        vim.wait(50)
        assert.are.same({}, notes)
    end)

    it("names the command, exit code and stderr at ERROR level", function()
        util.notify_failure({ "keyd", "reload" }, { code = 2, signal = 0, stderr = "denied\n", stdout = "" })
        vim.wait(50, function() return #notes > 0 end)
        assert.are.same({ { msg = "keyd reload failed (exit 2):\ndenied", level = vim.log.levels.ERROR } }, notes)
    end)

    it("works from a vim.system callback (fast context)", function()
        vim.system({ "sh", "-c", "echo oops >&2; exit 3" }, { text = true }, function(obj)
            util.notify_failure({ "sh" }, obj)
        end)
        vim.wait(2000, function() return #notes > 0 end)
        assert.are.same({ { msg = "sh failed (exit 3):\noops", level = vim.log.levels.ERROR } }, notes)
    end)
end)
