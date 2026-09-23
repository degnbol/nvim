---@diagnostic disable: undefined-global
local async = require("utils/async")
local vim_async = require("vim._async")

describe("async.system", function()
    it("returns the result on the main loop", function()
        local result, fast = vim_async.run(function()
            local result = async.system({ "echo", "hi" }, { text = true })
            return result, vim.in_fast_event()
        end):wait()
        assert.are.equal("hi\n", result.stdout)
        assert.is_false(fast)
    end)
end)

describe("async.shared", function()
    it("gives concurrent calls with one key one run", function()
        local n_runs = 0
        local shared = async.shared(function(key)
            n_runs = n_runs + 1
            vim_async.await(1, vim.schedule)
            return { key = key }
        end)
        local tasks = vim.tbl_map(function(key)
            return vim_async.run(function() return shared(key) end)
        end, { "a", "a", "b" })
        local results = vim.tbl_map(function(task) return task:wait() end, tasks)
        assert.are.equal(2, n_runs)
        assert.are.equal(results[1], results[2])
        assert.are.same({ key = "b" }, results[3])
    end)

    it("raises the run's error in every waiting call", function()
        local shared = async.shared(function()
            vim_async.await(1, vim.schedule)
            error("boom", 0)
        end)
        local tasks = { vim_async.run(function() return shared("k") end), vim_async.run(function() return shared("k") end) }
        for _, task in ipairs(tasks) do
            local ok, err = pcall(task.wait, task)
            assert.is_false(ok)
            assert.truthy(err:find("boom$"))
        end
    end)
end)
