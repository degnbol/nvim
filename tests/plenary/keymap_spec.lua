---@diagnostic disable: undefined-global
local map = require "utils/keymap"

describe("keymap mode wrappers", function()
    after_each(function() pcall(vim.keymap.del, "n", "<F13>") end)

    it("set desc without writing into the caller's opts", function()
        local opts = { silent = true }
        map.n("<F13>", "<Nop>", "test desc", opts)
        assert.are.same({ silent = true }, opts)
        assert.are.equal("test desc", vim.fn.maparg("<F13>", "n", false, true).desc)
    end)

    it("buf sets a buffer-local map without writing into the caller's opts", function()
        local opts = {}
        map.buf("n", "<F13>", "<Nop>", "test desc", opts)
        assert.are.same({}, opts)
        assert.are.equal(1, vim.fn.maparg("<F13>", "n", false, true).buffer)
        vim.keymap.del("n", "<F13>", { buffer = true })
    end)
end)
