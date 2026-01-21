---@diagnostic disable: undefined-global
-- Tests for ftplugin/tsv.lua buffer-local state handling.
-- Tests vim.b.tsv_* isolation and vim.NIL conversion.

-- Helper to create a scratch buffer (no filetype to avoid syntax dependencies).
local function create_buffer()
    local buf = vim.api.nvim_create_buf(true, false)
    return buf
end

-- Helper to switch to a buffer.
local function switch_to_buffer(buf)
    vim.api.nvim_set_current_buf(buf)
end

describe("tsv buffer-local state", function()
    local buf1, buf2

    before_each(function()
        buf1 = create_buffer()
        buf2 = create_buffer()
    end)

    after_each(function()
        pcall(vim.api.nvim_buf_delete, buf1, { force = true })
        pcall(vim.api.nvim_buf_delete, buf2, { force = true })
    end)

    describe("vim.b isolation", function()
        it("tsv_widths are isolated between buffers", function()
            switch_to_buffer(buf1)
            vim.b.tsv_widths = { 10, 20, 30 }

            switch_to_buffer(buf2)
            vim.b.tsv_widths = { 5, 15 }

            -- Verify buf1 state is unchanged.
            switch_to_buffer(buf1)
            local widths1 = vim.b.tsv_widths
            assert.are.equal(3, #widths1)
            assert.are.equal(10, widths1[1])

            -- Verify buf2 state.
            switch_to_buffer(buf2)
            local widths2 = vim.b.tsv_widths
            assert.are.equal(2, #widths2)
            assert.are.equal(5, widths2[1])
        end)

        it("tsv_maxwidths are isolated between buffers", function()
            switch_to_buffer(buf1)
            vim.b.tsv_maxwidths = { [1] = 5, [2] = 8 }

            switch_to_buffer(buf2)
            -- buf2 should not have buf1's maxwidths.
            local maxwidths2 = vim.b.tsv_maxwidths
            assert.is_nil(maxwidths2)

            -- Set different value for buf2.
            vim.b.tsv_maxwidths = { [1] = 99 }

            -- Verify buf1 unchanged.
            switch_to_buffer(buf1)
            local maxwidths1 = vim.b.tsv_maxwidths
            assert.are.equal(5, maxwidths1[1])
            assert.are.equal(8, maxwidths1[2])
        end)

        it("tsv_hidden are isolated between buffers", function()
            switch_to_buffer(buf1)
            vim.b.tsv_hidden = { [1] = { [2] = "hidden_in_buf1" } }

            switch_to_buffer(buf2)
            vim.b.tsv_hidden = { [1] = { [2] = "hidden_in_buf2" } }

            -- Verify buf1 state unchanged.
            switch_to_buffer(buf1)
            local hidden1 = vim.b.tsv_hidden
            assert.are.equal("hidden_in_buf1", hidden1[1][2])

            -- Verify buf2 has its own state.
            switch_to_buffer(buf2)
            local hidden2 = vim.b.tsv_hidden
            assert.are.equal("hidden_in_buf2", hidden2[1][2])
        end)

        it("nil state in one buffer does not affect another", function()
            switch_to_buffer(buf1)
            vim.b.tsv_maxwidths = { [1] = 10 }
            vim.b.tsv_hidden = { [1] = { [1] = "test" } }

            switch_to_buffer(buf2)
            -- buf2 has nil state.
            assert.is_nil(vim.b.tsv_maxwidths)
            assert.is_nil(vim.b.tsv_hidden)

            -- buf1 still has its state.
            switch_to_buffer(buf1)
            assert.is_not_nil(vim.b.tsv_maxwidths)
            assert.is_not_nil(vim.b.tsv_hidden)
        end)
    end)
end)

-- vim.NIL handling tests (denilify function behaviour)
describe("vim.NIL in buffer variables", function()
    it("vim.b pads sparse tables with vim.NIL", function()
        local buf = vim.api.nvim_create_buf(true, false)
        vim.api.nvim_set_current_buf(buf)

        -- Store sparse table: only key 3 is set
        vim.b.test_sparse = { [3] = 10 }
        local retrieved = vim.b.test_sparse

        -- Vim pads with vim.NIL
        assert.are.equal(vim.NIL, retrieved[1])
        assert.are.equal(vim.NIL, retrieved[2])
        assert.are.equal(10, retrieved[3])

        vim.api.nvim_buf_delete(buf, { force = true })
    end)

    it("vim.NIL is truthy but not indexable", function()
        local buf = vim.api.nvim_create_buf(true, false)
        vim.api.nvim_set_current_buf(buf)

        vim.b.test_sparse = { [3] = { [5] = "text" } }
        local retrieved = vim.b.test_sparse

        -- vim.NIL is truthy
        assert.is_truthy(retrieved[1])
        -- vim.NIL is not nil
        assert.is_not_nil(retrieved[1])
        assert.are_not.equal(nil, retrieved[1])

        -- vim.NIL cannot be indexed (would error)
        assert.has_error(function()
            local _ = retrieved[1][1]
        end, "attempt to index a userdata value")

        vim.api.nvim_buf_delete(buf, { force = true })
    end)

    it("nested sparse tables also get vim.NIL padding", function()
        local buf = vim.api.nvim_create_buf(true, false)
        vim.api.nvim_set_current_buf(buf)

        vim.b.test_nested = { [2] = { [3] = "value" } }
        local retrieved = vim.b.test_nested

        -- Outer table padded
        assert.are.equal(vim.NIL, retrieved[1])
        -- Inner table also padded
        assert.are.equal(vim.NIL, retrieved[2][1])
        assert.are.equal(vim.NIL, retrieved[2][2])
        assert.are.equal("value", retrieved[2][3])

        vim.api.nvim_buf_delete(buf, { force = true })
    end)
end)
