---@diagnostic disable: undefined-global
-- Tests for completion/mlr/blink_mlr_columns.lua

local mod = require("completion.mlr.blink_mlr_columns")

--- Create a buffer with the given lines, set filetype to zsh, and wait for
--- tree-sitter to parse. Returns the buffer handle.
-- Match the real config: zsh files use compound ft "sh.zsh", and
-- lua/autocmds/treesitter.lua registers "zsh" parser for "sh.zsh".
vim.treesitter.language.register("zsh", "sh.zsh")

local function create_zsh_buffer(lines)
    local buf = vim.api.nvim_create_buf(true, false)
    vim.api.nvim_set_current_buf(buf)
    vim.api.nvim_buf_set_lines(buf, 0, -1, false, lines)
    vim.bo[buf].filetype = "sh.zsh"
    vim.treesitter.get_parser(buf):parse()
    return buf
end

--- Run the provider at the given (1-indexed row, 0-indexed col) and return
--- the list of completion item labels.
local function get_labels(row, col)
    vim.api.nvim_win_set_cursor(0, { row, col })
    local labels = {}
    local provider = mod.new()
    provider:get_completions(nil, function(result)
        for _, item in ipairs(result.items) do
            labels[#labels + 1] = item.label
        end
    end)
    return labels
end

--- End-of-line col for a given buffer line (0-indexed row).
local function eol(row)
    return #vim.api.nvim_buf_get_lines(0, row, row + 1, false)[1]
end

local test_tsv = vim.fn.tempname() .. ".tsv"

local function setup_data()
    local fh = io.open(test_tsv, "w")
    fh:write("name\tscore\tcategory\n")
    fh:write("foo\t0.8\tA\n")
    fh:close()
end

setup_data()

describe("mlr column completion", function()
    local buf

    after_each(function()
        if buf then
            pcall(vim.api.nvim_buf_delete, buf, { force = true })
            buf = nil
        end
    end)

    describe("chain detection", function()
        it("returns empty outside mlr command", function()
            buf = create_zsh_buffer({ "echo hello" })
            local labels = get_labels(1, 5)
            assert.are.same({}, labels)
        end)

        it("finds columns with --from", function()
            local line = "mlr -t --from " .. test_tsv .. " cut -f x"
            buf = create_zsh_buffer({ line })
            local labels = get_labels(1, #line - 1)
            assert.is_true(vim.tbl_contains(labels, "name"))
            assert.is_true(vim.tbl_contains(labels, "score"))
            assert.is_true(vim.tbl_contains(labels, "category"))
        end)

        it("finds columns with trailing filename", function()
            local line = "mlr -t filter '$score > 0' " .. test_tsv
            buf = create_zsh_buffer({ line })
            -- cursor on the path at end
            local labels = get_labels(1, #line - 1)
            assert.is_true(vim.tbl_contains(labels, "name"))
        end)

        it("follows + chain backwards", function()
            buf = create_zsh_buffer({
                "mlr -t --from " .. test_tsv .. " filter '$score > 0' +\\",
                "    cut -f x",
            })
            local labels = get_labels(2, 11)
            assert.is_true(vim.tbl_contains(labels, "name"))
            assert.is_true(vim.tbl_contains(labels, "score"))
        end)
    end)

    describe("whitespace cursor", function()
        it("works with cursor in trailing whitespace", function()
            buf = create_zsh_buffer({
                "mlr -t --from " .. test_tsv .. " cut -f ",
            })
            local labels = get_labels(1, eol(0))
            assert.is_true(vim.tbl_contains(labels, "name"))
        end)

        it("works with cursor before pipe", function()
            buf = create_zsh_buffer({
                "mlr -t --from " .. test_tsv .. " cut -f  | less",
            })
            local line = vim.api.nvim_buf_get_lines(0, 0, 1, false)[1]
            local col = string.find(line, "| less") - 2
            local labels = get_labels(1, col)
            assert.is_true(vim.tbl_contains(labels, "name"))
        end)
    end)

    describe("join -f path", function()
        it("collects path from join -f", function()
            local line = "mlr -t --from " .. test_tsv .. " join -f " .. test_tsv .. " -j name x"
            buf = create_zsh_buffer({ line })
            local labels = get_labels(1, #line - 1)
            assert.is_true(vim.tbl_contains(labels, "name"))
        end)
    end)

    describe("caching", function()
        it("returns same results on second call (cache hit)", function()
            local line = "mlr -t --from " .. test_tsv .. " cut -f x"
            buf = create_zsh_buffer({ line })
            local col = #line - 1
            local labels1 = get_labels(1, col)
            local labels2 = get_labels(1, col)
            assert.are.same(labels1, labels2)
        end)
    end)
end)
