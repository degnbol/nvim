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

    describe("rename verb", function()
        it("shows completions at from position (position 0)", function()
            local line = "mlr -t --from " .. test_tsv .. " rename "
            buf = create_zsh_buffer({ line })
            local labels = get_labels(1, #line)
            assert.is_true(vim.tbl_contains(labels, "name"))
            assert.is_true(vim.tbl_contains(labels, "score"))
        end)

        it("suppresses completions at to position (position 1)", function()
            -- Cursor on 'x' = one past comma = insert-mode "after comma" position
            local line = "mlr -t --from " .. test_tsv .. " rename name,x"
            buf = create_zsh_buffer({ line })
            local labels = get_labels(1, #line - 1)
            assert.are.same({}, labels)
        end)

        it("shows completions at second from position (position 2)", function()
            -- Cursor on 'x' after second comma = from position (even)
            local line = "mlr -t --from " .. test_tsv .. " rename name,new_name,x"
            buf = create_zsh_buffer({ line })
            local labels = get_labels(1, #line - 1)
            assert.is_true(vim.tbl_contains(labels, "name"))
        end)

        it("suppresses completions at second to position (position 3)", function()
            -- Cursor on 'x' after third comma = to position (odd)
            local line = "mlr -t --from " .. test_tsv .. " rename name,new_name,score,x"
            buf = create_zsh_buffer({ line })
            local labels = get_labels(1, #line - 1)
            assert.are.same({}, labels)
        end)

        it("does not suppress in verb after rename (single command node)", function()
            -- All on one line: mlr ... rename ... + cut -f ...
            local line = "mlr -t --from " .. test_tsv
                .. " rename name,Name + cut -f "
            buf = create_zsh_buffer({ line })
            local labels = get_labels(1, #line)
            assert.is_true(#labels > 0)
        end)

        it("shows completions for non-rename verbs in same chain", function()
            buf = create_zsh_buffer({
                "mlr -t --from " .. test_tsv .. " rename name,new_name +\\",
                "    cut -f ",
            })
            local labels = get_labels(2, eol(1))
            assert.is_true(#labels > 0)
        end)

        it("applies rename mapping to column names", function()
            buf = create_zsh_buffer({
                "mlr -t --from " .. test_tsv .. " rename name,Name +\\",
                "    cut -f ",
            })
            local labels = get_labels(2, eol(1))
            -- "name" should be renamed to "Name"
            assert.is_true(vim.tbl_contains(labels, "Name"))
            assert.is_false(vim.tbl_contains(labels, "name"))
            -- Other columns unaffected
            assert.is_true(vim.tbl_contains(labels, "score"))
        end)

        it("applies rename mapping in single-line chain", function()
            local line = "mlr -t --from " .. test_tsv
                .. " rename name,Name + cut -f "
            buf = create_zsh_buffer({ line })
            local labels = get_labels(1, #line)
            assert.is_true(vim.tbl_contains(labels, "Name"))
            assert.is_false(vim.tbl_contains(labels, "name"))
        end)

        it("shows completions right before comma in rename (from position)", function()
            -- Cursor ON the comma in ",Element" — in insert mode this means
            -- cursor is between the space and comma, text before cursor has 0 commas
            local line = "mlr -t --from " .. test_tsv
                .. " rename ,Element + cut -f name"
            buf = create_zsh_buffer({ line })
            -- find returns 1-indexed; -1 for 0-indexed = position of comma char
            local comma_pos = line:find(",Element")
            local labels = get_labels(1, comma_pos - 1)
            assert.is_true(vim.tbl_contains(labels, "name"))
            assert.is_true(vim.tbl_contains(labels, "score"))
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
