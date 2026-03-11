---@diagnostic disable: undefined-global
-- Tests for completion/mlr/blink_mlr_columns.lua

local mod = require("completion.mlr.blink_mlr")

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
--- the list of completion items (label + kind).
--- Set insert_mode = true to allow cursor one past line end (simulates insert mode).
local function get_items(row, col, insert_mode)
    local old_ve = vim.o.virtualedit
    if insert_mode then vim.o.virtualedit = "onemore" end
    vim.api.nvim_win_set_cursor(0, { row, col })
    local result_items = {}
    local provider = mod.new()
    provider:get_completions(nil, function(result)
        for _, item in ipairs(result.items) do
            result_items[#result_items + 1] = item
        end
    end)
    vim.o.virtualedit = old_ve
    return result_items
end

local function get_labels(row, col, insert_mode)
    local items = get_items(row, col, insert_mode)
    local labels = {}
    for _, item in ipairs(items) do
        labels[#labels + 1] = item.label
    end
    return labels
end

local function get_kinds(row, col)
    local items = get_items(row, col)
    local kinds = {}
    for _, item in ipairs(items) do
        kinds[item.label] = item.kind
    end
    return kinds
end

local Kind = vim.lsp.protocol.CompletionItemKind

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

describe("mlr completion", function()
    local buf

    after_each(function()
        if buf then
            pcall(vim.api.nvim_buf_delete, buf, { force = true })
            buf = nil
        end
    end)

    describe("verb completion", function()
        it("offers verbs after mlr flags", function()
            buf = create_zsh_buffer({ "mlr -t --from " .. test_tsv .. " " })
            local labels = get_labels(1, eol(0))
            assert.is_true(vim.tbl_contains(labels, "cut"))
            assert.is_true(vim.tbl_contains(labels, "filter"))
            assert.is_true(vim.tbl_contains(labels, "sort"))
            assert.is_true(vim.tbl_contains(labels, "head"))
        end)

        it("returns Keyword kind for verbs", function()
            buf = create_zsh_buffer({ "mlr -t --from " .. test_tsv .. " " })
            local kinds = get_kinds(1, eol(0))
            assert.are.equal(Kind.Keyword, kinds["cut"])
        end)

        it("offers verbs after + chain", function()
            buf = create_zsh_buffer({
                "mlr -t --from " .. test_tsv .. " cut -f name +\\",
                "    ",
            })
            local labels = get_labels(2, eol(1))
            assert.is_true(vim.tbl_contains(labels, "sort"))
            assert.is_true(vim.tbl_contains(labels, "head"))
        end)

        it("offers verbs after + on same line", function()
            local line = "mlr -t --from " .. test_tsv .. " cut -f name + "
            buf = create_zsh_buffer({ line })
            local labels = get_labels(1, #line)
            assert.is_true(vim.tbl_contains(labels, "filter"))
        end)

        it("returns empty outside mlr command", function()
            buf = create_zsh_buffer({ "echo hello" })
            local labels = get_labels(1, 5)
            assert.are.same({}, labels)
        end)
    end)

    describe("flag completion", function()
        it("offers flags after verb name", function()
            local line = "mlr -t --from " .. test_tsv .. " cut "
            buf = create_zsh_buffer({ line })
            local labels = get_labels(1, #line)
            assert.is_true(vim.tbl_contains(labels, "-f"))
            assert.is_true(vim.tbl_contains(labels, "-x"))
            assert.is_true(vim.tbl_contains(labels, "+"))
        end)

        it("returns Property kind for flags", function()
            local line = "mlr -t --from " .. test_tsv .. " cut "
            buf = create_zsh_buffer({ line })
            local kinds = get_kinds(1, #line)
            assert.are.equal(Kind.Property, kinds["-f"])
            assert.are.equal(Kind.Operator, kinds["+"])
        end)

        it("does not include -h/--help", function()
            local line = "mlr -t --from " .. test_tsv .. " cut "
            buf = create_zsh_buffer({ line })
            local labels = get_labels(1, #line)
            assert.is_false(vim.tbl_contains(labels, "-h"))
            assert.is_false(vim.tbl_contains(labels, "--help"))
        end)

        it("textEdit replaces typed dash prefix to avoid --f", function()
            local line = "mlr -t --from " .. test_tsv .. " cut -"
            buf = create_zsh_buffer({ line })
            local items = get_items(1, #line, true)
            local f_item
            for _, item in ipairs(items) do
                if item.label == "-f" then f_item = item; break end
            end
            assert.is_not_nil(f_item)
            assert.is_not_nil(f_item.textEdit)
            -- textEdit range should start at the '-', not after it
            local start_char = f_item.textEdit.range.start.character
            assert.are.equal(string.find(line, "%-$") - 1, start_char)
            assert.are.equal("-f", f_item.textEdit.newText)
        end)

        it("offers flags for sort verb", function()
            local line = "mlr -t --from " .. test_tsv .. " sort "
            buf = create_zsh_buffer({ line })
            local labels = get_labels(1, #line)
            assert.is_true(vim.tbl_contains(labels, "-f"))
            assert.is_true(vim.tbl_contains(labels, "-nr"))
            assert.is_true(vim.tbl_contains(labels, "-nf"))
        end)
    end)

    describe("field completion (context-aware)", function()
        it("offers columns after -f flag", function()
            local line = "mlr -t --from " .. test_tsv .. " cut -f "
            buf = create_zsh_buffer({ line })
            local labels = get_labels(1, #line)
            assert.is_true(vim.tbl_contains(labels, "name"))
            assert.is_true(vim.tbl_contains(labels, "score"))
            assert.is_true(vim.tbl_contains(labels, "category"))
        end)

        it("offers columns after -g flag", function()
            local line = "mlr -t --from " .. test_tsv .. " count -g "
            buf = create_zsh_buffer({ line })
            local labels = get_labels(1, #line)
            assert.is_true(vim.tbl_contains(labels, "name"))
        end)

        it("does not offer columns after non-field flag", function()
            -- -n takes a number, not field names
            local line = "mlr -t --from " .. test_tsv .. " head -n "
            buf = create_zsh_buffer({ line })
            local labels = get_labels(1, #line)
            -- Should get flags, not columns
            assert.is_false(vim.tbl_contains(labels, "name"))
        end)

        it("finds columns with trailing filename", function()
            local line = "mlr -t filter '$score > 0' " .. test_tsv
            buf = create_zsh_buffer({ line })
            -- cursor on the path at end — context is still "field" because
            -- filter is not a positional_field_verb, so it falls to "flag"
            -- but let's check that paths are still collected
            local labels = get_labels(1, #line - 1)
            -- filter doesn't have field flags, so we get flag items with +
            assert.is_true(vim.tbl_contains(labels, "+"))
        end)

        it("follows + chain backwards for column data", function()
            buf = create_zsh_buffer({
                "mlr -t --from " .. test_tsv .. " filter '$score > 0' +\\",
                "    cut -f ",
            })
            local labels = get_labels(2, eol(1))
            assert.is_true(vim.tbl_contains(labels, "name"))
            assert.is_true(vim.tbl_contains(labels, "score"))
        end)

        it("works with cursor in trailing whitespace after -f", function()
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

        it("collects path from join -f", function()
            local line = "mlr -t --from " .. test_tsv
                .. " join -f " .. test_tsv .. " -j "
            buf = create_zsh_buffer({ line })
            local labels = get_labels(1, #line)
            assert.is_true(vim.tbl_contains(labels, "name"))
        end)
    end)

    describe("positional field verbs", function()
        it("offers columns for group-by args", function()
            local line = "mlr -t --from " .. test_tsv .. " group-by "
            buf = create_zsh_buffer({ line })
            local labels = get_labels(1, #line)
            assert.is_true(vim.tbl_contains(labels, "name"))
            assert.is_true(vim.tbl_contains(labels, "score"))
        end)

        it("offers columns for label args", function()
            local line = "mlr -t --from " .. test_tsv .. " label "
            buf = create_zsh_buffer({ line })
            local labels = get_labels(1, #line)
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
            local line = "mlr -t --from " .. test_tsv .. " rename name,x"
            buf = create_zsh_buffer({ line })
            local labels = get_labels(1, #line - 1)
            assert.are.same({}, labels)
        end)

        it("shows completions at second from position (position 2)", function()
            local line = "mlr -t --from " .. test_tsv .. " rename name,new_name,x"
            buf = create_zsh_buffer({ line })
            local labels = get_labels(1, #line - 1)
            assert.is_true(vim.tbl_contains(labels, "name"))
        end)

        it("suppresses completions at second to position (position 3)", function()
            local line = "mlr -t --from " .. test_tsv .. " rename name,new_name,score,x"
            buf = create_zsh_buffer({ line })
            local labels = get_labels(1, #line - 1)
            assert.are.same({}, labels)
        end)

        it("does not suppress in verb after rename (single command node)", function()
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
            assert.is_true(vim.tbl_contains(labels, "Name"))
            assert.is_false(vim.tbl_contains(labels, "name"))
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
            local line = "mlr -t --from " .. test_tsv
                .. " rename ,Element + cut -f name"
            buf = create_zsh_buffer({ line })
            local comma_pos = line:find(",Element")
            local labels = get_labels(1, comma_pos - 1)
            assert.is_true(vim.tbl_contains(labels, "name"))
            assert.is_true(vim.tbl_contains(labels, "score"))
        end)
    end)

    describe("DSL string completion", function()
        it("offers columns inside filter DSL string", function()
            local line = "mlr -t --from " .. test_tsv .. " filter '$sc'"
            buf = create_zsh_buffer({ line })
            local quote_start = line:find("'%$sc")
            local labels = get_labels(1, quote_start + 3)
            assert.is_true(vim.tbl_contains(labels, "score"))
            assert.is_true(vim.tbl_contains(labels, "name"))
        end)

        it("offers columns inside put DSL string", function()
            local line = "mlr -t --from " .. test_tsv .. " put '$new = $'"
            buf = create_zsh_buffer({ line })
            local labels = get_labels(1, #line - 1)
            assert.is_true(vim.tbl_contains(labels, "score"))
        end)

        it("returns empty inside non-DSL string", function()
            local line = "mlr -t --from " .. test_tsv .. " grep 'foo'"
            buf = create_zsh_buffer({ line })
            local quote_start = line:find("'foo")
            local labels = get_labels(1, quote_start + 2)
            assert.are.same({}, labels)
        end)
    end)

    describe("partial verb completion", function()
        it("offers verbs when typing partial verb name (normal mode)", function()
            local line = "mlr -t --from " .. test_tsv .. " filt"
            buf = create_zsh_buffer({ line })
            local labels = get_labels(1, #line)
            assert.is_true(vim.tbl_contains(labels, "filter"))
            assert.is_true(vim.tbl_contains(labels, "cut"))
        end)

        it("offers verbs when typing partial verb after +", function()
            local line = "mlr -t --from " .. test_tsv .. " cut -f name + so"
            buf = create_zsh_buffer({ line })
            local labels = get_labels(1, #line)
            assert.is_true(vim.tbl_contains(labels, "sort"))
        end)

        -- Insert mode: cursor one past line end (at tok.ec boundary)
        it("offers verbs when typing partial verb name (insert mode)", function()
            local line = "mlr -t --from " .. test_tsv .. " filt"
            buf = create_zsh_buffer({ line })
            local labels = get_labels(1, #line, true)
            assert.is_true(vim.tbl_contains(labels, "filter"))
            assert.is_false(vim.tbl_contains(labels, "-f"))
        end)

        it("offers flags not verbs after complete verb (insert mode)", function()
            local line = "mlr -t --from " .. test_tsv .. " cut "
            buf = create_zsh_buffer({ line })
            local labels = get_labels(1, #line, true)
            assert.is_true(vim.tbl_contains(labels, "-f"))
            assert.is_false(vim.tbl_contains(labels, "filter"))
        end)

        it("offers columns after -f in insert mode", function()
            local line = "mlr -t --from " .. test_tsv .. " cut -f "
            buf = create_zsh_buffer({ line })
            local labels = get_labels(1, #line, true)
            assert.is_true(vim.tbl_contains(labels, "name"))
            assert.is_true(vim.tbl_contains(labels, "score"))
        end)
    end)

    describe("caching", function()
        it("returns same results on second call (cache hit)", function()
            local line = "mlr -t --from " .. test_tsv .. " cut -f "
            buf = create_zsh_buffer({ line })
            local col = #line
            local labels1 = get_labels(1, col)
            local labels2 = get_labels(1, col)
            assert.are.same(labels1, labels2)
        end)
    end)
end)
