---@diagnostic disable: undefined-global
-- Tests for ftplugin/tsv.lua buffer-local state handling (vim.b.tsv_* isolation
-- and vim.NIL conversion), and for the smarts injection into a chemical column:
-- queries/tsv/*.scm and lua/chem/tsv.lua.

-- plugin/* isn't reliably sourced before PlenaryBustedFile runs, and the
-- injection needs the `#unquote!` directive and the `#chem-column?` predicate.
vim.cmd.source(vim.fn.getcwd() .. "/plugin/treesitter.lua")
vim.cmd.source(vim.fn.getcwd() .. "/plugin/chem.lua")

-- A missing tsv parser would otherwise surface as "language could not be
-- determined" in every injection case below.
assert(vim.treesitter.language.add("tsv"))

local chem_tsv = require "chem/tsv"
-- Same call the module makes, so the same id: nvim_create_namespace is a lookup
-- for a name already taken.
local ns = vim.api.nvim_create_namespace("chem.elements")

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

--- Run a function on a scratch TSV buffer.
--- @param lines string[]
--- @param fn fun(buf: integer, parser: vim.treesitter.LanguageTree): any
--- @return any
local function with_tsv(lines, fn)
    local buf = vim.api.nvim_create_buf(false, true)
    vim.api.nvim_buf_set_lines(buf, 0, -1, false, lines)
    vim.bo[buf].filetype = "tsv"
    local result = fn(buf, assert(vim.treesitter.get_parser(buf)))
    vim.api.nvim_buf_delete(buf, { force = true })
    return result
end

--- The text of every structure a TSV buffer injects, in buffer order.
--- @param buf integer
--- @param parser vim.treesitter.LanguageTree
--- @return string[]
local function structures(buf, parser)
    parser:parse(true)
    local regions = {}
    parser:for_each_tree(function(tree, ltree)
        if ltree:lang() ~= "smarts" then return end
        local row, col = tree:root():range()
        regions[#regions + 1] =
            { row, col, vim.treesitter.get_node_text(tree:root(), buf) }
    end)
    table.sort(regions, function(a, b)
        return a[1] < b[1] or (a[1] == b[1] and a[2] < b[2])
    end)
    return vim.tbl_map(function(region) return region[3] end, regions)
end

--- The structures the cells of a TSV buffer inject.
--- @param lines string[]
--- @param prepare fun(buf: integer)|nil buffer state to set before parsing
--- @return string[]
local function injected(lines, prepare)
    return with_tsv(lines, function(buf, parser)
        if prepare then prepare(buf) end
        return structures(buf, parser)
    end)
end

--- The element groups a buffer's marks name, sorted and deduplicated.
--- @param buf integer
--- @return string[]
local function element_groups(buf)
    local names = {}
    for _, mark in ipairs(
        vim.api.nvim_buf_get_extmarks(buf, ns, 0, -1, { details = true })) do
        names[assert(mark[4]).hl_group] = true
    end
    local sorted = vim.tbl_keys(names)
    table.sort(sorted)
    return sorted
end

describe("chemical column injection", function()
    -- The predicate reads the comment character off 'commentstring', which is
    -- the ftplugin's to set.
    it("finds the comment character on a tsv buffer", function()
        assert.are.equal("#%s",
            with_tsv({ "a\tb" }, function(buf) return vim.bo[buf].commentstring end))
    end)

    it("injects the column the header names, and no other", function()
        assert.are.same({ "CCO" }, injected { "smiles\tname", "CCO\tCCN" })
    end)

    it("matches a header name whatever its case", function()
        assert.are.same({ "CCO" }, injected { "SMILES", "CCO" })
    end)

    it("reads a reaction column as SMARTS", function()
        assert.are.same({ "[C:1]>>[C:1]" },
            injected { "reaction_smarts", "[C:1]>>[C:1]" })
    end)

    it("takes the header from the first row that is not a comment", function()
        assert.are.same({ "CCO" },
            injected { "# generated by rdkit", "name\tsmiles", "ethanol\tCCO" })
    end)

    it("leaves a comment row in a chemical column alone", function()
        assert.are.same({ "CCO" },
            injected { "smiles\tname", "# CCO\tskipped", "CCO\tethanol" })
    end)

    it("leaves the header itself alone", function()
        -- A header named after a structure would otherwise inject as one.
        assert.are.same({}, injected { "smiles", "1.5" })
    end)

    -- `.bed` resolves to filetype tsv and has no header, and a number in the
    -- first row is the evidence for that.
    it("reads a numeric first row as data rather than a header", function()
        assert.are.same({ "CCO", "CCN" },
            injected({ "CCO\t100", "CCN\t200" }, function(buf)
                vim.b[buf].chem_columns = { true }
            end))
    end)

    -- The grammar drops an empty cell rather than noding it, so a column index
    -- counted in `field` siblings would shift left at every missing value.
    it("counts past an empty cell", function()
        assert.are.same({ "CCO", "CO" }, injected {
            "name\tnote\tsmiles",
            "ethanol\t\tCCO",
            "methanol\tx\tCO",
        })
    end)

    it("refuses a header with fewer fields than its data rows", function()
        assert.are.same({}, injected { "smiles", "CCO\tethanol" })
    end)

    it("drops the quotes of a quoted cell", function()
        assert.are.same({ "C(C)O" }, injected { "smiles", '"C(C)O"' })
    end)

    -- A doubled quote stands for one quote, which no range can undouble.
    it("leaves a cell holding a doubled quote alone", function()
        assert.are.same({}, injected { "smiles", '"C(""C)O"' })
    end)

    it("injects a hand-marked column", function()
        assert.are.same({ "CCO", "CCN" },
            injected({ "name\tstructure", "ethanol\tCCO", "amine\tCCN" },
                function(buf)
                    -- Sparse: the hole at index 1 comes back as the truthy
                    -- vim.NIL, which would mark the name column too.
                    vim.b[buf].chem_columns = { [2] = true }
                end))
    end)

    it("keeps injecting a column left of a hidden one", function()
        assert.are.same({ "CCO" },
            injected({ "smiles\tname", "CCO\tethanol" }, function(buf)
                -- The hole at index 1 comes back as the truthy vim.NIL.
                vim.b[buf].tsv_hidden = { [2] = {} }
            end))
    end)

    -- Two chemical columns in one row are two regions on one line, and the
    -- element marks of each are cleared by range so they do not wipe each other.
    it("colours the elements of two chemical columns in one row", function()
        with_tsv({ "smiles\tsmarts", "CCO\t[Fe]" }, function(buf, parser)
            assert.are.same({ "CCO", "[Fe]" }, structures(buf, parser))
            assert.are.same({
                "@chem.element.carbon",
                "@chem.element.iron",
                "@chem.element.oxygen",
            }, element_groups(buf))
        end)
    end)

    -- Rewriting a line collapses every mark on it to the start of the next line,
    -- where no region's columns reach: without 'invalidate' to tell those from a
    -- live mark, each edit would leave the row's marks behind for good.
    it("leaves one mark per atom after the row is rewritten", function()
        with_tsv({ "smiles\tsmarts", "CCO\t[Fe]" }, function(buf, parser)
            structures(buf, parser)
            for _ = 1, 3 do
                vim.api.nvim_buf_set_lines(buf, 1, 2, false, { "CCO\t[Fe]" })
                parser:parse(true)
            end
            local marks = vim.tbl_map(function(mark)
                return { mark[2], mark[3] }
            end, vim.api.nvim_buf_get_extmarks(buf, ns, 0, -1, {}))
            assert.are.same(
                { { 1, 0 }, { 1, 1 }, { 1, 2 }, { 1, 5 } }, marks)
        end)
    end)

    it("stops injecting a column the ftplugin hides, marks and all", function()
        with_tsv({ "smiles\tname", "[Fe]\tferrocene" }, function(buf, parser)
            assert.are.same({ "[Fe]" }, structures(buf, parser))
            assert.are.same({ "@chem.element.iron" }, element_groups(buf))
            vim.b[buf].tsv_hidden = { [1] = {} }
            parser:invalidate(true)
            assert.are.same({}, structures(buf, parser))
            assert.are.same({}, element_groups(buf))
        end)
    end)
end)

describe("chemical column mark", function()
    --- Toggle the column of a cursor position, in the current window.
    --- @param buf integer
    --- @param lnum integer 1-based
    --- @param col integer 0-based byte column
    local function toggle_at(buf, lnum, col)
        vim.api.nvim_win_set_buf(0, buf)
        vim.api.nvim_win_set_cursor(0, { lnum, col })
        chem_tsv.toggle_column()
    end

    it("makes a column its header does not name inject", function()
        with_tsv({ "structure\tname", "CCO\tethanol", "CCN\tamine" },
            function(buf, parser)
                assert.are.same({}, structures(buf, parser))
                toggle_at(buf, 2, 0)
                assert.are.same({ "CCO", "CCN" }, structures(buf, parser))
            end)
    end)

    it("survives an edit to an unrelated row", function()
        with_tsv({ "structure\tname", "CCO\tethanol", "CCN\tamine" },
            function(buf, parser)
                toggle_at(buf, 2, 0)
                vim.api.nvim_buf_set_lines(buf, 2, 3, false, { "CCN\tethylamine" })
                assert.are.same({ "CCO", "CCN" }, structures(buf, parser))
            end)
    end)

    it("counts the column in tabs, so an empty cell does not shift it", function()
        with_tsv({ "name\tnote\tstructure", "ethanol\t\tCCO" },
            function(buf, parser)
                toggle_at(buf, 2, 10)
                assert.are.same({ "CCO" }, structures(buf, parser))
            end)
    end)

    it("unmarks on a second toggle", function()
        with_tsv({ "structure\tname", "CCO\tethanol" }, function(buf, parser)
            toggle_at(buf, 2, 0)
            assert.are.same({ "CCO" }, structures(buf, parser))
            toggle_at(buf, 2, 0)
            assert.are.same({}, structures(buf, parser))
        end)
    end)
end)

describe("tsv highlights", function()
    -- The shipped query gives every cell @string, which would paint over the
    -- injected structure; this config's replaces it and keeps the rest.
    it("captures the numeric cells and no text cell", function()
        local captures = with_tsv({ "n\tx\tflag\tname", "1\t1.5\ttrue\tethanol" },
            function(buf, parser)
                local query = assert(vim.treesitter.query.get("tsv", "highlights"))
                local found = {}
                for id, node in query:iter_captures(parser:parse()[1]:root(), buf) do
                    found[vim.treesitter.get_node_text(node, buf)] = query.captures[id]
                end
                return found
            end)
        assert.are.same({
            ["1"] = "number",
            ["1.5"] = "number.float",
            ["true"] = "boolean",
        }, captures)
    end)
end)
