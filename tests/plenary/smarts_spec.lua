---@diagnostic disable: undefined-global
-- Tests for the smarts grammar's nvim wiring: plugin/ftdetect.lua,
-- plugin/treesitter.lua, plugin/chem.lua and modules/tree-sitter-smarts.

-- plugin/* isn't reliably sourced before PlenaryBustedFile runs.
vim.cmd.source(vim.fn.getcwd() .. "/plugin/ftdetect.lua")
vim.cmd.source(vim.fn.getcwd() .. "/plugin/treesitter.lua")
vim.cmd.source(vim.fn.getcwd() .. "/plugin/chem.lua")

local chem = require "chem/highlight"
-- Same call the module makes, so the same id: nvim_create_namespace is a lookup
-- for a name already taken.
local ns = vim.api.nvim_create_namespace("chem.elements")

--- Run a function on a scratch buffer holding lines of one filetype.
--- @param lines string[]
--- @param filetype string
--- @param fn fun(buf: integer, parser: vim.treesitter.LanguageTree): any
--- @return any
local function with_buffer(lines, filetype, fn)
    local buf = vim.api.nvim_create_buf(false, true)
    vim.api.nvim_buf_set_lines(buf, 0, -1, false, lines)
    vim.bo[buf].filetype = filetype
    local result = fn(buf, assert(vim.treesitter.get_parser(buf)))
    vim.api.nvim_buf_delete(buf, { force = true })
    return result
end

--- Parse text as the given filetype and collect the ranges of its ERROR nodes.
--- @param src string
--- @param filetype string
--- @return string[] descriptions one per ERROR or MISSING node
local function parse_errors(src, filetype)
    return with_buffer(vim.split(src, "\n", { plain = true }), filetype, function(buf, parser)
        local errors = {}
        for node in parser:parse()[1]:root():iter_children() do
            if node:has_error() then
                local sr, sc = node:range()
                errors[#errors + 1] = ("%s at %d,%d: %s")
                    :format(node:type(), sr, sc, vim.treesitter.get_node_text(node, buf))
            end
        end
        return errors
    end)
end

--- The members of one symbol set the grammar's generated alphabet declares.
--- @param name string the `const` the set is declared as
--- @return string the bracketed member list, verbatim
local function alphabet(name)
    local source = table.concat(vim.fn.readfile(
        vim.fn.getcwd() .. "/modules/tree-sitter-smarts/symbols.js"), "\n")
    return assert(source:match("const " .. name .. " = %[(.-)%]"))
end

--- The distinct @chem.* groups each glyph of a buffer's structures is given, by
--- either mechanism: a capture name for what a node's type says, an extmark for
--- what its text says. Distinct because a glyph recurs across a structure and
--- every occurrence of it takes the same groups.
--- @param buf integer buffer the trees were parsed from
--- @param roots TSNode[] roots of the buffer's smarts trees
--- @return table<string, string[]> glyph text -> its groups
local function chem_groups(buf, roots)
    local groups = {}
    local function add(text, name)
        groups[text] = groups[text] or {}
        if not vim.list_contains(groups[text], name) then
            table.insert(groups[text], name)
        end
    end

    local query = assert(vim.treesitter.query.get("smarts", "highlights"))
    for _, root in ipairs(roots) do
        for id, node in query:iter_captures(root, buf) do
            local name = query.captures[id]
            if vim.startswith(name, "chem.") then
                add(vim.treesitter.get_node_text(node, buf), name)
            end
        end
    end
    for _, mark in ipairs(vim.api.nvim_buf_get_extmarks(buf, ns, 0, -1, { details = true })) do
        local row, col, details = mark[2], mark[3], assert(mark[4])
        local text = vim.api.nvim_buf_get_text(
            buf, row, col, assert(details.end_row), assert(details.end_col), {})[1]
        -- Capture names have no leading `@`, so drop it and compare either way.
        add(text, details.hl_group:sub(2))
    end
    return groups
end

--- The @chem.* groups of a one-line structure.
--- @param src string
--- @param filetype string
--- @return table<string, string[]>
local function structure_groups(src, filetype)
    return with_buffer({ src }, filetype, function(buf, parser)
        chem.attach(buf)
        return chem_groups(buf, { parser:parse()[1]:root() })
    end)
end

--- The @chem.* groups of a structure inside a markdown fence.
--- @param src string
--- @return table<string, string[]>
local function fenced_groups(src)
    return with_buffer({ "```smiles", src, "```" }, "markdown", function(buf, parser)
        chem.attach(buf)
        parser:parse(true)
        local roots = {}
        for _, tree in pairs(assert(parser:children().smarts):trees()) do
            roots[#roots + 1] = tree:root()
        end
        return chem_groups(buf, roots)
    end)
end

describe("chemical line notation filetypes", function()
    -- `.smi` resolves to the built-in `mib` without the ftdetect entry.
    local extensions = {
        smi = "smiles",
        smiles = "smiles",
        cxsmiles = "smiles",
        rsmi = "smiles",
        smarts = "smarts",
        sma = "smarts",
        smirks = "smarts",
    }

    for extension, filetype in pairs(extensions) do
        it("detects ." .. extension .. " as " .. filetype, function()
            assert.are.equal(filetype, vim.filetype.match { filename = "mol." .. extension })
        end)
    end

    it("resolves both filetypes to the smarts parser", function()
        assert.are.equal("smarts", vim.treesitter.language.get_lang("smiles"))
        assert.are.equal("smarts", vim.treesitter.language.get_lang("smarts"))
    end)

    -- A fence's info string is its injection language, so the alias reaches
    -- embedded blocks with no injection query of our own. Typst raw blocks name
    -- their language the same way and resolve through the same table.
    it("injects the parser into a markdown fence", function()
        local injected = with_buffer({ "```smiles", "CCO ethanol", "```" }, "markdown",
            function(_, parser)
                parser:parse(true)
                return parser:children().smarts
            end)
        assert.is_not_nil(injected)
    end)
end)

describe("smarts parser", function()
    it("parses a .smi file without error", function()
        local src = table.concat({
            "SMILES\tName",
            "# a comment",
            "",
            "CCO ethanol",
            "c1ccccc1\tbenzene",
            "[O-]C(=O)c1ccccc1 |$_R1;;;;;;;;$| benzoate",
            "[$([N;H2,H1;!$(N-[C,S]=[O,S,N])])]",
            "[C:1][O:2]>>[C:1][O:2]",
            "N->[Pt]<-N platinum diammine",
            "CCO |(-1.29904,-0.25,;0,0.5,;1.29904,-0.25,)|",
            "CCCC |SgD:3,2,1,0:pH:7:>=:mol/L:t:(1.,1.),atomProp:0.k.v|",
            "*-CCCN-* |$star_e;;;;;star_e$,Sg:n:1,2,3::hh&#44;f|",
            "[CH2]C |^1:0,rb:0:*,wU:1.0|",
        }, "\n")
        assert.are.same({}, parse_errors(src, "smiles"))
    end)

    -- nvim's parser reads buffer lines, so the last line arrives without a
    -- trailing newline; `source_file` has to accept a final unterminated line.
    it("parses a final line with no trailing newline", function()
        assert.are.same({}, parse_errors("CCO ethanol", "smiles"))
    end)
end)

describe("element layer", function()
    local elements = require("chem/elements")

    it("gives a symbol and its atomic number the same group", function()
        local groups = structure_groups("[Fe][#26]", "smiles")
        assert.are.same({ "chem.element.iron" }, groups["Fe"])
        assert.are.same({ "chem.element.iron" }, groups["#26"])
    end)

    -- `Fe` is not the example: bare, it lexes as `F` plus an ERROR, and on line
    -- one as a header field.
    it("reads a two-letter symbol outside brackets as two atoms", function()
        local groups = structure_groups("Sn", "smiles")
        assert.are.same({ "chem.element.sulfur" }, groups["S"])
        assert.are.same({ "chem.element.nitrogen" }, groups["n"])
        assert.is_nil(groups["Sn"])
    end)

    -- The two-letter bracket aromatics (`si`, `as`, `se`, `te`) have no bare
    -- spelling, so they reach the layer through (aromatic_symbol) alone.
    it("reads a bracket aromatic as its element", function()
        assert.are.same({ "chem.element.selenium" },
            structure_groups("c1cc[se]c1", "smiles")["se"])
    end)

    it("warns on an element no structure can contain", function()
        local groups = structure_groups("[Og][#118]", "smiles")
        assert.are.same({ "chem.suspect" }, groups["Og"])
        assert.are.same({ "chem.suspect" }, groups["#118"])
    end)

    -- `[H]` is the hydrogen atom and `[CH3]` an H-count, both spelled `H`, and
    -- neither is an element_symbol — so the atom comes from the query layer's
    -- structural patterns rather than from the text lookup.
    it("colours the hydrogen atom but not a hydrogen count", function()
        assert.are.same({ "chem.element.hydrogen" }, structure_groups("[H]", "smiles")["H"])
        assert.is_nil(structure_groups("[CH3]", "smiles")["H"])
    end)

    it("colours the atomic-number spelling of hydrogen", function()
        assert.are.same({ "chem.element.hydrogen" }, structure_groups("[#1]", "smiles")["#1"])
    end)

    -- The two mechanisms have to stay disjoint: an `H` reached by both would be
    -- an element colour depending on which extmark was set last.
    it("gives every glyph at most one element group", function()
        local groups = structure_groups("[Fe]c1ccccc1[O-]S(=O)[H][#1]", "smiles")
        for glyph, names in pairs(groups) do
            local elemental = vim.tbl_filter(function(name)
                return vim.startswith(name, "chem.element.")
            end, names)
            assert(#elemental <= 1, ("%s takes %s"):format(glyph, vim.inspect(elemental)))
        end
    end)

    -- Marks persist between redraws, so a reparse has to clear its range before
    -- re-adding, or every edit leaves another mark on the line.
    it("replaces the marks of an edited structure", function()
        with_buffer({ "CCO" }, "smiles", function(buf, parser)
            chem.attach(buf)
            parser:parse()
            vim.api.nvim_buf_set_lines(buf, 0, 1, false, { "CCN" })
            parser:parse()
            assert.are.equal(3, #vim.api.nvim_buf_get_extmarks(buf, ns, 0, -1, {}))
            assert.are.same({ "chem.element.nitrogen" }, chem_groups(buf,
                { parser:parse()[1]:root() })["N"])
        end)
    end)

    -- A root smarts tree covers the whole document as one range whose end row is
    -- a sentinel far past the buffer, and a reparse of it still has to clear
    -- what it is about to repaint.
    it("replaces the marks of a reparsed root structure", function()
        with_buffer({ "CCO" }, "smiles", function(buf, parser)
            chem.attach(buf)
            parser:parse()
            parser:invalidate(true)
            parser:parse()
            assert.are.equal(3, #vim.api.nvim_buf_get_extmarks(buf, ns, 0, -1, {}))
        end)
    end)

    -- A fence has no filetype of its own, so an ftplugin would not reach it.
    it("reaches a markdown fence, groups and all", function()
        assert.are.same({ "chem.element.oxygen" }, fenced_groups("CCO")["O"])
        assert.is_not.same({}, vim.api.nvim_get_hl(0, { name = "@chem.element.oxygen" }))
    end)

    -- The ordinary state of a markdown file, not an edge case: no fence, and so
    -- no smarts tree for the recursive callback to have been inherited by.
    it("leaves a markdown buffer with no structure alone", function()
        with_buffer({ "# heading", "prose" }, "markdown", function(buf, parser)
            chem.attach(buf)
            parser:parse(true)
            assert.are.same({}, vim.api.nvim_buf_get_extmarks(buf, ns, 0, -1, {}))
        end)
    end)

    -- Both what elements.lua states and the alphabet it is generated from.
    it("covers exactly the symbols the grammar can spell", function()
        local spellable = {}
        for symbol in alphabet("bracketElement"):gmatch('"(%a+)"') do
            spellable[symbol] = true
        end
        for symbol in pairs(elements.suspect) do
            assert(spellable[symbol], symbol .. " is warned about but cannot be spelled")
            spellable[symbol] = nil
        end
        for _, element in ipairs(elements) do
            assert(spellable[element.symbol], element.symbol .. " is coloured but not spellable")
            spellable[element.symbol] = nil
        end
        assert.are.same({}, vim.tbl_keys(spellable))
    end)

    -- The one thing the lookup needs from the grammar beyond the symbol itself:
    -- only an element with a legal lowercase spelling gets a lowercase key.
    it("records exactly the aromatic spellings the grammar accepts", function()
        local aromatic = {}
        for _, set in ipairs { "bareAromatic", "bracketAromatic" } do
            for spelling in alphabet(set):gmatch('"(%a+)"') do
                aromatic[spelling] = true
            end
        end
        for _, element in ipairs(elements) do
            assert.are.equal(aromatic[element.symbol:lower()], element.aromatic,
                element.name .. " disagrees with the grammar about an aromatic spelling")
        end
    end)
end)

--- Whether setting a filetype makes plugin/treesitter.lua announce that it
--- started treesitter highlighting — the event the element layer attaches on.
--- @param filetype string
--- @param claimed boolean whether another plugin marks the buffer as its own first
--- @return boolean
local function announces(filetype, claimed)
    local buffers = {}
    local id = vim.api.nvim_create_autocmd("User", {
        pattern = "TSHighlightStart",
        callback = function(args) buffers[#buffers + 1] = args.data.buf end,
    })
    local buf = vim.api.nvim_create_buf(false, true)
    if claimed then vim.b[buf].ts_highlight = true end
    vim.bo[buf].filetype = filetype
    vim.api.nvim_del_autocmd(id)
    local announced = vim.list_contains(buffers, buf)
    vim.api.nvim_buf_delete(buf, { force = true })
    return announced
end

describe("treesitter highlighting", function()
    it("announces the buffers it starts", function()
        assert.is_true(announces("smiles", false))
    end)

    -- Most filetypes have no grammar, and the scratch filetypes of a plugin's
    -- own windows are the ones that reach here most often.
    it("stays quiet about a filetype with no parser", function()
        assert.is_false(announces("snacks_layout_box", false))
    end)

    -- snacks.win sets 'ts_highlight' before 'filetype' precisely so that
    -- FileType handlers leave the buffer to it.
    it("stays quiet about a buffer another plugin has claimed", function()
        assert.is_false(announces("smiles", true))
    end)
end)
