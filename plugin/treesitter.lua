-- use the zsh parser for all shell filetypes. tree-sitter-bash is unmaintained
-- and has cascading error bugs (parameter expansion, extglob, heredocs).
-- zsh is a superset of bash, so the zsh parser handles bash files fine.
vim.treesitter.language.register("zsh", "sh")
vim.treesitter.language.register("zsh", "bash")
vim.treesitter.language.register("zsh", "sh.zsh")
-- ```scm fenced code blocks in markdown should use the query parser
vim.treesitter.language.register("query", "scm")
-- jsonl/ndjson is newline-delimited json; the json grammar's document rule
-- accepts repeated top-level values, so it highlights every record.
vim.treesitter.language.register("json", "jsonl")
-- SMILES is very nearly a syntactic subset of SMARTS, so one grammar covers
-- both filetypes. Also makes ```smiles fences resolve, since language.get_lang
-- consults the same alias table.
vim.treesitter.language.register("smarts", "smiles")
-- typc = typst code mode. tinymist tags hover/completion code fences ```typc;
-- no typc parser exists, so map it to the typst parser for markdown injection
-- (e.g. LSP hover floats). Caveat: the typst grammar parses markup-first, so
-- bare code-mode lines (`let x = …`) stay plain — only `#`-prefixed
-- expressions and the hex constant get captured.
vim.treesitter.language.register("typst", "typc")

local ts_utils = require "utils/treesitter"
vim.treesitter.query.add_directive(
    "trim!", ts_utils.trim_directive, { force = true })
vim.treesitter.query.add_directive(
    "head!", ts_utils.head_directive, { force = true })
vim.treesitter.query.add_directive(
    "tail!", ts_utils.tail_directive, { force = true })
vim.treesitter.query.add_directive(
    "unquote!", ts_utils.unquote_directive, { force = true })
vim.treesitter.query.add_directive(
    "inject-by-ext!", ts_utils.inject_by_ext_directive, { force = true })
vim.treesitter.query.add_directive(
    "inject-interp!", ts_utils.inject_interp_directive, { force = true })
vim.treesitter.query.add_directive(
    "inject-interp-cmd!", ts_utils.inject_interp_cmd_directive, { force = true })
vim.treesitter.query.add_predicate(
    "command-is?", ts_utils.command_is, { force = true })
vim.treesitter.query.add_predicate(
    "arg-after?", ts_utils.arg_after, { force = true })
vim.treesitter.query.add_predicate(
    "regex-pattern?", ts_utils.is_regex_pattern, { force = true })
-- Registered here with the rest and not in plugin/chem.lua, because a query
-- naming a predicate no one has registered does not degrade: the parse throws
-- and the buffer goes unhighlighted.
vim.treesitter.query.add_predicate(
    "chem-column?", require("chem/tsv").is_chem_column, { force = true })

local disabled = {
    -- messes with vimtex in lots of ways, e.g. conceal, detection of mathzone, cycling with ts$,
    "latex", "plaintex", "tex",
}
local additional_vim_regex_highlighting = {
    "vimdoc",   -- treesitter version doesn't contain useful colors from :h group-name
    "sh", "bash", "sh.zsh",
    "markdown", -- my custom comment syntax matches in after/syntax/markdown.vim
    -- Semicolon isn't currently highlighted in all cases by TS so we want to incl vim regex hl for jl.
    -- However, jl can get slowed down a lot in certain files from the syntax hl. The solution:
    -- We enable it, but avoid any default syntax hl and only set custom syntax hl in syntax/julia.vim.
    "julia",
    "sql",  -- custom postgres highlight in syntax/sql.vim
    "wgsl", -- custom in syntax/wgsl.vim
}

--- Start treesitter highlighting on a buffer, and announce that it did.
--- @param buf integer
local function start_treesitter(buf)
    -- Most filetypes have no grammar, which is a capability and not an error, so
    -- ask for the parser rather than catch a throw. With one built, start cannot
    -- fail: it is get_parser plus highlighter.new.
    if not vim.treesitter.get_parser(buf) then return end
    vim.treesitter.start(buf)
    vim.api.nvim_exec_autocmds("User", {
        pattern = "TSHighlightStart", data = { buf = buf }, modeline = false,
    })
end

-- start treesitter for each new filetype
vim.api.nvim_create_autocmd("FileType", {
    pattern = "*",
    group = vim.api.nvim_create_augroup("start_treesitter", { clear = false }),
    callback = function(args)
        if vim.b[args.buf].largefile then return end
        -- A plugin that highlights a buffer itself claims it by setting
        -- 'ts_highlight' before it sets 'filetype' — see the OPTIM comment in
        -- snacks.nvim's win.lua, which goes on to highlight the language it means
        -- rather than the one its scratch filetype names. Leave the buffer to it,
        -- vim regex syntax included: setting 'syntax' here would stop it too.
        if vim.b[args.buf].ts_highlight then return end
        if vim.list_contains(disabled, args.match) then return end
        -- WORKAROUND: nvim 0.12 bundled markdown parser crashes during initial load.
        -- Delay treesitter.start for markdown to after buffer is fully set up.
        if args.match == "markdown" then
            vim.schedule(function()
                -- The buffer can be gone by the time this runs, and get_parser
                -- throws on an invalid one.
                if vim.api.nvim_buf_is_valid(args.buf) then start_treesitter(args.buf) end
            end)
        else
            start_treesitter(args.buf)
        end
        if vim.list_contains(additional_vim_regex_highlighting, args.match) then
            vim.bo[args.buf].syntax = 'on'
        end
    end
})
