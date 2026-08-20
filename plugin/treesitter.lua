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
-- Beside the parser wiring rather than in an ftplugin: a ```smiles fence in a
-- markdown buffer is coloured by the buffer's own attach below but has no
-- filetype of its own, so an ftplugin would leave the groups undefined.
local chem = require "chem/highlight"
chem.setup()
-- typc = typst code mode. tinymist tags hover/completion code fences ```typc;
-- no typc parser exists, so map it to the typst parser for markdown injection
-- (e.g. LSP hover floats). Caveat: the typst grammar parses markup-first, so
-- bare code-mode lines (`let x = …`) stay plain — only `#`-prefixed
-- expressions and the hex constant get captured.
vim.treesitter.language.register("typst", "typc")

local ts_utils = require("utils.treesitter")
vim.treesitter.query.add_directive(
    "trim!", ts_utils.trim_directive, { force = true })
vim.treesitter.query.add_directive(
    "head!", ts_utils.head_directive, { force = true })
vim.treesitter.query.add_directive(
    "tail!", ts_utils.tail_directive, { force = true })
vim.treesitter.query.add_directive(
    "inject-by-ext!", ts_utils.inject_by_ext_directive, { force = true })
vim.treesitter.query.add_directive(
    "inject-interp!", ts_utils.inject_interp_directive, { force = true })
vim.treesitter.query.add_directive(
    "inject-interp-cmd!", ts_utils.inject_interp_cmd_directive, { force = true })
vim.treesitter.query.add_predicate(
    "any-basename-of?", ts_utils.any_basename_of, { force = true })

-- Element colours ride on the highlighter's own decision rather than restating
-- it: 'ts_highlight' is false for a largefile buffer and for the filetypes
-- disabled below, and the initial repaint covers the whole buffer.
local function start_treesitter(buf)
    pcall(vim.treesitter.start, buf)
    if vim.b[buf].ts_highlight then chem.attach(buf) end
end

-- start treesitter for each new filetype
vim.api.nvim_create_autocmd("FileType", {
    pattern = "*",
    group = vim.api.nvim_create_augroup("start_treesitter", { clear = false }),
    callback = function(args)
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
        if vim.b[args.buf].largefile then return end
        if not vim.list_contains(disabled, vim.bo.filetype) then
            -- WORKAROUND: nvim 0.12 bundled markdown parser crashes during initial load.
            -- Delay treesitter.start for markdown to after buffer is fully set up.
            if vim.bo.filetype == "markdown" then
                vim.schedule(function() start_treesitter(args.buf) end)
            else
                start_treesitter(args.buf)
            end
            if vim.list_contains(additional_vim_regex_highlighting, vim.bo.filetype) then
                vim.bo[args.buf].syntax = 'on'
            end
        end
    end
})
