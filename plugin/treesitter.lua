-- use the zsh parser for all shell filetypes. tree-sitter-bash is unmaintained
-- and has cascading error bugs (parameter expansion, extglob, heredocs).
-- zsh is a superset of bash, so the zsh parser handles bash files fine.
vim.treesitter.language.register("zsh", "sh")
vim.treesitter.language.register("zsh", "bash")
vim.treesitter.language.register("zsh", "sh.zsh")
-- ```scm fenced code blocks in markdown should use the query parser
vim.treesitter.language.register("query", "scm")

local ts_utils = require("utils.treesitter")
vim.treesitter.query.add_directive(
    "trim!", ts_utils.trim_directive, { force = true })
vim.treesitter.query.add_directive(
    "head!", ts_utils.head_directive, { force = true })
vim.treesitter.query.add_directive(
    "tail!", ts_utils.tail_directive, { force = true })

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
                vim.schedule(function()
                    pcall(vim.treesitter.start, args.buf)
                end)
            else
                pcall(vim.treesitter.start, args.buf)
            end
            if vim.list_contains(additional_vim_regex_highlighting, vim.bo.filetype) then
                vim.bo[args.buf].syntax = 'on'
            end
        end
    end
})
