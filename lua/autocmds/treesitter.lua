-- compound filetype "sh.zsh" defaults to first component for treesitter parser.
-- register so zsh files use the dedicated zsh parser instead of bash.
vim.treesitter.language.register("zsh", "sh.zsh")

-- start treesitter for each new filetype
vim.api.nvim_create_autocmd("FileType", {
    pattern = "*",
    group = vim.api.nvim_create_augroup("start_treesitter", { clear = false }),
    callback = function(args)
        local disabled = {
            -- not perfect
            "vim",
            -- messes with vimtex in lots of ways, e.g. conceal, detection of mathzone, cycling with ts$,
            "latex", "plaintex", "tex",
            -- broken
            "sh", "bash",
        }
        local additional_vim_regex_highlighting = {
            "vimdoc",   -- treesitter version doesn't contain useful colors from :h group-name
            "bash",     -- spending too much time writing treesitter query.
            "sh.zsh",
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
