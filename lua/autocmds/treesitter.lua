-- start treesitter for each new filetype
vim.api.nvim_create_autocmd("FileType", {
    pattern = "*",
    group = vim.api.nvim_create_augroup("start_treesitter", { clear = false }),
    callback = function(args)
        local disabled = {
            "vim",                                   -- not perfect
            "latex",                                 -- messes with vimtex in lots of ways, e.g. conceal, detection of mathzone, cycling with ts$
            "sh", "bash", "zsh", "sh.zsh", "zsh.sh", -- broken
        }
        local additional_vim_regex_highlighting = {
            "vimdoc",   -- treesitter version doesn't contain useful colors from :h group-name
            "bash",     -- spending too much time writing treesitter query. Also covers zsh.
            "markdown", -- my custom comment syntax matches in after/syntax/markdown.vim
            -- Semicolon isn't currently highlighted in all cases by TS so we want to incl vim regex hl for jl.
            -- However, jl can get slowed down a lot in certain files from the syntax hl. The solution:
            -- We enable it, but avoid any default syntax hl and only set custom syntax hl in syntax/julia.vim.
            "julia",
            "sql",  -- custom postgres highlight in syntax/sql.vim
            "wgsl", -- custom in syntax/wgsl.vim
        }
        if not vim.list_contains(disabled, vim.bo.filetype) then
            -- Ignore errors for buffers where treesitter parser is not installed.
            pcall(vim.treesitter.start, args.buf)
            if additional_vim_regex_highlighting[vim.bo.filetype] then
                vim.bo[args.buf].syntax = 'on'
            end
        end
    end
})
