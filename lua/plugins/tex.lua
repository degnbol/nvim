local cmd = vim.cmd
local fn = vim.fn
local g = vim.g

-- :VimtexCompile. Adds so much more good stuff, e.g. dse, cse to delete or change surrounding env
return {
    -- more aggressive conceal
    {
        "KeitaNakamura/tex-conceal.vim",
        ft="tex",
        config=function ()
            -- https://github.com/gillescastel/latex-snippets
            g.tex_conceal='abdmg'
        end,
    },
    -- AirLatex is a plugin for editing overleaf locally where you can even see the cursors of other users.
    -- Currently doesn't work for me so I use my git auto update approach instead for now.
    -- is detected as a robot by overleaf and looking up the cookie string isn't convenient and doesn't seem to work either.
    {
        "da-h/AirLatex.vim",
        -- cmd = "AirLatex", -- breaks it
        enabled = false,
        build = {
            "pip install 'keyring[completion]' tornado requests pynvim",
            -- if password is not already set then you can use
            -- keyring set airlatex_www.overleaf.com cdmadsen@student.unimelb.edu.au
            ":UpdateRemotePlugins",
        },
        init = function ()
            vim.g.AirLatexAllowInsecure = 0
            -- vim.g.AirLatexUsername = "cdmadsen@student.unimelb.edu.au"
            vim.g.AirLatexUsername = "cookies"
        end,
    },
    {
        "lervag/vimtex",
        init=function()
            -- instead of <localleader>l
            -- We want <localleader>l for LSP and double leader for filetype specific mappings.
            g.vimtex_mappings_prefix = "<localleader><localleader>"

            -- g.vimtex_view_method = 'sioyek' -- not working
            g.vimtex_view_method = 'skim'
            -- inspiration from https://dr563105.github.io/blog/skim-vimtex-setup/
            -- forward search after every successful compilation
            g.vimtex_view_skim_sync = true
            -- change focus to skim after command `:VimtexView` is given
            g.vimtex_view_skim_activate = true
            g.vimtex_view_skim_reading_bar = true

            g.vimtex_log_ignore = { "Underfull", "Overfull" }
            g.vimtex_quickfix_ignore_filters = { "Underfull", "Overfull" }
            g.vimtex_quickfix_open_on_warning = false

            -- :h vimtex-af-enhanced-matchparen
            -- wrong matching for some {} with this on:
            -- vim.g.matchup_override_vimtex = true

            -- defaults at :h vimtex_compiler_latexmk
            g.vimtex_compiler_latexmk = {
                aux_dir = "aux",
            }
            -- set default latex engine to the modern lualatex over pdflatex
            g.vimtex_compiler_latexmk_engines = {
                _="-lualatex",
            }

            -- formatter when calling gq, but note that autoformatting calls the 
            -- builtin vim formatter which you can call with gw.
            -- In order to not autoformat in math mode I made an autocmd to 
            -- disable/enable auto formatoption based on vimtex's detection in_mathzone.
            -- See ftplugin/tex.vim
            g.vimtex_format_enabled = true
            -- trying to avoid indent after \item\n. Why? when the "\item" is 
            -- concealed the text doesn't align nicely which must have been the 
            -- original intent. Also, my snippet expanding a '-' for a new item becomes more complicated.
            -- vimtex solution doesn't work: https://github.com/lervag/vimtex/issues/2599
            -- use builtin indent function GetTexIndent since the solution for that does work:
            -- https://github.com/vim/vim/blob/eb3dc87f01391bb075d97aef3d00f91b4e08a25c/runtime/indent/tex.vim#L70-L122
            g.vimtex_indent_enabled = false
            -- builtin vim setting for tex files, not a vimtex setting
            g.tex_indent_items = 0

            g.vimtex_toc_config = {
                indent_levels = 1, -- has no effect but would be nice
                mode = 4, -- the only mode that works
                -- layers = {'content', 'todo', 'include'}, -- don't list labels
            }

            -- if in a subfile, by default we compile only that.
            -- <leader><leader>m to toggle compiling main instead.
            g.vimtex_subfile_start_local = true -- doesn't seem to work.
            -- Also toggling is inconsistent so just use it and look at the message.

            -- Since we use function as bold func def and function.call as unbold, we relink:
            vim.defer_fn(function ()
                vim.api.nvim_set_hl(0, "texCmd", {link="@function.call", force=true})
            end, 1000)

            local grp = vim.api.nvim_create_augroup("TexMain", {clear=true})
            -- Set tex main file by looking for a file upwards named "main.tex"
            vim.api.nvim_create_autocmd("BufReadPre", {
                pattern = "*.tex",
                group = grp,
                callback = function ()
                    local mainfile = vim.fs.find("main.tex", {type="file", upward=true})[1]
                    if mainfile ~= nil then vim.b.vimtex_main = mainfile end
                end
            })
        end,
    },
}

