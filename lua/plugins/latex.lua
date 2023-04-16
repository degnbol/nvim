-- :VimtexCompile. Adds so much more good stuff, e.g. dse, cse to delete or change surrounding env
return {
    -- more aggressive conceal
    {"KeitaNakamura/tex-conceal.vim", config=function ()
        g = vim.g
        -- https://github.com/gillescastel/latex-snippets
        g.tex_conceal='abdmg'
    end},
    {"lervag/vimtex", config=function() 
        g = vim.g

        g.vimtex_view_method = 'skim'
        -- inspiration from https://dr563105.github.io/blog/skim-vimtex-setup/
        -- forward search after every successful compilation
        g.vimtex_view_skim_sync = true
        -- change focus to skim after command `:VimtexView` is given
        g.vimtex_view_skim_activate = true

        -- don't show log for every compilation.
        g.vimtex_quickfix_mode = 0
        g.vimtex_log_ignore = { "Underfull", "Overfull", }

        -- :h vimtex-af-enhanced-matchparen
        -- wrong matching for some {} with this on:
        -- vim.g.matchup_override_vimtex = true

        -- defaults at :h vimtex_compiler_latexmk
        g.vimtex_compiler_latexmk = { build_dir = "build", }

        g.vimtex_syntax_conceal = {
            sections = true, -- all other conceals are enabled by default
        }
        
        -- formatter when calling gq, but note that autoformatting calls the 
        -- builtin vim formatter which you can call with gw.
        -- In order to not autoformat in math mode I made an autocmd to 
        -- disable/enable auto formatoption based on vimtex's detection in_mathzone.
        -- See ftplugin/tex.vim
        g.vimtex_format_enabled = true
    end},
}

