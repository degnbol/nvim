local cmd = vim.cmd
local fn = vim.fn
local g = vim.g

-- :VimtexCompile. Adds so much more good stuff, e.g. dse, cse to delete or change surrounding env
return {
    -- more aggressive conceal
    {"KeitaNakamura/tex-conceal.vim", config=function ()
        -- https://github.com/gillescastel/latex-snippets
        g.tex_conceal='abdmg'
    end},
    {"lervag/vimtex", config=function() 
        g.vimtex_view_method = 'skim'
        -- inspiration from https://dr563105.github.io/blog/skim-vimtex-setup/
        -- forward search after every successful compilation
        g.vimtex_view_skim_sync = true
        -- change focus to skim after command `:VimtexView` is given
        g.vimtex_view_skim_activate = true
        
        g.vimtex_log_ignore = { "Underfull", "Overfull" }
        g.vimtex_quickfix_open_on_warning = false

        -- :h vimtex-af-enhanced-matchparen
        -- wrong matching for some {} with this on:
        -- vim.g.matchup_override_vimtex = true

        -- defaults at :h vimtex_compiler_latexmk
        g.vimtex_compiler_latexmk = {
            build_dir = "build",
        }
        -- set default latex engine to the modern lualatex over pdflatex
        g.vimtex_compiler_latexmk_engines = { _="-lualatex" }

        g.vimtex_syntax_conceal = {
            sections = true, -- all other conceals are enabled by default
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
        -- use builtin indent function GetTexIndent() since the solution for that does work:
        -- https://github.com/vim/vim/blob/eb3dc87f01391bb075d97aef3d00f91b4e08a25c/runtime/indent/tex.vim#L70-L122
        g.vimtex_indent_enabled = false
        -- builtin vim setting for tex files, not a vimtex setting
        g.tex_indent_items = 0
        
        g.vimtex_toc_config = {
            indent_levels = 1, -- has no effect but would be nice
            mode = 4, -- the only mode that works
            -- layers = {'content', 'todo', 'include'}, -- don't list labels
        }

    end},
}

