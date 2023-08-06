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

        vim.defer_fn(function ()
            require "utils/keymap"
            set_keymap_desc('n', '<leader>la', "Context menu")
            -- align table see ftplugin/tex.lua. ma ... `a to not move cursor.
            vim.keymap.set('n', '<leader>lA', "ma<plug>AlignTable<CR>`a", { desc="Align table" })
            vim.keymap.set("n", "<leader>lb", "<Cmd>Telescope bibtex<CR>", { desc="Cite" })
            set_keymap_desc('n', '<leader>lc', "Clean")
            set_keymap_desc('n', '<leader>lC', "Clean all")
            set_keymap_desc('n', '<leader>le', "Errors")
            set_keymap_desc('n', '<leader>lg', "Status")
            set_keymap_desc('n', '<leader>lG', "Status all")
            set_keymap_desc('n', '<leader>li', "Info")
            set_keymap_desc('n', '<leader>lI', "Info full")
            set_keymap_desc('n', '<leader>lq', "Log")
            -- j for jump. Not using p for preamble since I want to use it for pasting tables.
            vim.keymap.set("n", "<leader>lj", "<Plug>TexJumpPre", { desc="goto/from preamble (table)" })
            set_keymap_desc('n', '<leader>ls', "Toggle main")
            set_keymap_desc('n', '<leader>lt', "TOC open")
            set_keymap_desc('n', '<leader>lT', "TOC toggle")
            vim.keymap.set({"n", "x"}, "<leader>lu", "<Plug>Latex2Unicode", { desc="TeX -> unicode" })
            vim.keymap.set({"n", "x"}, "<leader>lU", "<Plug>Unicode2Latex", { desc="Unicode -> TeX" })
            set_keymap_desc('n', '<leader>lv', "View")
            set_keymap_desc('n', '<leader>lx', "Reload")
            set_keymap_desc('n', '<leader>lX', "Reload state")
            vim.keymap.set("n", "<leader>ly", "<plug>YankTable", { desc="Yank as TSV" })
            set_keymap_desc('n', '<leader>lk', "Stop")
            set_keymap_desc('n', '<leader>lK', "Stop all")
            set_keymap_desc('n', '<leader>ll', "Compile")
            set_keymap_desc('n', '<leader>lL', "Compile selected")
            -- TODO: maybe take inspo from these insert mode mappings and make 
            -- snippet equivalents then disable them? Or use them if they are useful.
            set_keymap_desc('n', '<leader>lm', "imaps list")
            set_keymap_desc('n', '<leader>lo', "Raw compl output")

            -- e.g. \section*{}
            set_keymap_desc('n', 'tsc', "Cmd/Star")
            set_keymap_desc('n', 'tse', "Env/Star")
            -- e.g. with(out) \left 
            set_keymap_desc('n', 'tsd', "Delim")
            -- same as d, but looks through g:vimtex_delim_toggle_mod_list in reverse
            set_keymap_desc('n', 'tsD', "Delim rev")
            -- toggle / <-> \frac
            set_keymap_desc('n', 'tsf', "Fraction")
            -- change surrounding ...
            set_keymap_desc('n', 'csc', "Cmd")
            set_keymap_desc('n', 'cse', "Env")
            set_keymap_desc('n', 'csm', "Math")
            -- in/around ...
            -- TODO: P is currently conflicting with paste from unimpaired
            set_keymap_desc({'o', 'x'}, 'id', "Delim")
            set_keymap_desc({'o', 'x'}, 'iP', "Section")
            vim.keymap.set("n", "ts$", "<Plug>(vimtex-env-toggle-math)", { desc="Inline <-> display" })
            vim.keymap.set("n", "ts4", "<Plug>(vimtex-env-toggle-math)", { desc="Inline <-> display" })
            vim.keymap.set("n", "tsm", "<plug>(vimtex-env-toggle-math)", { desc="Inline <-> display"})
            vim.keymap.set("n", "dsm", "<plug>(vimtex-env-delete-math)", { desc="Delete math"})
            vim.keymap.set("n", "csm", "<plug>(vimtex-env-change-math)", { desc="Change math"})
            vim.keymap.set("n", "xad", "yaddad", {remap=true, desc="Cut a delim"})
            vim.keymap.set("n", "xid", "yiddid", {remap=true, desc="Cut in delim"})
            -- item with i instead of m and math with m
            vim.keymap.set({"o", "x"}, "ai", "<Plug>(vimtex-am)", {desc="An item"})
            vim.keymap.set({"o", "x"}, "ii", "<Plug>(vimtex-im)", {desc="In item"})
            vim.keymap.set({"o", "x"}, "am", "<Plug>(vimtex-a$)", {desc="An eq"})
            vim.keymap.set({"o", "x"}, "im", "<Plug>(vimtex-i$)", {desc="In eq"})
            -- shorthand to $ just using 4 ($ without shift)
            vim.keymap.set({"o", "x"}, "a4", "<Plug>(vimtex-a$)", {desc="An eq"})
            vim.keymap.set({"o", "x"}, "i4", "<Plug>(vimtex-i$)", {desc="In eq"})
            -- next/prev start/end of ...
            -- TODO: the m ones are conflicting with something
            vim.keymap.set("n", "[m", "<Plug>(vimtex-[n)", { desc="Math start" })
            vim.keymap.set("n", "[M", "<Plug>(vimtex-[N)", { desc="Math end" })
            vim.keymap.set("n", "[4", "<Plug>(vimtex-[n)", { desc="Math start" })
            vim.keymap.set("n", "[$", "<Plug>(vimtex-[N)", { desc="Math end" })
            vim.keymap.set("n", "]m", "<Plug>(vimtex-]n)", { desc="Math start" })
            vim.keymap.set("n", "]M", "<Plug>(vimtex-]N)", { desc="Math end" })
            vim.keymap.set("n", "]4", "<Plug>(vimtex-]n)", { desc="Math start" })
            vim.keymap.set("n", "]$", "<Plug>(vimtex-]N)", { desc="Math end" })
        end, 0)
    end},
}

