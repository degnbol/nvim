#!/usr/bin/env lua
g = vim.g

-- try out skim, otherwise preview would open.
g.vimtex_view_method = 'skim'
-- inspiration from https://dr563105.github.io/blog/skim-vimtex-setup/
-- forward search after every successful compilation
g.vimtex_view_skim_sync = true
-- change focus to skim after command `:VimtexView` is given
g.vimtex_view_skim_activate = true

-- don't show log for every compilation.
g.vimtex_quickfix_mode = 0
g.vimtex_log_ignore = {
    "Underfull",
    "Overfull",
}

-- :h vimtex-af-enhanced-matchparen
-- wrong matching for some {} with this on:
-- vim.g.matchup_override_vimtex = true
