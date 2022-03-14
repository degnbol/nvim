#!/usr/bin/env lua

-- try out skim, otherwise preview would open.
vim.g.vimtex_view_method = 'skim'
-- don't show log for every compilation.
vim.g.vimtex_quickfix_mode = 0
vim.g.vimtex_log_ignore = {
    "Underfull",
    "Overfull",
}

-- :h vimtex-af-enhanced-matchparen
-- wrong matching for some {} with this on:
-- vim.g.matchup_override_vimtex = true
