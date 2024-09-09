#!/usr/bin/env lua

-- quickfix specific keymaps.
-- Don't change buffer in quickfix window, first go to other win.
vim.keymap.set('n', '[b', "<C-w><C-w>[b", { desc="", buffer=true, remap=true, nowait=true, silent=true })
vim.keymap.set('n', ']b', "<C-w><C-w>]b", { desc="", buffer=true, remap=true, nowait=true, silent=true })
