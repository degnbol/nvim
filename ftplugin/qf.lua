local map = require "utils/keymap"

map.n('<S-CR>', "<CR><Cmd>:ccl<CR>", "Goto qf entry and close qf", { buffer=true })

-- quickfix specific keymaps.
-- Don't change buffer in quickfix window, first go to other win.
map.n('[b', "<C-w><C-w><Cmd>lua MiniBracketed.buffer('forward')<CR>", "Change buf outside qf", { buffer=true, silent=true })
map.n(']b', "<C-w><C-w><Cmd>lua MiniBracketed.buffer('forward')<CR>", "Change buf outside qf", { buffer=true, silent=true })

