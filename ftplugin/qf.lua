local map = require "utils/keymap"

map.n('<S-CR>', "<CR><Cmd>:ccl<CR>", "Goto qf entry and close qf", { buf=0 })

-- quickfix specific keymaps.
-- Don't change buffer in quickfix window, first go to other win.
map.n('[b', "<C-w><C-w><Cmd>lua MiniBracketed.buffer('forward')<CR>", "Change buf outside qf", { buf=0, silent=true })
map.n(']b', "<C-w><C-w><Cmd>lua MiniBracketed.buffer('forward')<CR>", "Change buf outside qf", { buf=0, silent=true })

