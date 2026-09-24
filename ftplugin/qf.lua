local map = require "utils/keymap"

-- Fit to contents. The height from `:copen [height]` stays the upper bound.
vim.api.nvim_win_set_height(0, math.min(vim.api.nvim_win_get_height(0), vim.api.nvim_buf_line_count(0)))

map.n('<S-CR>', "<CR><Cmd>:ccl<CR>", "Goto qf entry and close qf", { buf=0 })

-- quickfix specific keymaps.
-- Don't change buffer in quickfix window, first go to other win.
map.n('[b', "<C-w><C-w><Cmd>lua MiniBracketed.buffer('forward')<CR>", "Change buf outside qf", { buf=0, silent=true })
map.n(']b', "<C-w><C-w><Cmd>lua MiniBracketed.buffer('forward')<CR>", "Change buf outside qf", { buf=0, silent=true })

