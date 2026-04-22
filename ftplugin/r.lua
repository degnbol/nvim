vim.cmd([[syntax region @variable matchgroup=Operator start=/`/ end=/`/]])

local map = require "utils/keymap"

-- Errors.
vim.diagnostic.enable(false)

map.buf('n', '<leader>cc', '<Cmd>!Rscript %<CR>', "Run this script")
