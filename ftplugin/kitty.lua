local map = require "utils/keymap"

map.buf('n', '<leader>cc', function() vim.fn.system("kitty @ load-config") end, "Reload kitty config")

vim.api.nvim_create_autocmd("BufWritePost", {
	buffer = 0,
	callback = function() vim.fn.system("kitty @ load-config") end,
})
