local util = require "utils/init"
local map = require "utils/keymap"

local function load_config()
	local cmd = { "kitty", "@", "load-config" }
	vim.system(cmd, { text = true }, function(obj) util.notify_failure(cmd, obj) end)
end

map.buf('n', '<leader>cc', load_config, "Reload kitty config")

local group = vim.api.nvim_create_augroup("kitty_reload", { clear = false })
vim.api.nvim_clear_autocmds({ group = group, buffer = 0 })
vim.api.nvim_create_autocmd("BufWritePost", {
	group = group,
	buffer = 0,
	callback = load_config,
})
