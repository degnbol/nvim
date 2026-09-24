local util = require "utils/init"

-- Reload settings on write.
local group = vim.api.nvim_create_augroup("keyd_reload", { clear = false })
vim.api.nvim_clear_autocmds({ group = group, buf = 0 })
vim.api.nvim_create_autocmd("BufWritePost", {
	group = group,
	buf = 0,
	callback = function()
		local cmd = { "keyd", "reload" }
		vim.system(cmd, { text = true }, function(obj) util.notify_failure(cmd, obj) end)
	end,
})
