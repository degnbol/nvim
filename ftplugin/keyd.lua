-- Reload settings on write.
vim.api.nvim_create_autocmd("BufWritePost", {
	buffer = 0,
	callback = function() vim.fn.system("keyd reload") end,
})
