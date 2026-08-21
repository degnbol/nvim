-- https://github.com/ribru17/ts_query_ls
local function LSP_start(cmdpath)
	-- :InspectTree scratch buffers are filetype=query (vim/treesitter/dev.lua).
	if vim.bo.buftype == "nofile" then
		return
	end
	vim.lsp.start({
		name = "ts_query_ls",
		cmd = { cmdpath },
		root_dir = vim.fs.root(0, { "queries" }),
		settings = {
			-- Resolved off the rtp so it survives a change of plugin manager,
			-- and so first-hit order matches nvim's own parser resolution.
			parser_install_directories = vim.api.nvim_get_runtime_file("parser", true),
			parser_aliases = {
				ecma = "javascript",
			},
			language_retrieval_patterns = {
				"languages/src/([^/]+)/[^/]+\\.scm$",
			},
		},
	})
end

local lsppath = vim.fn.stdpath("config") .. "/lsp_ext/ts_query_ls"
if vim.fn.executable(lsppath) == 1 then
	LSP_start(lsppath)
end
