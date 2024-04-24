-- textobj iP and aP is for between pipes in e.g. zsh by default but also map 
-- i| and a|. Note that textobj "|" means column down, which is useful for e.g. 
-- <C-v>|<C-q> where <C-q> is my custom blockim keymap.
vim.keymap.set({ "o", "x" }, "i|",
	"<cmd>lua require('various-textobjs').shellPipe(true)<CR>",
	{ buffer = true }
)
vim.keymap.set({ "o", "x" }, "a|",
	"<cmd>lua require('various-textobjs').shellPipe(false)<CR>",
	{ buffer = true }
)

vim.opt_local.list = false

-- was linked to Conditional which is italic.
-- It matches == and ! in a conditional, which should not be italic.
vim.api.nvim_set_hl(0, "shTestOpr", {link="Operator", default=true})

-- stop continuing comments with o/O
vim.opt_local.formatoptions:remove('o')
