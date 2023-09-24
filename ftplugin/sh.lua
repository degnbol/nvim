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
