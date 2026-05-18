-- Minimal init for plenary tests.
-- Add the lua path so utils modules can be found.
vim.opt.runtimepath:append(vim.fn.getcwd())
-- Vendored plenary — clone with: make test-deps
vim.opt.runtimepath:append("rtps/plenary.nvim")
-- Specs create many [No Name] buffers; swap files collide and trigger E303.
vim.opt.swapfile = false
