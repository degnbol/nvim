-- Minimal init for plenary tests.
-- Add the lua path so utils modules can be found.
vim.opt.runtimepath:append(vim.fn.getcwd())
-- Add plenary for testing.
vim.opt.runtimepath:append(vim.fn.expand("~/.local/share/nvim/lazy/plenary.nvim"))
