-- Minimal init for plenary tests.
-- Add the lua path so utils modules can be found.
vim.opt.runtimepath:append(vim.fn.getcwd())
-- Local plugins carry the grammars' queries; init.lua puts them on rtp too.
require("dev_plugins").add(vim.fn.getcwd())
-- Vendored plenary — clone with: make test-deps
vim.opt.runtimepath:append("rtps/plenary.nvim")
-- Specs create many [No Name] buffers; swap files collide and trigger E303.
vim.opt.swapfile = false
-- PlenaryBustedDirectory runs one nvim per spec file, each with this init, and
-- each would otherwise write the real shada: fixture paths into the file marks
-- and history, plus a leftover main.shada.tmp.X whenever a spec aborts. Exhaust
-- those tmp names and an interactive nvim can no longer write its own (E138).
vim.o.shadafile = "NONE"
-- Specs that cache (pkg-config index, generated compile databases) start cold and
-- leave the real cache alone.
vim.env.XDG_CACHE_HOME = vim.fn.tempname()
