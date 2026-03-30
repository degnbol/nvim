-- Check for large files before plugins load
require("largefile").check_argv()

-- Leaders must be set before plugin loading
require "options"

-- Dev/local plugins: add to rtp before anything else
local cfg = vim.fn.stdpath("config")
for _, mod in ipairs({ "agentic.nvim", "kitty-conf.nvim", "kittyREPL.nvim", "nvim-revJ.lua" }) do
    vim.opt.runtimepath:append(cfg .. "/modules/" .. mod)
    local after_dir = cfg .. "/modules/" .. mod .. "/after"
    if vim.uv.fs_stat(after_dir) then vim.opt.runtimepath:append(after_dir) end
end

-- Build hooks (register BEFORE vim.pack.add so install hooks fire)
vim.api.nvim_create_autocmd("PackChanged", {
    callback = function(ev) require("pack_hooks").on_changed(ev) end,
})

-- Install all remote plugins (opt/ only, nothing loaded yet)
require("pack_specs")

-- Bootstrap lz.n (installed by vim.pack, packadd makes it requireable)
vim.cmd.packadd("lz.n")
require("lz.n").load("plugins")

