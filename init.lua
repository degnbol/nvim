-- Check for large files before plugins load
require("largefile").check_argv()

-- Leaders must be set before plugin loading
require "options"

-- Dev/local plugins: add to rtp before anything else
local cfg = vim.fn.stdpath("config")
for mod, kind in vim.fs.dir(cfg .. "/modules") do
    if kind == "directory" then
        local base = cfg .. "/modules/" .. mod
        vim.opt.runtimepath:append(base)
        local after_dir = base .. "/after"
        if vim.uv.fs_stat(after_dir) then vim.opt.runtimepath:append(after_dir) end
        local doc_dir = base .. "/doc"
        if vim.uv.fs_stat(doc_dir) then vim.cmd.helptags(doc_dir) end
    end
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

