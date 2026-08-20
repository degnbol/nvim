-- Check for large files before plugins load
require("largefile").check_argv()

-- Leaders must be set before plugin loading
require "options"

-- Dev/local plugins: add to rtp before anything else
for _, base in ipairs(require("dev_plugins").add(vim.fn.stdpath("config"))) do
    local doc_dir = base .. "/doc"
    if vim.uv.fs_stat(doc_dir) then vim.cmd.helptags(doc_dir) end
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

require("path_highlight").setup()
