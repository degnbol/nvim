-- Update vim.pack plugins, tree-sitter parsers and Mason packages. Blocks until done.
-- Needs the user config loaded, see update.sh.

vim.pack.update(nil, { force = true })

local ok, err = require("nvim-treesitter").update(nil, { summary = true }):pwait()
assert(ok, err)

require("lz.n").trigger_load("mason.nvim")
local registry = require("mason-registry")
local a = require("mason-core.async")
local updated, result = a.run_blocking(function() return a.wait(registry.update) end)
assert(updated, result)

local outdated = vim.iter(registry.get_installed_packages())
    :filter(function(pkg) return pkg:get_installed_version() ~= pkg:get_latest_version() end)
    :map(function(pkg) return pkg.name end)
    :totable()
-- Blocks when headless and exits 1 if any install fails.
if #outdated > 0 then vim.cmd.MasonInstall(outdated) end
