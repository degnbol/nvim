-- use the zsh parser for all shell filetypes. tree-sitter-bash is unmaintained
-- and has cascading error bugs (parameter expansion, extglob, heredocs).
-- zsh is a superset of bash, so the zsh parser handles bash files fine.
-- In after/ because nvim-treesitter's plugin/filetypes.lua registers bash for
-- sh, and is sourced again in the plugin-loading step after plugin/.
vim.treesitter.language.register("zsh", "sh")
vim.treesitter.language.register("zsh", "bash")
vim.treesitter.language.register("zsh", "sh.zsh")
