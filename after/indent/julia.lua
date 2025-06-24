#!/usr/bin/env lua
-- Undo julia runtime setting a broken julia indentexpr GetJuliaIndent()
vim.bo.indentexpr = "v:lua.require'nvim-treesitter'.indentexpr()"
