-- insert real tab by default
vim.opt_local.expandtab = false
-- there was some annoying issue with indents inserted in regular text
vim.opt_local.smartindent = false
-- for text we want to see tabs and spaces by default
vim.opt_local.list = true

-- I copied the default ftplugin code here and rm a line unsetting
-- commentstring
if vim.b.did_ftplugin then return end
vim.b.did_ftplugin = 1

vim.b.undo_ftplugin = "setlocal comments< commentstring<"

-- Pseudo comment leaders to indent bulleted lists with '-' and '*'.  And allow
-- for Mail quoted text with '>'.
vim.opt_local.comments = "fb:-,fb:*,n:>"
