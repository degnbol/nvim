local utils = require('utils')
local cmd = vim.cmd
local opt = vim.opt


cmd 'syntax on'
utils.opt('b', 'expandtab', true)
cmd 'autocmd BufRead,BufNewFile *.tsv,*.tab,*.txt setlocal noexpandtab'  -- Only for text files insert real tab by default.
utils.opt('b', 'shiftwidth', 4)
utils.opt('b', 'smartindent', true)
utils.opt('b', 'copyindent', true)
utils.opt('b', 'smartindent', true) -- round off indents to multiple of indent size (4)
utils.opt('b', 'tabstop', 4) -- how many spaces does a tab correspond to?
utils.opt('o', 'hidden', true)
utils.opt('o', 'ignorecase', true) -- search ignoring case. use \c \C anywhere in search pattern to force case-sensitivity.
utils.opt('o', 'smartcase', true) -- only match case-insensitively is query is all lowercase
utils.opt('o', 'scrolloff', 4) -- number of lines of context to always keep above and below the cursorline
utils.opt('o', 'splitbelow', true)
utils.opt('o', 'splitright', true)
utils.opt('o', 'wildmode', 'longest:full,full') -- settings for how to show completion on command line
utils.opt('w', 'number', true) -- show line numbering
-- utils.opt('w', 'relativenumber', true) -- should the line numbering be shown relative to current line?
utils.opt('o', 'clipboard','unnamed,unnamedplus') -- share clipboard between copy paste and yank
opt.wrap = false -- don't wrap lines
utils.opt("o", "numberwidth", 2) -- reduce default numbering from starting as 3 characters wide to 2
utils.opt("o", "mouse", "a") -- activate the mouse, i.e. click to move cursor, drag to visual select and scroll to scroll window instead of cursor
utils.opt("o", "termguicolors", true)
-- utils.opt("w", "cul", true) -- highlight current line
utils.opt("w", "signcolumn", "yes") -- show a column left of the numbering that can be used to add signs to lines for git etc.
-- http://stackoverflow.com/questions/2490227/how-does-vims-autoread-work#20418591
-- when regaining focus, reload file if it was changed somewhere else
cmd 'autocmd FocusGained,BufEnter * :silent! !'


-- Highlight on yank, e.g. press Y to yank line which will highlight the line for a moment
cmd 'autocmd TextYankPost * lua vim.highlight.on_yank {on_visual = false}'
-- yoink integration with cutlass
vim.g.yoinkIncludeDeleteOperations = 1
-- add yanks to numbered register
vim.g.yoinkSyncNumberedRegisters = 1 
-- move cursor to end instead of start of multi-line paste
vim.g.yoinkMoveCursorToEndOfPaste = 1

