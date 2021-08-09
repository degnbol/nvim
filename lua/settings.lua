local utils = require('utils')
local cmd = vim.cmd
local opt = vim.opt

local indent = 4

cmd 'syntax enable'
cmd 'filetype plugin indent on'
utils.opt('b', 'expandtab', true)
cmd 'autocmd BufRead,BufNewFile *.tsv,*.tab,*.txt setlocal noexpandtab'  -- Only for text files insert real tab by default.
utils.opt('b', 'shiftwidth', indent)
utils.opt('b', 'smartindent', true)
utils.opt('b', 'tabstop', indent) -- how many spaces does a tab correspond to?
utils.opt('o', 'hidden', true)
utils.opt('o', 'ignorecase', true)
utils.opt('o', 'scrolloff', 4)
utils.opt('o', 'shiftround', true)  -- round indent
utils.opt('o', 'smartcase', true)
utils.opt('o', 'splitbelow', true)
utils.opt('o', 'splitright', true)
utils.opt('o', 'wildmode', 'list:longest')
utils.opt('w', 'number', true) -- show line numbering
utils.opt('w', 'relativenumber', false) -- should the line numbering be shown relative to current line?
utils.opt('o', 'clipboard','unnamed,unnamedplus') -- share clipboard between copy paste and yank
opt.wrap = false -- don't wrap lines
utils.opt("o", "numberwidth", 2) -- reduce default numbering from starting as 3 characters wide to 2
utils.opt("o", "mouse", "a") -- activate the mouse, i.e. click to move cursor, drag to visual select and scroll to scroll window instead of cursor
utils.opt("o", "termguicolors", true)
-- utils.opt("w", "cul", true) -- highlight current line
utils.opt("w", "signcolumn", "yes") -- show a column left of the numbering that can be used to add signs to lines for git etc.
utils.opt("o", "updatetime", 2000) -- update interval can be decreased for gitsigns, default = 4000
-- utils.opt("o", "timeoutlen", 500) -- time in ms to wait for mapped sequence to complete. Default = 1000. Reduced for quicker suggestion pop-up

-- Highlight on yank, e.g. press Y to yank line which will highlight the line for a moment
cmd 'au TextYankPost * lua vim.highlight.on_yank {on_visual = false}'


-- yoink integration with cutlass
vim.g.yoinkIncludeDeleteOperations = 1
-- add yanks to numbered register
vim.g.yoinkSyncNumberedRegisters = 1 
-- move cursor to end instead of start of multi-line paste
vim.g.yoinkMoveCursorToEndOfPaste = 1

