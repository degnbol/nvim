local utils = require('utils')
local cmd = vim.cmd
local opt = vim.opt
local g = vim.g

opt.syntax = "enable"
opt.expandtab = true
cmd 'autocmd BufRead,BufNewFile *.tsv,*.tab,*.txt setlocal noexpandtab'  -- Only for text files insert real tab by default.
opt.shiftwidth = 4
opt.smartindent = true
opt.copyindent = true
opt.smartindent = true -- round off indents to multiple of indent size (4)
opt.tabstop = 4 -- how many spaces does a tab correspond to?
opt.hidden = true
opt.ignorecase = true -- search ignoring case. use \c \C anywhere in search pattern to force case-sensitivity.
opt.smartcase = true -- only match case-insensitively is query is all lowercase
opt.scrolloff = 4 -- number of lines of context to always keep above and below the cursorline
opt.splitbelow = true
opt.splitright = true
opt.wildmode = 'longest:full,full' -- settings for how to show completion on command line
-- opt.number = true -- show line numbering
-- opt.relativenumber true -- should the line numbering be shown relative to current line?
opt.clipboard = 'unnamed,unnamedplus' -- share clipboard between copy paste and yank
opt.wrap = false -- don't wrap lines
opt.numberwidth = 2 -- reduce default numbering from starting as 3 characters wide to 2
opt.mouse = "a" -- activate the mouse, i.e. click to move cursor, drag to visual select and scroll to scroll window instead of cursor
opt.termguicolors = true
-- opt.cul = true -- highlight current line
opt.signcolumn = "yes" -- show a column left of the numbering that can be used to add signs to lines for git etc.
-- http://stackoverflow.com/questions/2490227/how-does-vims-autoread-work#20418591
-- when regaining focus, reload file if it was changed somewhere else
cmd 'autocmd FocusGained,BufEnter * :silent! !'
opt.showmode = false
opt.showcmd = false

-- Highlight on yank, e.g. press Y to yank line which will highlight the line for a moment
cmd 'autocmd TextYankPost * lua vim.highlight.on_yank {on_visual = false}'
-- yoink integration with cutlass
g.yoinkIncludeDeleteOperations = 1
-- add yanks to numbered register
g.yoinkSyncNumberedRegisters = 1 
-- move cursor to end instead of start of multi-line paste
g.yoinkMoveCursorToEndOfPaste = 1

