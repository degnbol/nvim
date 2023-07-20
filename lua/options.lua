local cmd = vim.cmd
local opt = vim.opt
local g = vim.g
local api = vim.api

local rtp = vim.opt.runtimepath:get()[1]

g.mapleader = ' '
g.maplocalleader = ' '

opt.expandtab = true
opt.shiftwidth = 4
opt.smartindent = true
opt.breakindent = true -- when wrapping line, match indent on the wrapped line.
-- opt.breakindentopt = "shift:2" -- indent to show line was wrapped.
-- OR show "> "
opt.showbreak = "> "
opt.copyindent = true
-- indent after words in cinwords (^for,^while,...) and stuff with {}. Should def not be active for normal text docs.
opt.smartindent = false -- has to be set to false explicitly even though it is default probs because some plugin changes it.
opt.tabstop = 4 -- how many spaces does a tab correspond to?
-- Has to be set for bufferline to work by hiding an open buffer when switching to another
opt.hidden = true 
opt.ignorecase = true -- search ignoring case. use \c \C anywhere in search pattern to force case-sensitivity.
opt.smartcase = true -- only match case-insensitively is query is all lowercase
opt.scrolloff = 4 -- number of lines of context to always keep above and below the cursorline
opt.sidescrolloff = 12 -- number of blocks to keep to sides of cursor
opt.splitbelow = true
opt.splitright = true
opt.wildmode = 'longest:full,full' -- settings for how to show completion on command line
-- opt.number = true -- show line numbering by default. yon toggles
-- opt.relativenumber = true -- should the line numbering be shown relative to current line?
opt.clipboard = 'unnamed,unnamedplus' -- share clipboard between copy paste and yank
opt.wrap = false -- something run before init.lua is changing the default so we change it back here.
opt.smoothscroll = true -- if we wrap lines, then show partial start
opt.linebreak = true -- if I sometimes were to wrap lines, do it at whitespaces and other characters indicated in breakat
opt.numberwidth = 2 -- reduce default numbering from starting as 3 characters wide to 2
opt.mouse = "a" -- activate the mouse, i.e. click to move cursor, drag to visual select and scroll to scroll window instead of cursor
opt.termguicolors = true
opt.cursorline = true -- highlight current line
opt.cursorlineopt = "number" -- only highlight cursorline number
-- show a column that can be used to add signs to lines showing git changes and LSP diagnostics.
-- "number" means it replaces line numbering rather than e.g. "yes" where it is a column left of numbering.
opt.signcolumn = "no" -- "number"
-- opt.cmdheight = 0 -- hide cmdline when not in use. Messes with search currently, by asking for confirm after a search.
-- http://stackoverflow.com/questions/2490227/how-does-vims-autoread-work#20418591
-- when regaining focus, reload file if it was changed somewhere else
api.nvim_create_autocmd({"FocusGained", "BufEnter"}, { command=':silent! !' })
opt.showmode = false
opt.showcmd = false
-- t=use textwidth for formatting. a=auto format. w=respect explicit newline. r=continue comment leader with newline in insert mode.
-- tcqj is default, so only adding w which is relevant when autoformatting with set fo+=a or manually with gq
opt.formatoptions='tcqjwr'
-- why is this not default. Persistent undo history.
opt.undofile = true
-- set term title based on file being edited.
opt.title = true
-- ms of wait before keybinding times out, default 1000
-- with 500 I'm sometimes too slow
opt.timeoutlen = 750

-- yoink integration with cutlass
g.yoinkIncludeDeleteOperations = 1
-- add yanks to numbered register
g.yoinkSyncNumberedRegisters = 1 
-- move cursor to end instead of start of multi-line paste
g.yoinkMoveCursorToEndOfPaste = 1
-- preserve yank between neovim sessions
g.yoinkSavePersistently = 1

-- go between lines with left/right arrow keys only in insert mode
opt.whichwrap = '[,]'

-- hide the statusline that shows current file being edited
opt.laststatus = 0
-- hide the location in file by default
opt.ruler = false
-- hide ~ tilde at end of buffer.
opt.fillchars = "eob: "
-- add space symbol to whitespace chars
opt.listchars:append( 'space:⋅')
opt.listchars:append( 'tab:▏ ' )

-- Danglish support
opt.keymap = "danglish"
-- use ctrl+6 to toggle
opt.iminsert = 0

-- 'set spell' to show red underline for spelling errors.
-- As of writing, spell is on for markdown (ftplugin/markdown.lua)
vim.opt.spelllang = { 'en', 'da' }
-- never complain about sentence starting with lowercase word
vim.opt.spellcapcheck = ""
-- synonyms <C-xt>
vim.opt.thesaurus = rtp .. "/thesaurus/english.txt"
-- complete word spelling <C-xk>
vim.opt.dictionary = rtp .. "/spell/en.dic"
-- custom words. add under cursor: zg, remove: zw. temp: z{G,W}. undo: zu{g,w,G,W}
-- visual mode also works.
-- It's possible to have multiple spellfiles and use a preceding count.
vim.opt.spellfile = rtp .. "/spell/custom.utf8.add"


