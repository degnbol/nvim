local opt = vim.opt
local g   = vim.g
local api = vim.api

local rtp = vim.opt.runtimepath:get()[1]

opt.path:append("./src")

g.mapleader = ' '
-- could also use default \ for single key more conveniently allowing for sub categories
g.maplocalleader = '  '

opt.expandtab = true
opt.tabstop = 4        -- how many spaces does a tab correspond to?
opt.shiftwidth = 0     -- use tabstop number of spaces for indentation
opt.smartindent = true
opt.breakindent = true -- when wrapping line, match indent on the wrapped line.
-- opt.breakindentopt = "shift:2" -- indent to show line was wrapped.
-- OR show "> "
-- opt.showbreak = " "
-- opt.showbreak = " "
opt.showbreak = " "
-- opt.showbreak = " "
-- opt.showbreak = " "
-- opt.showbreak = "▒"
-- opt.showbreak = "▋"
-- opt.showbreak = "▉"
opt.copyindent = true
-- indent after words in cinwords (^for,^while,...) and stuff with {}. Should def not be active for normal text docs.
opt.smartindent = false -- has to be set to false explicitly even though it is default probs because some plugin changes it.
-- note vimscript indentexpr
-- https://github.com/JuliaEditorSupport/julia-vim/blob/master/indent/julia.vim
-- is terrible for julia, so definitely use
-- treesitter indent, if for no other lang.
vim.bo.indentexpr = "v:lua.require'nvim-treesitter'.indentexpr()"

-- s: Don't say "search hit BOTTOM, continuing at TOP"
-- A: Don't make ATTENTION warning when a swap file exists, i.e. file is open somewhere else.
vim.opt.shortmess:append("sA")

-- Has to be set for bufferline to work by hiding an open buffer when switching to another
opt.hidden = true
opt.ignorecase = true -- search ignoring case. use \c \C anywhere in search pattern to force case-sensitivity.
opt.smartcase = true  -- only match case-insensitively is query is all lowercase
opt.scrolloff = 4     -- number of lines of context to always keep above and below the cursorline
-- Use default 0 since otherwise sideways scrolling will be stopped in annoying ways.
-- E.g. there is one long line and cursor is on another short line. We then can't scroll along the long line without moving cursor to it first.
-- opt.sidescrolloff = 12 -- number of blocks to keep to sides of cursor
opt.splitbelow = true
opt.splitright = true
opt.wildmode = 'longest:full,full' -- settings for how to show completion on command line
opt.wildcharm = 9                  -- enables cmdline tab completion when recording macro. 9 is the ascii code (and UTF8 code) for tab.
-- opt.number = true -- show line numbering by default. yon toggles
-- opt.relativenumber = true -- should the line numbering be shown relative to current line?
opt.clipboard = 'unnamed,unnamedplus' -- share clipboard between copy paste and yank
opt.wrap = false                      -- something run before init.lua is changing the default so we change it back here.
opt.smoothscroll = true               -- if we wrap lines, then show partial start
-- Whether to break lines when wrapping at whitespace etc. instead of middle of word.
-- It might seem like that would be nice but breaking in the middle of word
-- means only breaking at the end of the screen, so it is much more clear that
-- a linebreak has been forced.
-- opt.linebreak = true
opt.numberwidth = 2 -- reduce default numbering from starting as 3 characters wide to 2
-- Mouse click navigation even in cmdline mode. Default is not in cmdline mode but otherwise.
opt.mouse = "a"
opt.mousescroll = "ver:1,hor:1"
opt.termguicolors = true
opt.cursorline = true        -- highlight current line
opt.cursorlineopt = "number" -- only highlight cursorline number
-- show a column that can be used to add signs to lines showing git changes and LSP diagnostics.
-- "number" means it replaces line numbering rather than e.g. "yes" where it is a column left of numbering.
opt.signcolumn = "no" -- "number"
-- opt.cmdheight = 0 -- hide cmdline when not in use. Messes with search currently, by asking for confirm after a search.
-- http://stackoverflow.com/questions/2490227/how-does-vims-autoread-work#20418591
-- when regaining focus, reload file if it was changed somewhere else
api.nvim_create_autocmd({ "FocusGained", "BufEnter" }, { command = ':silent! !' })
opt.showmode = false
opt.showcmd = false
-- t=use textwidth for formatting. a=auto format. w=respect explicit newline. r=continue comment leader with newline in insert mode.
-- tcqj is default, so only adding w which is relevant when autoformatting with set fo+=a or manually with gq
opt.formatoptions = 'tcqjwr'
-- why is this not default. Persistent undo history.
opt.undofile = true
-- set term title based on file being edited.
opt.title = true
-- ms of wait before keybinding times out, default 1000
-- with 500 I'm sometimes too slow
opt.timeoutlen = 750

-- go between lines with left/right arrow keys only in insert mode
opt.whichwrap = '[,]'

-- Hide the statusline when there's only one file open.
opt.laststatus = 1
-- hide the location in file by default
opt.ruler = false
-- hide ~ tilde at end of buffer.
opt.fillchars = "eob: "
-- add space symbol to whitespace chars
opt.listchars:append('space:⋅')
opt.listchars:append('tab:▏ ')
-- change from showing default ---- to ╱╱╱╱ for deleted lines in git diff
opt.fillchars:append('diff:╱')
-- Danglish support
opt.keymap = "danglish"
-- use ctrl+6 to toggle
opt.iminsert = 0

-- 'set spell' to show red underline for spelling errors.
-- As of writing, spell is on for markdown (ftplugin/markdown.lua)
vim.opt.spelllang = { 'en', 'da' }
-- never complain about sentence starting with lowercase word
vim.opt.spellcapcheck = ""
-- look for spelling in a camelCase word as multiple distinct words
vim.opt.spelloptions = "camel"
-- synonyms <C-xt>
vim.opt.thesaurus = rtp .. "/thesaurus/english.txt"
-- complete word spelling <C-xk>
vim.opt.dictionary = rtp .. "/spell/en.dic"
-- custom words. add under cursor: zg, remove: zw. temp: z{G,W}. undo: zu{g,w,G,W}
-- visual mode also works.
-- It's possible to have multiple spellfiles and use a preceding count.
vim.opt.spellfile = rtp .. "/spell/custom.utf8.add"
-- set a default commentstring
vim.opt.commentstring = "#%s"

vim.wo.foldlevel = 99 -- so we don't fold from the start
-- Fallback if treesitter folding doesn't work:
-- vim.wo.foldmethod = 'indent'
vim.wo.foldmethod = 'expr'
vim.wo.foldexpr = 'v:lua.vim.treesitter.foldexpr()'

-- Reduce content jumping when splitting and unsplitting screens.
-- Seems cleaner when doing e.g. goto-ref, then close qf.
vim.opt.splitkeep = "screen"

-- vim.opt.messagesopt='wait:200,history:500'

-- blinking cursor would be nice but only after jumps
-- vim.opt.guicursor:append("n-v-sm:blinkon150")

-- only show error for virtual_text since it is often incorrect and is distracting.
vim.diagnostic.config {
    -- virtual_text = {severity = vim.diagnostic.severity.ERROR}
}

-- custom foldtext function that shows <first line> … <lines hidden> … <last line>
-- FIXME: when we scroll and are not seeing the beginning of the line the text moves
function FoldText()
    local linestart = vim.fn.getline(vim.v.foldstart)
    local lineend = vim.fn.getline(vim.v.foldend)
    local line_count = vim.v.foldend - vim.v.foldstart + 1
    return linestart .. " … " .. line_count .. " … " .. lineend:match("%s*(.*)%s*")
end

vim.opt.foldtext = "v:lua.FoldText()"
vim.opt.fillchars:append('fold: ')

-- Enable second pass hunk visual that aligns lines for a git diff better giving nicer overview of changes.
-- https://old.reddit.com/r/neovim/comments/1ihpvaf/the_linematch_diffopt_makes_builtin_diff_so_sweat/
vim.opt.diffopt:append("linematch:60")

-- disable coding ligatures for some filetypes that aren't code, e.g. we don't want == ligature when it is used for heading levels.
local kitty = require "utils/kitty"
kitty.ligatures_pattern(false, { "*.typ", "*.adoc" })

-- replace the default lsp diagnostic letters with prettier symbols
vim.fn.sign_define("LspDiagnosticsSignError", { text = "", numhl = "LspDiagnosticsDefaultError" })
vim.fn.sign_define("LspDiagnosticsSignWarning", { text = "", numhl = "LspDiagnosticsDefaultWarning" })
vim.fn.sign_define("LspDiagnosticsSignInformation", { text = "", numhl = "LspDiagnosticsDefaultInformation" })
vim.fn.sign_define("LspDiagnosticsSignHint", { text = "", numhl = "LspDiagnosticsDefaultHint" })
