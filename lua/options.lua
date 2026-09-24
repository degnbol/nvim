local util = require "utils/init"
local o   = vim.o
local g   = vim.g
local api = vim.api

local config = vim.fn.stdpath("config")

vim.opt.path:append("./src")

g.mapleader = ' '
g.maplocalleader = '\\'

-- disable netrw
g.loaded_netrw = 1
g.loaded_netrwPlugin = 1

o.expandtab = true
o.tabstop = 4        -- how many spaces does a tab correspond to?
o.shiftwidth = 0     -- use tabstop number of spaces for indentation
o.breakindent = true -- when wrapping line, match indent on the wrapped line.
-- o.breakindentopt = "shift:2" -- indent to show line was wrapped.
-- OR show "> "
-- o.showbreak = " "
-- o.showbreak = " "
o.showbreak = " "
-- o.showbreak = " "
-- o.showbreak = " "
-- o.showbreak = "▒"
-- o.showbreak = "▋"
-- o.showbreak = "▉"
o.copyindent = true
-- indent after words in cinwords (^for,^while,...) and stuff with {}. Should def not be active for normal text docs.
o.smartindent = false -- set explicitly because some plugin changes it
-- note vimscript indentexpr
-- https://github.com/JuliaEditorSupport/julia-vim/blob/master/indent/julia.vim
-- is terrible for julia, so definitely use
-- treesitter indent, if for no other lang.
vim.o.indentexpr = "v:lua.require'nvim-treesitter'.indentexpr()"

-- s: Don't say "search hit BOTTOM, continuing at TOP"
-- A: Don't make ATTENTION warning when a swap file exists, i.e. file is open somewhere else.
vim.opt.shortmess:append("sA")

o.ignorecase = true -- search ignoring case. use \c \C anywhere in search pattern to force case-sensitivity.
o.smartcase = true  -- only match case-insensitively is query is all lowercase
o.scrolloff = 4     -- number of lines of context to always keep above and below the cursorline
-- Use default 0 since otherwise sideways scrolling will be stopped in annoying ways.
-- E.g. there is one long line and cursor is on another short line. We then can't scroll along the long line without moving cursor to it first.
-- o.sidescrolloff = 12 -- number of blocks to keep to sides of cursor
o.splitbelow = true
o.splitright = true
o.wildmode = 'longest:full,full' -- settings for how to show completion on command line
o.wildcharm = 9                  -- enables cmdline tab completion when recording macro. 9 is the ascii code (and UTF8 code) for tab.
-- o.number = true -- show line numbering by default. yon toggles
-- o.relativenumber = true -- should the line numbering be shown relative to current line?
o.clipboard = 'unnamed,unnamedplus' -- share clipboard between copy paste and yank
o.wrap = false                      -- something run before init.lua is changing the default so we change it back here.
o.smoothscroll = true               -- if we wrap lines, then show partial start
-- Whether to break lines when wrapping at whitespace etc. instead of middle of word.
-- It might seem like that would be nice but breaking in the middle of word
-- means only breaking at the end of the screen, so it is much more clear that
-- a linebreak has been forced.
-- o.linebreak = true
o.numberwidth = 2 -- reduce default numbering from starting as 3 characters wide to 2
-- Touchpad on Mac should scroll slower.
if util.is_mac() then
    o.mousescroll = "ver:1,hor:1"
end
o.cursorline = true        -- highlight current line
o.cursorlineopt = "number" -- only highlight cursorline number
-- show a column that can be used to add signs to lines showing git changes and LSP diagnostics.
-- "number" means it replaces line numbering rather than e.g. "yes" where it is a column left of numbering.
o.signcolumn = "no" -- "number"
-- o.cmdheight = 0 -- hide cmdline when not in use. Messes with search currently, by asking for confirm after a search.
-- when regaining focus, check if files changed on disk and reload silently.
-- Only on FocusGained (returning from another app), not BufEnter — checktime
-- on BufEnter causes recursive autocmd cascades with plugins that do
-- programmatic buffer switches (agentic, fzf, etc.): checktime → reload →
-- BufReadPost → buffer switch → BufEnter → checktime, bouncing between buffers.
-- nested = true: checktime triggers FileChangedShell, which needs to fire
-- from within this handler. Without nested, neovim's autocmd anti-nesting
-- guard (autocmd_busy && !autocmd_nested) silently skips FileChangedShell,
-- causing W12 to display instead of silent reload. With ui2, the W12 prompt
-- crashed neovim when concurrent async events fired (fixed in 0.12.1).
api.nvim_create_autocmd("FocusGained", { nested = true, command = 'silent! checktime' })
-- Force-reload buffers when files change on disk (e.g. agent edits).
-- Suppresses W12 prompt so modified buffers reload silently.
api.nvim_create_autocmd("FileChangedShell", {
    pattern = "*",
    callback = function()
        vim.v.fcs_choice = "reload"
    end,
})
o.showmode = false
o.showcmd = false
-- t=use textwidth for formatting. a=auto format. w=respect explicit newline. r=continue comment leader with newline in insert mode.
-- tcqj is default, so only adding w which is relevant when autoformatting with set fo+=a or manually with gq
o.formatoptions = 'tcqjwr'
-- why is this not default. Persistent undo history.
o.undofile = true
-- set term title based on file being edited.
o.title = true
-- ms of wait before keybinding times out, default 1000
-- with 500 I'm sometimes too slow
o.timeoutlen = 750

-- go between lines with left/right arrow keys only in insert mode
o.whichwrap = '[,]'

-- Hide the statusline when there's only one file open.
o.laststatus = 1
-- hide the location in file by default
o.ruler = false
-- hide ~ tilde at end of buffer.
o.fillchars = "eob: ,vert:│,horiz:─,horizup:┴,horizdown:┬,vertleft:┤,vertright:├,verthoriz:┼"
-- add space symbol to whitespace chars
vim.opt.listchars:append('space:⋅')
vim.opt.listchars:append('tab:▏ ')
-- change from showing default ---- to ╱╱╱╱ for deleted lines in git diff
vim.opt.fillchars:append('diff:╱')
-- Danglish support
o.keymap = "danglish"
-- use ctrl+6 to toggle
o.iminsert = 0

-- 'set spell' to show red underline for spelling errors.
-- As of writing, spell is on for markdown (ftplugin/markdown.lua)
vim.opt.spelllang = { 'en', 'da' }
-- never complain about sentence starting with lowercase word
vim.o.spellcapcheck = ""
-- look for spelling in a camelCase word as multiple distinct words
vim.o.spelloptions = "camel"
-- synonyms <C-xt>
vim.o.thesaurus = config .. "/thesaurus/english.txt"
-- complete word spelling <C-xk>
vim.o.dictionary = config .. "/spell/en.dic"
-- custom words. add under cursor: zg, remove: zw. temp: z{G,W}. undo: zu{g,w,G,W}
-- visual mode also works.
-- It's possible to have multiple spellfiles and use a preceding count.
vim.o.spellfile = config .. "/spell/custom.utf8.add"
-- set a default commentstring
vim.o.commentstring = "#%s"

o.foldlevel = 99 -- don't close folds from the start
o.foldminlines = 10 -- don't close trivial folds (automatically)
o.foldmethod = 'expr'
o.foldexpr = 'v:lua.vim.treesitter.foldexpr()'

-- Setting it to "screen" reduces content jumping when splitting and unsplitting screens.
-- Seems cleaner when doing e.g. goto-ref, then close qf.
-- However it's just not as practical. If I split the cursor might be out of view.
vim.o.splitkeep = "cursor"

-- vim.o.messagesopt = 'wait:200,history:500'

-- blinking cursor would be nice but only after jumps
-- vim.opt.guicursor:append("n-v-sm:blinkon150")
-- bar cursor in terminal mode (e.g. agent window)
vim.o.guicursor = "n-v-c-sm:block,i-ci-ve:ver25,r-cr-o:hor20,t:ver25-TermCursor"

-- only show error for virtual_text since it is often incorrect and is distracting.
vim.diagnostic.config {
    -- virtual_text = {severity = vim.diagnostic.severity.ERROR}
    severity_sort = true, -- higher-severity extmark wins on overlap (ERROR underline over WARN)
    signs = {
        text = {
            [vim.diagnostic.severity.ERROR] = "",
            [vim.diagnostic.severity.WARN]  = "",
            [vim.diagnostic.severity.INFO]  = "",
            [vim.diagnostic.severity.HINT]  = "",
        },
        numhl = {
            [vim.diagnostic.severity.ERROR] = "DiagnosticError",
            [vim.diagnostic.severity.WARN]  = "DiagnosticWarn",
            [vim.diagnostic.severity.INFO]  = "DiagnosticInfo",
            [vim.diagnostic.severity.HINT]  = "DiagnosticHint",
        },
    },
}

-- <first line> … <lines hidden> … <last line>
function FoldText()
    local linestart = vim.api.nvim_buf_get_lines(0, vim.v.foldstart - 1, vim.v.foldstart, false)[1]
    local lineend = vim.api.nvim_buf_get_lines(0, vim.v.foldend - 1, vim.v.foldend, false)[1]
    local line_count = vim.v.foldend - vim.v.foldstart + 1
    return linestart .. " … " .. line_count .. " … " .. lineend:match("%s*(.*)%s*")
end

vim.o.foldtext = "v:lua.FoldText()"
vim.opt.fillchars:append('fold: ')

-- Enable second pass hunk visual that aligns lines for a git diff better giving nicer overview of changes.
-- https://old.reddit.com/r/neovim/comments/1ihpvaf/the_linematch_diffopt_makes_builtin_diff_so_sweat/
-- Only :set+= merges key:value items; vim.opt:append would add a second linematch.
vim.cmd("set diffopt+=linematch:60")

-- disable coding ligatures for some filetypes that aren't code, e.g. we don't want == ligature when it is used for heading levels.
local kitty = require "utils/kitty"
kitty.ligatures_pattern(false, { "*.typ", "*.adoc" })

-- Enable experimental cmdline with syntax hl (eliminates "Press ENTER" prompts).
require('vim._core.ui2').enable {}

