require "guitar.vnewFurther"
require "guitar.chords"
require "guitar.tabs"
local map = require "utils/keymap"

-- smarttab and autoformat are annoying when writing the chord lines that have
-- lots of empty spaces.
vim.opt_local.smarttab = false
vim.opt_local.formatoptions:remove('a')

map.buf('n', '<localleader>v', vnewFurther,     "vnew with one window height offset and scrollbind enabled. Accepts a count.")
map.buf('n', '<localleader>f', forwardScreen,   "Like <C-f> but jumps foward one screen visually in the split regardless of scrolloff.")
map.buf('n', '<localleader>b', backwardScreen,  "Like <C-b> but jumps backwards one screen visually in the split regardless of scrolloff.")
map.buf('n', '<localleader>c', prettyChordLine, "Write a pretty print of chord at current line.")
map.buf('n', '<localleader>t', tabNew,          "New tab")
-- map.buf('i', '-',             tabDash,        "Potentially extend tab")
-- sometimes you want to simply BS normally, e.g. while writing anything other than additional -.
-- the <leader><C-v> mapping is useful for deleting parts of tab.
-- map.buf('i', '<BS>',          tabBS,          "Potentially delete tab")
map.buf('n', '<leader><C-V>',   tabSelect,      "Select all 6 strings in tab", { expr = true })

-- TODO: disable completion suggestions inside syntax "Tab"
vim.api.nvim_set_hl(0, "Tab", {link="String", default=true})

