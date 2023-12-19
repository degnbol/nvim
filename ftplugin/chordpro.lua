require "guitar.vnewFurther"
require "guitar.chords"
require "guitar.tabs"

-- smarttab and autoformat are annoying when writing the chord lines that have 
-- lots of empty spaces.
vim.opt_local.smarttab = false
vim.opt_local.formatoptions:remove('a')

vim.keymap.set('n', '<leader><leader>v', vnewFurther,     { buffer=true, desc="vnew with one window height offset and scrollbind enabled. Accepts a count." })
vim.keymap.set('n', '<leader><leader>f', forwardScreen,   { buffer=true, desc="Like <C-f> but jumps foward one screen visually in the split regardless of scrolloff." })
vim.keymap.set('n', '<leader><leader>b', backwardScreen,  { buffer=true, desc="Like <C-b> but jumps backwards one screen visually in the split regardless of scrolloff." })
vim.keymap.set('n', '<leader><leader>c', prettyChordLine, { buffer=true, desc="Write a pretty print of chord at current line." })
vim.keymap.set('n', '<leader><leader>t', tabNew,          { buffer=true, desc="New tab" })
vim.keymap.set('i', '-',                 tabDash,         { buffer=true, desc="Potentially extend tab" })
-- sometimes you want to simply BS normally, e.g. while writing anything other than additional -.
-- the <leader><C-v> mapping is useful for deleting parts of tab.
-- vim.keymap.set('i', '<BS>',              tabBS,           { buffer=true, desc="Potentially delete tab" })
vim.keymap.set('n', '<leader><C-V>',     tabSelect,       { buffer=true, expr=true, desc="Select all 6 strings in tab" })

-- TODO: disable completion suggestions inside syntax "Tab"
vim.api.nvim_set_hl(0, "Tab", {link="String", default=true})

