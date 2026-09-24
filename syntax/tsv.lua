vim.cmd [[syntax match Comment "^#.*$"]]

-- Load the runtime csv.vim with tab delimiter (same as runtime tsv.vim would),
-- then set b:current_syntax to prevent the runtime tsv.vim from re-sourcing it.
vim.b.csv_delimiter = "\t"
vim.cmd.runtime { "syntax/csv.vim", bang = true }

-- Clear rainbow column highlights set by runtime csv.vim.
-- The colorscheme generator also clears these at load time, but csv.vim re-applies
-- them via `hi def link` each time a TSV buffer's syntax is loaded.
for i = 1, 11 do
    vim.api.nvim_set_hl(0, "csvCol" .. i, {})
end

vim.b.current_syntax = "tsv"
