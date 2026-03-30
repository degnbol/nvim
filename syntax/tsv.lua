-- Clear rainbow column highlights set by runtime csv.vim (sourced by runtime tsv.vim).
-- The colorscheme generator also clears these at load time, but csv.vim re-applies
-- them via `hi def link` each time a TSV buffer's syntax is loaded.
for i = 1, 11 do
    vim.api.nvim_set_hl(0, "csvCol" .. i, {})
end
