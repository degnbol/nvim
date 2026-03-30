syntax match Comment "^#.*$"

" Load the runtime csv.vim with tab delimiter (same as runtime tsv.vim would),
" then set b:current_syntax to prevent the runtime tsv.vim from re-sourcing it.
let b:csv_delimiter = '\t'
runtime! syntax/csv.vim
let b:current_syntax = 'tsv'
