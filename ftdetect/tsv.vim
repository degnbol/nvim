" create a neovim filetype called tsv which is recognizes by file endings .tab and .tsv
au BufRead,BufNewFile *.tsv,*.tab,*.bed setfiletype tsv
