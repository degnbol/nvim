" remove line numbers from terminals (REPLs)
au TermOpen * setlocal listchars= nonumber norelativenumber
" close terminal buffer without showing the exit status of the shell
au TermClose * call feedkeys("\<cr>")

" keymaps
" in ther terminal map escape to changing from terminal mode (insert mode) to
" normal terminal mode <C-\><C-n> then change window to the left assuming that term is on the right
tnoremap <Esc> <C-\><C-n><C-w>h
