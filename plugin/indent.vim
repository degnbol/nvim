" indentation can also be done in insert mode with <C+D>, <C+T> and <C+F> (see :h i_CTRL-T etc.)

" copy pasted from https://vi.stackexchange.com/questions/18310/keep-relative-cursor-position-after-indenting-with
" modified to handle vim-repeat based on http://vimcasts.org/episodes/creating-repeatable-mappings-with-repeat-vim/

" indent cursor according to shiftwidth
func! IndentInc()
    exe "norm!" . (virtcol('.') + shiftwidth()) . '|'
endfunc
func! IndentDcr()
    exe "norm!" . (virtcol('.') - shiftwidth()) . '|'
endfunc

" create a <Plug> command and setting it as last edit with vim-repeat#set
nnoremap <silent> <Plug>IndentInc >>:call IndentInc()<CR>:call repeat#set("\<Plug>IndentInc")<CR>
nnoremap <silent> <Plug>IndentDcr <<:call IndentDcr()<CR>:call repeat#set("\<Plug>IndentDcr")<CR>
" use the regular >> and << to call this modified version
nmap >> <Plug>IndentInc
nmap << <Plug>IndentDcr

