
cnoreabbrev Gs Git status -uno
cnoreabbrev Gp Git push
cnoreabbrev Gl Git pull
cnoreabbrev Ga Git add %
cnoreabbrev Gc Git commit -m

" substitute motion with clipboard content with https://github.com/svermeulen/vim-subversive
nmap s <plug>(SubversiveSubstitute)
nmap ss <plug>(SubversiveSubstituteLine)
nmap S <plug>(SubversiveSubstituteToEndOfLine)


" yank, cut, delete behavior
" a logical choice that was supposed to be default? https://github.com/neovim/neovim/issues/416
map Y y$
" back space and delete buttons delete should delete also in normal mode.
" Delete to black hole register.
nnoremap <BS> "_X
" cutlass map so x is cut
nnoremap x d
xnoremap x d
nnoremap xx dd
nnoremap X D
" yoink keymaps
nmap <c-n> <plug>(YoinkPostPasteSwapBack)
nmap <c-p> <plug>(YoinkPostPasteSwapForward)
nmap p <plug>(YoinkPaste_p)
nmap P <plug>(YoinkPaste_P)
nmap gp <plug>(YoinkPaste_gp)
nmap gP <plug>(YoinkPaste_gP)
nmap [y <plug>(YoinkRotateBack)
nmap ]y <plug>(YoinkRotateForward)
nmap <c-=> <plug>(YoinkPostPasteToggleFormat)
nmap y <plug>(YoinkYankPreserveCursorPosition)
xmap y <plug>(YoinkYankPreserveCursorPosition)

" one scroll signal is one line change instead of 3
map <ScrollWheelUp>   <c-y>
map <2-ScrollWheelUp> <c-y>
map <3-ScrollWheelUp> <c-y>
map <4-ScrollWheelUp> <c-y>
map <ScrollWheelDown>   <c-e>
map <2-ScrollWheelDown> <c-e>
map <3-ScrollWheelDown> <c-e>
map <4-ScrollWheelDown> <c-e>
" insert mode requires going to normal mode temporarily
imap <ScrollWheelUp>   <c-\><c-o><c-y>
imap <2-ScrollWheelUp> <c-\><c-o><c-y>
imap <3-ScrollWheelUp> <c-\><c-o><c-y>
imap <4-ScrollWheelUp> <c-\><c-o><c-y>
imap <ScrollWheelDown>   <c-\><c-o><c-e>
imap <2-ScrollWheelDown> <c-\><c-o><c-e>
imap <3-ScrollWheelDown> <c-\><c-o><c-e>
imap <4-ScrollWheelDown> <c-\><c-o><c-e>

" shift scroll scrolls horizontally
noremap <S-ScrollWheelUp>   zl
noremap <S-2-ScrollWheelUp> zl
noremap <S-3-ScrollWheelUp> zl
noremap <S-4-ScrollWheelUp> zl
noremap <S-ScrollWheelDown>   zh
noremap <S-2-ScrollWheelDown> zh
noremap <S-3-ScrollWheelDown> zh
noremap <S-4-ScrollWheelDown> zh
" insert mode requires going to normal mode temporarily
inoremap <S-ScrollWheelUp>   <c-\><c-o>zl
inoremap <S-2-ScrollWheelUp> <c-\><c-o>zl
inoremap <S-3-ScrollWheelUp> <c-\><c-o>zl
inoremap <S-4-ScrollWheelUp> <c-\><c-o>zl
inoremap <S-ScrollWheelDown>   <c-\><c-o>zh
inoremap <S-2-ScrollWheelDown> <c-\><c-o>zh
inoremap <S-3-ScrollWheelDown> <c-\><c-o>zh
inoremap <S-4-ScrollWheelDown> <c-\><c-o>zh

" completion popup keys
inoremap <silent><expr> <Esc>   pumvisible() ? "\<C-e><Esc>" : "\<Esc>"
inoremap <silent><expr> <C-c>   pumvisible() ? "\<C-e><C-c>" : "\<C-c>"
inoremap <silent><expr> <BS>    pumvisible() ? "\<C-e><BS>"  : "\<BS>"
inoremap <silent><expr> <CR>    pumvisible() ? "\<C-e><CR>" : "\<CR>"
inoremap <silent><expr> <TAB>   pumvisible() ? "\<C-y>" : "\<Tab>"
inoremap <silent><expr> <S-TAB> pumvisible() ? "\<C-p>" : "\<C-h>"

" move up/down on display lines instead of logical lines with arrows in insert mode
inoremap <Up> <c-\><c-o>gk
inoremap <Down> <c-\><c-o>gj

" ## coc setup ##

" Use tab for trigger completion with characters ahead and navigate.
" NOTE: Use command ':verbose imap <tab>' to make sure tab is not mapped by
" other plugin before putting this into your config.
" NOTE: not recognizing check_backspace for some reason but we can just use
" ctrl+space to manually show pum since it appears while typing anyways.
" function! s:check_backspace() abort
"   let col = col('.') - 1
"   return !col || getline('.')[col - 1]  =~# '\s'
" endfunction
" inoremap <silent><expr> <TAB>
"       \ pumvisible() ? "\<C-n>" :
"       \ <SID>check_backspace() ? "\<TAB>" :
"       \ coc#refresh()


" Use <c-space> to trigger completion.
inoremap <silent><expr> <c-space> coc#refresh()
