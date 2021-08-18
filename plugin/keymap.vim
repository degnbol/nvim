
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
" back space and delete buttons delete should delete also in normal mode
nnoremap <BS> X
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
nmap <ScrollWheelUp>   <c-y>
nmap <2-ScrollWheelUp> <c-y>
nmap <3-ScrollWheelUp> <c-y>
nmap <4-ScrollWheelUp> <c-y>
nmap <ScrollWheelDown>   <c-e>
nmap <2-ScrollWheelDown> <c-e>
nmap <3-ScrollWheelDown> <c-e>
nmap <4-ScrollWheelDown> <c-e>

" shift scroll scrolls horizontally
nnoremap <S-ScrollWheelUp>   zl
nnoremap <S-2-ScrollWheelUp> zl
nnoremap <S-3-ScrollWheelUp> zl
nnoremap <S-4-ScrollWheelUp> zl
nnoremap <S-ScrollWheelDown>   zh
nnoremap <S-2-ScrollWheelDown> zh
nnoremap <S-3-ScrollWheelDown> zh
nnoremap <S-4-ScrollWheelDown> zh

" completion popup keys
ino <silent><expr> <Esc>   pumvisible() ? "\<C-e><Esc>" : "\<Esc>"
ino <silent><expr> <C-c>   pumvisible() ? "\<C-e><C-c>" : "\<C-c>"
ino <silent><expr> <BS>    pumvisible() ? "\<C-e><BS>"  : "\<BS>"
ino <silent><expr> <CR>    pumvisible() ? "\<C-e><CR>" : "\<CR>"
ino <silent><expr> <Tab>   pumvisible() ? "\<C-y>" : "\<Tab>"



