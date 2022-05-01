" Fugitive: git on commandline
cnoreabbrev Gs Git status -uno
cnoreabbrev Gp Git push
cnoreabbrev Gl Git pull
cnoreabbrev Ga Git add %
cnoreabbrev Gc Git commit -m


" substitute motion with clipboard content with https://github.com/svermeulen/vim-subversive
nmap s  <plug>(SubversiveSubstitute)
nmap ss <plug>(SubversiveSubstituteLine)
nmap S  <plug>(SubversiveSubstituteToEndOfLine)
" Integration with yoink
xmap s <plug>(SubversiveSubstitute)
xmap p <plug>(SubversiveSubstitute)
xmap P <plug>(SubversiveSubstitute)
" example: <leader>siwip to replace all instances of the current word under the cursor that exist within the paragraph under the cursor. 
" example: <leader>sl_ to replace all instances of the character under the cursor on the current line.
" example: <leader>ssip to replace the word under cursor in the current paragraph. Matches complete words so is different from <leader>siwip
nmap <leader>s <plug>(SubversiveSubstituteRange)
xmap <leader>s <plug>(SubversiveSubstituteRange)
nmap <leader>ss <plug>(SubversiveSubstituteWordRange)

" Yank, cut, delete behavior
" a logical choice that was supposed to be default? https://github.com/neovim/neovim/issues/416
map Y y$
" back space and delete buttons delete should delete also in normal mode.
" Delete to black hole register.
nnoremap <BS> "_X
" Yoink:
nmap <c-n> <plug>(YoinkPostPasteSwapBack)
nmap <c-p> <plug>(YoinkPostPasteSwapForward)
nmap p <plug>(YoinkPaste_p)
nmap P <plug>(YoinkPaste_P)
nmap gp <plug>(YoinkPaste_gp)
nmap gP <plug>(YoinkPaste_gP)
nmap [y <plug>(YoinkRotateBack)
nmap ]y <plug>(YoinkRotateForward)
" toggle if pasting with == formatting or not
nmap <c-=> <plug>(YoinkPostPasteToggleFormat)
nmap y <plug>(YoinkYankPreserveCursorPosition)
xmap y <plug>(YoinkYankPreserveCursorPosition)

" when adding new line below or above, write something (a space) and delete it
" so the indent aren't removed on ESC.
nmap o o<c-\><c-o> <BS>
nmap O O<c-\><c-o> <BS>

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
noremap <S-ScrollWheelUp>   zh
noremap <S-2-ScrollWheelUp> <NOP>
noremap <S-3-ScrollWheelUp> <NOP>
noremap <S-4-ScrollWheelUp> <NOP>
noremap <S-ScrollWheelDown> zl
noremap <S-2-ScrollWheelDown> <NOP>
noremap <S-3-ScrollWheelDown> <NOP>
noremap <S-4-ScrollWheelDown> <NOP>
" insert mode requires going to normal mode temporarily
inoremap <S-ScrollWheelUp>   <c-\><c-o>zh
inoremap <S-2-ScrollWheelUp> <NOP>
inoremap <S-3-ScrollWheelUp> <NOP>
inoremap <S-4-ScrollWheelUp> <NOP>
inoremap <S-ScrollWheelDown>   <c-\><c-o>zl
inoremap <S-2-ScrollWheelDown> <NOP>
inoremap <S-3-ScrollWheelDown> <NOP>
inoremap <S-4-ScrollWheelDown> <NOP>

" completion popup keys. <C-e> == reject and close pum. <C-y> == accept.
inoremap <silent><expr> <Esc>   pumvisible() ? "\<C-y><Esc>" : "\<Esc>"
inoremap <silent><expr> <C-c>   pumvisible() ? "\<C-e><C-c>" : "\<C-c>"
inoremap <silent><expr> <BS>    pumvisible() ? "\<C-y><BS>"  : "\<BS>"
" Enter makes newline unless pum is visible AND an item has been selected
" https://vi.stackexchange.com/questions/15092/auto-complete-popup-menu-make-enter-trigger-newline-if-no-item-was-selected
" Also fix unindenting blankline by writing something (a space) and deleting
" it again.
inoremap <silent><expr> <CR>    pumvisible() ? (complete_info().selected == -1 ? "<C-y><CR> <BS>" : "<C-y>") : "<CR> <BS>"
inoremap <silent><expr> <TAB>   pumvisible() ? "\<C-n>" : "\<Tab>"
inoremap <silent><expr> <S-TAB> pumvisible() ? "\<C-p>" : "\<C-h>"

" move up/down on display lines instead of logical lines with arrows in insert mode
" inoremap <Up> <c-\><c-o>gk
" inoremap <Down> <c-\><c-o>gj

" Danglish support. For when Danglish keyboard is selected, 
" generally you should instead stay in code keyboard and use iminsert=2
" This can also be done with langmap but since these are 
" never used in normal mode then it doesn't hurt.
nnoremap æ ;
nnoremap Æ :
nnoremap ø '
nnoremap Ø "
nnoremap å [
nnoremap Å {
xnoremap æ ;
xnoremap Æ :
xnoremap ø '
xnoremap Ø "
xnoremap å [
xnoremap Å {
" In Danglish I moved : and ; to the }] button
" But this messes with things when I'm not in Danglish.

" ## coc setup ##

" Use <c-space> to trigger completion.
inoremap <silent><expr> <c-space> coc#refresh()

" default = 4000 [ms]
set updatetime=1000


