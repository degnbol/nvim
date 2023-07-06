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
" Maybe too dangerous sometimes with backspace, e.g. when leaving commandline
" nnoremap <BS> "_X
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

" For remote with problematic clipboard we replace the cutlass x that works 
" locally with a mapping where we call y (which is mapped to smartyank) then 
" d.
nmap xx yydd
xmap x ygvd
function CutOperator(type, ...)
    if a:type == 'line'
        normal `[V`]ygvd
    else
        normal `[v`]ygvd
    endif
endfunction
" Operator is a convenience function in plugin/operator.vim
nnoremap <expr> x Operator('CutOperator')

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
" uncomment to make shift have no effect on scroll
" map <S-ScrollWheelUp>   <c-y>
map <S-2-ScrollWheelUp> <c-y>
map <S-3-ScrollWheelUp> <c-y>
map <S-4-ScrollWheelUp> <c-y>
" map <S-ScrollWheelDown>   <c-e>
map <S-2-ScrollWheelDown> <c-e>
map <S-3-ScrollWheelDown> <c-e>
map <S-4-ScrollWheelDown> <c-e>
" insert mode requires going to normal mode temporarily
imap <ScrollWheelUp>   <c-\><c-o><c-y>
imap <2-ScrollWheelUp> <c-\><c-o><c-y>
imap <3-ScrollWheelUp> <c-\><c-o><c-y>
imap <4-ScrollWheelUp> <c-\><c-o><c-y>
imap <ScrollWheelDown>   <c-\><c-o><c-e>
imap <2-ScrollWheelDown> <c-\><c-o><c-e>
imap <3-ScrollWheelDown> <c-\><c-o><c-e>
imap <4-ScrollWheelDown> <c-\><c-o><c-e>
" make shift have no effect on scroll
imap <S-ScrollWheelUp>   <c-\><c-o><c-y>
imap <S-2-ScrollWheelUp> <c-\><c-o><c-y>
imap <S-3-ScrollWheelUp> <c-\><c-o><c-y>
imap <S-4-ScrollWheelUp> <c-\><c-o><c-y>
imap <S-ScrollWheelDown>   <c-\><c-o><c-e>
imap <S-2-ScrollWheelDown> <c-\><c-o><c-e>
imap <S-3-ScrollWheelDown> <c-\><c-o><c-e>
imap <S-4-ScrollWheelDown> <c-\><c-o><c-e>

" horizontal scroll with sideways trackpad 1 column instead of default 6
map <ScrollWheelLeft> zh
map <2-ScrollWheelLeft> <NOP>
map <3-ScrollWheelLeft> <NOP>
map <4-ScrollWheelLeft> <NOP>
map <ScrollWheelRight>   zl
map <2-ScrollWheelRight> <NOP>
map <3-ScrollWheelRight> <NOP>
map <4-ScrollWheelRight> <NOP>
" insert mode requires going to normal mode temporarily
imap <ScrollWheelLeft>   <c-\><c-o>zh
imap <2-ScrollWheelLeft> <NOP>
imap <3-ScrollWheelLeft> <NOP>
imap <4-ScrollWheelLeft> <NOP>
imap <ScrollWheelRight>   <c-\><c-o>zl
imap <2-ScrollWheelRight> <NOP>
imap <3-ScrollWheelRight> <NOP>
imap <4-ScrollWheelRight> <NOP>
" shift should have no effect on side scroll
map <S-ScrollWheelLeft> zh
map <S-2-ScrollWheelLeft> <NOP>
map <S-3-ScrollWheelLeft> <NOP>
map <S-4-ScrollWheelLeft> <NOP>
map <S-ScrollWheelRight>   zl
map <S-2-ScrollWheelRight> <NOP>
map <S-3-ScrollWheelRight> <NOP>
map <S-4-ScrollWheelRight> <NOP>
" insert mode requires going to normal mode temporarily
imap <S-ScrollWheelLeft>   <c-\><c-o>zh
imap <S-2-ScrollWheelLeft> <NOP>
imap <S-3-ScrollWheelLeft> <NOP>
imap <S-4-ScrollWheelLeft> <NOP>
imap <S-ScrollWheelRight>   <c-\><c-o>zl
imap <S-2-ScrollWheelRight> <NOP>
imap <S-3-ScrollWheelRight> <NOP>
imap <S-4-ScrollWheelRight> <NOP>

" typical modifier to arrow key behavior. 
" TODO decide what to do with shift: select things like outside vim, move 
" further than alt, etc.
nnoremap <A-Up> {
nnoremap <A-Down> }
nnoremap <A-Left> b
nnoremap <A-Right> w

xnoremap <A-Up> {
xnoremap <A-Down> }
xnoremap <A-Left> b
xnoremap <A-Right> w

inoremap <A-Up> <c-\><c-o>{
inoremap <A-Down> <c-\><c-o>}
inoremap <A-Left> <c-\><c-o>b
inoremap <A-Right> <c-\><c-o>w
" change to default behaviour: don't exit to normal mode.
inoremap <A-h> <c-\><c-o>h
inoremap <expr> <A-j> v:count == 0 ? "<c-\><c-o>gj" : "j"
inoremap <expr> <A-k> v:count == 0 ? "<c-\><c-o>gk" : "k"
inoremap <A-l> <c-\><c-o>l
inoremap <A-b> <c-\><c-o>b
inoremap <A-w> <c-\><c-o>w
inoremap <A-e> <c-\><c-o>e
inoremap <A-B> <c-\><c-o>B
inoremap <A-W> <c-\><c-o>W
inoremap <A-E> <c-\><c-o>E

" like J(oin) but selects the current container, e.g. brackets cursor is in, 
" and joins it all. Uses https://github.com/andymass/vim-matchup
nmap gJ va%J

" completion popup keys. <C-e> == reject and close pum. <C-y> == accept.
" seems unnecessary with cmp
" inoremap <silent><expr> <Esc>   pumvisible() ? "\<C-y><Esc>" : "\<Esc>"
" inoremap <silent><expr> <C-c>   pumvisible() ? "\<C-e><C-c>" : "\<C-c>"
" inoremap <silent><expr> <BS>    pumvisible() ? "\<C-y><BS>"  : "\<BS>"
" Enter makes newline unless pum is visible AND an item has been selected
" https://vi.stackexchange.com/questions/15092/auto-complete-popup-menu-make-enter-trigger-newline-if-no-item-was-selected
" Also fix unindenting blankline by writing something (a space) and deleting
" it again.
" NOTE: very important, the <c-\><c-o> is needed. I removed it at some point 
" because I thought it was redudant but it makes sure that when typing 
" newline, followed by e.g. else (from indentkeys) that indentexpr is called. 
" If <c-\><c-o> is removed, then indent will no longer be called 
" automatically. <c-\><c-o> goes to normal mode temporarily, <ESC> is to do 
" nothing there, but somehow that means that something has been done so the 
" whitespace line is not reduced to empty line.
" inoremap <silent><expr> <CR>    pumvisible() ? (complete_info().selected == -1 ? "<C-y><CR><c-\><c-o><ESC>" : "<C-y>") : "<CR><c-\><c-o><ESC>"
" above version seems unnecessary with cmp
inoremap <silent><expr> <CR>    "<CR><c-\><c-o><ESC>"
" inoremap <silent><expr> <TAB>   pumvisible() ? "\<C-n>" : "\<Tab>"
" inoremap <silent><expr> <S-TAB> pumvisible() ? "\<C-p>" : "\<C-h>"

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

" when changing to and from danglish with ctrl+6,
" also switch cmp-dictionary language (cmp_dict.lua)
inoremap <silent> <C-6> <C-6><C-\><C-o>:lua CmpDictUpdate()<CR>

" When language is set to danish we can get ;:'" with alt the same way as we 
" do in danglish keyboard input, however that is with the left alt which is 
" deeper in the OS and the mapping below is for when esc+key (^[) is detected which 
" is what is sent to the terminal, e.g. with kitty's setting 'macos_option_as_alt right'.
" Note that they are also available with the ]} and \| keys.
inoremap <A-;> ;
inoremap <A-S-;> :
inoremap <A-'> '
inoremap <A-S-'> "
" And in case danglish keyboard is active:
inoremap <A-æ> ;
inoremap <A-S-æ> :
inoremap <A-ø> '
inoremap <A-S-ø> "

" completion mapping for function keys available on mechanical keyboard
imap <F13> <C-space>
imap <F14> <C-p>
imap <F15> <C-n>

" use the following two commands to enable spelling
" setlocal spell
" set spelllang=en_us
" ctrl+s (in insert mode) to fix last misspelled word
" credit: https://castel.dev/post/lecture-notes-1/
" with git: https://github.com/gillescastel/latex-snippets
inoremap <C-s> <c-g>u<Esc>[s1z=`]a<c-g>u

