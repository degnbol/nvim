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
" instead of <c-\><c-o>l since then we can't go to the last char of the line
inoremap <A-l> <right>
inoremap <A-b> <c-\><c-o>b
inoremap <A-w> <c-\><c-o>w
inoremap <A-e> <c-\><c-o>e
inoremap <A-B> <c-\><c-o>B
inoremap <A-W> <c-\><c-o>W
inoremap <A-E> <c-\><c-o>E
" delete back word (insert/cmdmode)
noremap! <A-backspace> <c-w>
inoremap <A-Left> <c-\><c-o>b
" a little movement in cmdmode
cnoremap <C-a> <Home>
cnoremap <expr> <A-left> husk#left()
cnoremap <expr> <A-right> husk#right()
cnoremap <expr> <A-b> husk#left()
cnoremap <expr> <A-w> husk#right()

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

" typo help
command Q q
command X x
command WQ wq
command Wq wq

