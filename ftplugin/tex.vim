" insert real tab by default
setlocal noexpandtab
" there was some annoying issue with indents inserted in regular text
set nosmartindent
" lists doesn't start with number on first column but uses \item
set fo-=n
" we will often be writing prose that wraps at textwidth columns but not code.
" Let's try with sidescrolloff so the window doesn't scroll when we get close 
" to the edge.
set sidescrolloff=0

" vimtex has a lot of nice default conceals, e.g. greek in maths, \textbf, etc.
set conceallevel=2

" save on calls to this with event InsertEnter but if we sometimes move to and 
" from math mode without going in and out of normal mode then use CursorMovedI
" as well.
autocmd InsertEnter *.tex if vimtex#syntax#in_mathzone() | set fo-=a | else | set fo+=a | endif

