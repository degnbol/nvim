" insert real tab by default
setlocal noexpandtab
" there was some annoying issue with indents inserted in regular text
set nosmartindent
" let's try autoformat, and if it messes with code, environments or indented code, 
" then we should find a plugin or something to limit it to simple text.
" If there are problems then look into
" https://www.vim.org/scripts/script.php?script_id=2307
" https://www.vim.org/scripts/script.php?script_id=2187
set fo+=a
" lists doesn't start with number on first column but uses \item
set fo-=n
" we will often be writing prose that wraps at textwidth columns but not code.
" Let's try with sidescrolloff so the window doesn't scroll when we get close 
" to the edge.
set sidescrolloff=0

