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

