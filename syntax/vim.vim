" captures comments starting with # which isn't a valid comment in neovim 
" vimscript.
syntax match VimGroupName /@[a-z.]*/

" Function definition
hi def link vimFunction @function
" Was also missing
hi def link vimSubscriptBracket Delimiter
