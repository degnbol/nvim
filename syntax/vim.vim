" captures comments starting with # which isn't a valid comment in neovim 
" vimscript.
hi! default link vim9comment Error

syntax match VimGroupName /@[a-z.]*/
