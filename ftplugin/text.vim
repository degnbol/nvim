" insert real tab by default
setlocal noexpandtab
" there was some annoying issue with indents inserted in regular text
setlocal nosmartindent
" for text we want to see tabs and spaces by default
setlocal list

" I copied the default ftplugin code here and rm a line unsetting 
" commentstring
if exists('b:did_ftplugin')
  finish
endif
let b:did_ftplugin = 1

let b:undo_ftplugin = 'setlocal comments< commentstring<'

" Pseudo comment leaders to indent bulleted lists with '-' and '*'.  And allow
" for Mail quoted text with '>'.
setlocal comments=fb:-,fb:*,n:>
