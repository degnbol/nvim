#!/usr/bin/env vimscript
runtime after/syntax/sh.lua
hi link shLoop Repeat

" it is still colored according to @function.call since treesitter takes 
" precidence over regex, however it would be too messy to clear @function.call 
" so we accept this.

" treesitter has @comment but I had to make a regex version that covers the 
" whole line since e.g. a single ' inside a comment would break the highlight 
" for the rest of the file.
" syn match Comment '\s*#.*'

" stop bash from saying then is wrong when used in miller by modifying the 
" shCondError https://github.com/vim/vim/blob/master/runtime/syntax/sh.vim
syn clear shCondError
syn keyword shCondError elif else

