#!/usr/bin/env vimscript
source after/syntax/sh.vim

" add sudo to precommands and remove -
" https://github.com/vim/vim/blob/master/runtime/syntax/zsh.vim
syn clear zshPrecommand
syn keyword zshPrecommand noglob nocorrect exec command builtin time sudo nextgroup=zshAfterPrecommand
" things after sudo or time are also @function.call
syn match zshAfterPrecommand ' [a-zA-Z0-9]\+' contained
" zshPrecommand is linked to Special in zsh.vim in vim core but it always 
" colored according to @function.call in the examples where I have found it so 
" far.
" hi link zshPrecommand @function.call
hi def link zshAfterPrecommand @function.call

syn match @function.call ':[rth]' containedin=zshSubst,zshString

hi link zshString String

" slightly more subtle color for things in `` that aren't recognized as 
" function calls, flag arguments, etc.
hi link zshOldSubst String

" in vim core there are keywords such as "clone" that will be recognized in 
" all contexts and highlighted as a keyword, which we color differently from 
" other function calls and italize. Color them according to functions instead.
hi link zshCommands @function.call
