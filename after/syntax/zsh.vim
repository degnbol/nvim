#!/usr/bin/env vimscript
runtime after/syntax/sh.lua

" add sudo to precommands and remove -
" https://github.com/vim/vim/blob/master/runtime/syntax/zsh.vim
" syn clear zshPrecommand
" syn keyword zshPrecommand noglob nocorrect exec command builtin time sudo nextgroup=zshAfterPrecommand
" things after sudo or time are also @function.call
" syn match zshAfterPrecommand ' [a-zA-Z0-9]\+' contained
" zshPrecommand is linked to Special in zsh.vim in vim core but it always 
" colored according to @function.call in the examples where I have found it so 
" far.
" hi link zshPrecommand @function.call
" hi def link zshAfterPrecommand @function.call

syn match @parameter ':[rth]' containedin=zshSubst,zshString

hi link zshString String

" slightly more subtle color for things in `` that aren't recognized as 
" function calls, flag arguments, etc.
hi link zshOldSubst String

" in vim core there are keywords such as "clone" that will be recognized in 
" all contexts and highlighted as a keyword, which we color differently from 
" other function calls and italize. Color them according to functions instead.
hi link zshCommands @function.call

" let variable look the same when defined and expanded, the $ is enough to 
" make it clear
hi link zshDeref @variable
hi link zshSubst @variable

" works for bash but not for zsh for some reason
syn match Operator '\\$'

