runtime syntax/sh.vim

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
hi def link zshCommands @function.builtin

" let variable look the same when defined and expanded, the $ is enough to 
" make it clear
hi def link zshVariable @parameter
hi def link zshDeref @parameter
hi def link zshSubstQuoted Special
hi def link zshSubstDelim Special
hi def link zshOperator Operator

" works for bash but not for zsh for some reason
syn match Operator /\$/ contained containedin=zshDeref

" They are never actually numbers right? Just strings that may be interpreted 
" as numbers by other programs.
hi def link zshNumber None

syn match @path.zshShortDeref /$0/ containedin=zshShortDeref,ZshPathOp
" hi defined in lua/highlights.lua

syn match ZshPathOp /:/ containedin=zshString,zshSubstQuoted contains=@parameter nextgroup=zshPathOpArg
hi def link ZshPathOp Operator
syn match zshPathOpArg '[rth]' contained
hi def link ZshPathOpArg Function

