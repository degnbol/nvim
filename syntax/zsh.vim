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

hi link zshString String

" in vim core there are keywords such as "clone" that will be recognized in 
" all contexts and highlighted as a keyword, which we color differently from 
" other function calls and italize. Color them according to functions instead.
hi def link zshCommands @function.builtin

hi def link zshPrecommand @function.builtin

" let variable look the same when defined and expanded, the $ is enough to 
" make it clear
hi def link zshVariable @parameter
hi def link zshDeref @parameter
hi def link zshSubst zshDeref
hi def link zshSubstQuoted Special
hi def link zshSubstDelim Delimiter
hi def link zshOperator Operator
hi def link zshParentheses Delimiter

syn match Delimiter /\$/ contained containedin=zshDeref,zshShortDeref,zshSubstQuoted,zshSubstDelim,zshSubst,zshArraySubscript,zshString nextgroup=zshPathOp

" They are never actually numbers right? Just strings that may be interpreted 
" as numbers by other programs.
hi def link zshNumber None

syn match zshShortDerefPath /$0/ containedin=zshShortDeref nextgroup=zshPathOp
" hi defined in lua/highlights.lua

" The hs=e-1 after the pattern means highlight starts at end of match minus 1.
" I tried first with \zs but it didn't work for some reason.
syn match zshPathOp /\w*:[rth]/hs=e-1 contained containedin=zshDeref,zshShortDeref,zshShortDerefPath,zshSubstQuoted,zshString contains=zshPathOpArg
hi def link zshPathOp Delimiter
syn match zshPathOpArg '[rth]' contained nextgroup=zshPathOp
hi def link zshPathOpArg @function.builtin

" Remove colouring everything within `...` with a single simple colour.
hi def link zshOldSubst None
" First word within `...` is probably a cmd name like above.
syn match Function /`\@<=[A-Za-z_-]\+/ contained containedin=zshOldSubst

syn region zshArraySubscript matchgroup=Delimiter start=/\[/ end=/\]/
" E.g. in ${@:1:$#-1}
syn match Delimiter /:/ contained containedin=zshSubst
" E.g. numbers in ${@:1:$#-1}, however also number in ${0}
syn match Number /-\?\d\+/ contained containedin=zshArraySubscript,zshSubst
" Script arguments
syn match PreProc /[@#]/ contained containedin=zshSubst
" Unpacking operator e.g. "=" in "{=@}"
syn match Operator /\(${\)\@<==/ contained containedin=zshSubst

" Highlight = in variable assignments.
" It is checked to be preceded by a valid variable name using \@<=.
" \h == head of word, i.e. [A-Za-z_], and \w == word, i.e. the same plus 
" digits.
syn match Operator /\(\h\+\w*\)\@<==/
" [ ... ] used for test expressions is matched as zshArraySubscript.
" Match = in e.g. `[ "$var" = "value" ]`
syn match Operator /=/ contained containedin=zshArraySubscript
syn region zshString matchgroup=zshStringDelimiter start=/"/ skip=/\\"/ end=/"/ contained containedin=zshArraySubscript

syn match Operator /!/
" By default matches as an operator.
syn match Delimiter /;/ containedin=zshOperator

" zshString is both "" and '', we need to distinguish to avoid highlighting 
" e.g. ":" in literal.
syn region zshLiteral start="'" end="'"
hi def link zshLiteral zshString

syn match Delimiter /\[\[/
syn match Delimiter /\]\]/
