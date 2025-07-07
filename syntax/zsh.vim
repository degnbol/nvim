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

" in vim core there are keywords such as "clone" that will be recognized in 
" all contexts and highlighted as a keyword, which we color differently from 
" other function calls and italize. Color them according to functions instead.
hi def link zshCommands @function.builtin

" let variable look the same when defined and expanded, the $ is enough to 
" make it clear
hi def link zshVariable @parameter
hi def link zshDeref @parameter
hi def link zshSubst zshDeref
hi def link zshSubstQuoted Special
hi def link zshSubstDelim Delimiter
hi def link zshOperator Operator
hi def link zshParentheses Delimiter

" works for bash but not for zsh for some reason
syn match Operator /\$/ contained containedin=zshDeref,zshShortDeref,zshSubstQuoted,zshSubstDelim,zshSubst

" They are never actually numbers right? Just strings that may be interpreted 
" as numbers by other programs.
hi def link zshNumber None

syn match @path.zshShortDeref /$0/ containedin=zshShortDeref,ZshPathOp
" hi defined in lua/highlights.lua

syn match ZshPathOp /:/ containedin=zshString,zshSubstQuoted contains=@parameter nextgroup=zshPathOpArg
hi def link ZshPathOp Operator
syn match zshPathOpArg '[rth]' contained
hi def link ZshPathOpArg Function

" First word within $(...) is probably a cmd name. Only first word.
" This was hard to figure out. With contained and containedin we start looking 
" for matches at the first position within $(...), i.e. checking at ? in 
" $(?...). Using ^ doesn't work since the begining inside the brackets is 
" never start-of-line. \zs doesn't work since it only works by starting 
" highlighting late within a match. \@<= works since it looks backwards after 
" we start matching at a position. It requires a match, the opposite is 
" accomplished with \@<!
" The required text is then $( made into a single "atom": \( ... \)
syn match Function /\(\$(\)\@<=[A-Za-z_-]\+/ contained containedin=zshSubstQuoted

" Remove colouring everything within `...` with a single simple colour.
hi def link zshOldSubst None
" First word within `...` is probably a cmd name like above.
syn match Function /`\@<=[A-Za-z_-]\+/ contained containedin=zshOldSubst

syn region zshArraySubscript matchgroup=Delimiter start=/\[/ end=/\]/
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
