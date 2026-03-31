" neovim treats bash more or less as sh, but has a variable vim.b.is_bash or 
" in vimscript b:is_bash which is 1 if the shell is bash.
" if exists("b:is_bash") && b:is_bash
    runtime syntax/bash.vim
" endif

syn keyword Keyword in
syn keyword @function.builtin mkdir sed tr gzip gunzip rm cd cat mv

" Final \ at end of line for line continuation.
syn match Comment /\\$/ containedin=shCommandSub,shFunctionTwo
syn match Wildcard /*/

" First word on a line or only preceded by an env var assignment (e.g. `LC_ALL=C tr ...`)
" '+' appears in e.g. `g++`
syn match Function /^\s*\zs[A-Za-z_+-]\+/ containedin=shIf,shCommandSub
syn match Function /\(^\h\+\w*=\w\+\)\@<= [A-Za-z_+-]\+/ containedin=shIf,shCommandSub
" ... or after pipe
syn match Function /|\@<= *[A-Za-z_-]\+/ containedin=shIf,shCommandSub
" ... or after `exec`
syn match Function /\(exec\)\@<= *[A-Za-z_-]\+/

" First word within $(...) is probably a cmd name. Only first word.
" This was hard to figure out. With contained and containedin we start looking 
" for matches at the first position within $(...), i.e. checking at ? in 
" $(?...). Using ^ doesn't work since the begining inside the brackets is 
" never start-of-line. \zs doesn't work since it only works by starting 
" highlighting late within a match. \@<= works since it looks backwards after 
" we start matching at a position. It requires a match, the opposite is 
" accomplished with \@<!
" The required text is then $( made into a single "atom": \( ... \)
syn match Function /\(\$(\)\@<=[A-Za-z_-]\+/ contained containedin=zshSubstQuoted,shCommandSub

" Change from default shOperator match.
syn match Delimiter /;/ contained containedin=shIf,shOperator

" {a,b} which expands to multiple alternatives.
syn region ExpandAlts matchgroup=Delimiter start=/{/ end=/}/
syn match ExpandAltsComma /,/ contained containedin=ExpandAlts
hi def link ExpandAltsComma Delimiter

