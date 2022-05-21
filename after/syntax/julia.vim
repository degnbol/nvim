" conceal multiplication operator * as cdot
" clear juliaOperator means we don't clash so .* isn't concealed.
" It also clears the coloring which is a color I can't tell apart from the 
" Normal white.
syntax clear juliaOperator
syntax match juliaOperator '*' conceal cchar=⋅
" This removes the highlight 
hi clear Conceal
" We get back the color with a link for Conceal
hi link Conceal Operator

" a bit of a hack but it has to be put in this folder to be invoked late 
" enough apparently. See ftplugin/julia.lua for context.
" set indentexpr=
