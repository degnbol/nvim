" conceal multiplication operator * as cdot
hi clear Conceal
syntax match juliaOperator '*' conceal cchar=⋅

" a bit of a hack but it has to be put in this folder to be invoked late 
" enough apparently. See ftplugin/julia.lua for context.
" set indentexpr=
