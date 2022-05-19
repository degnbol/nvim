" conceal multiplication operator * as cdot
hi clear Conceal
" otherwise we clash and .* isn't concealed. This removes the highlight 
" coloring for julia operators but I cannot tell it apart from the default 
" white textcolor.
syntax clear juliaOperator
syntax match juliaOperator '*' conceal cchar=⋅

" a bit of a hack but it has to be put in this folder to be invoked late 
" enough apparently. See ftplugin/julia.lua for context.
" set indentexpr=
