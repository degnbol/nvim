" conceal multiplication operator * as cdot
" clear juliaOperator means we don't clash so .* isn't concealed.
" It also clears the coloring which is a color I can't tell apart from the 
" Normal white.
syntax clear juliaOperator
syntax match juliaOperator '*' conceal cchar=⋅
" This removes the highlight 
" hi clear Conceal
" Instead I color it according to operator, copied from :hi operator, 
" since :hi link didn't work
hi Conceal ctermfg=5 guifg=#ae94f9

" a bit of a hack but it has to be put in this folder to be invoked late 
" enough apparently. See ftplugin/julia.lua for context.
" set indentexpr=


" treesitter colors in as @keyword.operator, but I would rather it be considered @repeat to match 'for'
" The correct way to do this would probably be to make a custom capture group 
" in treesitter but I couldn't figure that out.
hi link @keyword.operator @repeat
syn keyword Repeat in
syn match Repeat '∈'

