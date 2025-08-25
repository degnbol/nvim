syn match Operator '||'
syn match Operator '&&'
" stdin/stdout dash
syn match Special '\<-\>'
" pipe
syn match Operator '|'

hi def link shVarAssign Operator
hi def link shTestOpr Operator
hi def link shQuote String
hi def link shOption @flag
hi def link shCommandSub None
hi def link shDerefSimple @parameter
hi def link bashStatement @function.builtin
hi def link kshStatement @function.builtin
hi def link shVariable Variable
hi def link shCmdSubRegion Delimiter
" Function definition
hi def link shFunction @function

syn match Delimiter /\$/ contained containedin=shDeref,shDerefSimple
