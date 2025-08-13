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

syn match Delimiter /\$/ contained containedin=shDeref,shDerefSimple
