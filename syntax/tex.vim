" Fix missing number recognition in e.g. width=.5\textwidth, \vspace*{-2cm}, etc.
" \ze ends matched region. Then we check for certain chars after match, so we don't match e.g. 214M written in figure caption text.
syntax match Number /[^0-9a-zA-Z_ ]\zs-\?[0-9]*\.\?[0-9]\+\(cm\|mm\|pt\|in\)\?\ze[\\}\]>]/ containedin=texArg,texConcealedArg,texFileOpt,texBeamerOpt,texGroup contained

" Fix missing recognition of beamer arg in e.g. \uncover<5->
syntax region texBeamerOpt matchgroup=texBeamerDelim start='<' end='>' containedin=texGroup

" hack for fixing \subref from subcaption package and custom cmds such as \see
syntax match texRefArg /sec:[A-Za-z0-9_]\+/
syntax match texRefArg /fig:[A-Za-z0-9_]\+/
syntax match texRefArg /tab:[A-Za-z0-9_]\+/
syntax match texRefArg /eq:[A-Za-z0-9_]\+/

hi def link texItemLabelDelim Delimiter
