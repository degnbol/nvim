" Fix missing number recognition in e.g. width=.5\textwidth and beamer's \only<1-3> etc
" Matching not space nor number [^ 0-9] after the capture (\ze) is to avoid when 
" numbers are in an arg as part of regular text.
syntax match Number /-\?[0-9]*\.[0-9]\+\(cm\|mm\|pt\|in\)\?\ze[^ 0-9]/ containedin=texArg,texConcealedArg,texFileOpt,texBeamerOpt,texGroup contained

" Fix missing recognition of beamer arg in e.g. \uncover<5->
syntax region texBeamerOpt matchgroup=texBeamerDelim start='<' end='>' containedin=texGroup

" hack for fixing \subref from subcaption package and custom cmds such as \see
syntax match texRefArg /sec:[A-Za-z0-9_]\+/
syntax match texRefArg /fig:[A-Za-z0-9_]\+/
syntax match texRefArg /tab:[A-Za-z0-9_]\+/
syntax match texRefArg /eq:[A-Za-z0-9_]\+/

hi def link texItemLabelDelim Delimiter
