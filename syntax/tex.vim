" Fix missing number recognition in e.g. width=.5\textwidth and beamer's \only<1-3> etc
syntax match Number /-\?[0-9]*.\?[0-9]\+/ containedin=texArg,texConcealedArg,texFileOpt,texBeamerOpt contained

" Fix missing recognition of beamer arg in e.g. \uncover<5->
syntax region texBeamerOpt matchgroup=texBeamerDelim start='<' end='>' containedin=texGroup

