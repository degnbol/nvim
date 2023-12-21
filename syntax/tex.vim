" Fix missing number recognition in e.g. width=.5\textwidth
syntax match Number /[0-9.-]\+/ containedin=texArg,texConcealedArg,texFileOpt,texBeamerOpt contained

" Fix missing recognition of beamer arg in e.g. \uncover<5->
syntax region texBeamerOpt matchgroup=texBeamerDelim start='<' end='>' containedin=texGroup

