" Fix missing number recognition in e.g. width=.5\textwidth and beamer's \only<1-3> etc
" Pattern starts with \< which means beginning of word. So it will match '.6', ^6, ' 6'
" but not 'B6'
syntax match Number /\<-\?[0-9]*\.\?[0-9]\+/ containedin=texArg,texConcealedArg,texFileOpt,texBeamerOpt contained

" Fix missing recognition of beamer arg in e.g. \uncover<5->
syntax region texBeamerOpt matchgroup=texBeamerDelim start='<' end='>' containedin=texGroup

