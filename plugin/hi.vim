" Italics are not default since terms aren't assumed to support it,
" however bold is sometimes set for highlight groups for some languages.
" italic for builtin reserved words makes sense.
hi @include gui=italic
hi @keyword gui=italic
hi @keyword.operator gui=italic
" Not setting `hi Keyword gui=italic` since the @ groups are treesitter which are more accurate. 
" E.g. = and other operations are Keyword in python.
hi Repeat gui=italic
hi Conditional gui=italic
hi Exception gui=italic

" using Operator and not treesitter @operator since the latter replaces rather 
" than combines, so we would lose the coloring. With Operator it is both 
" colored (purple) and bold.
hi Operator gui=bold
