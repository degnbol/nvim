" Italics are not default since terms aren't assumed to support it,
" however bold is sometimes set for highlight groups for some languages.
" italic for builtin reserved words makes sense.
hi @include gui=italic
hi @keyword gui=italic
hi @keyword.operator gui=italic
" not italic Keyword since the @ groups are treesitter which is more accurate. 
" E.g. = and other operations are Keyword in python.
hi @repeat gui=italic
hi @conditional gui=italic
hi @exception gui=italic

hi @operator gui=bold
