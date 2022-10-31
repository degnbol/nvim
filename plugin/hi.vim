" Italics are not default since terms aren't assumed to support it,
" however bold is sometimes set for highlight groups for some languages.
" italic for builtin reserved words makes sense.
" highlight groups starting with @ are from treesitter and others from regex 
" syntax groups. Which one is edited matters, since editing one with color 
" will combine the color with the italic, while editing the other may replace 
" the color with italic. Before making changes look at test files in 
" testfiles/
hi @include gui=italic
hi Keyword gui=italic
" effectively clears the guifg on @keyword.operator which let's the regex 
" guifg shine through, fixing "in"
hi @keyword.operator gui=italic
hi Repeat gui=italic
hi Conditional gui=italic
hi Exception gui=italic

" using Operator and not treesitter @operator since the latter replaces rather 
" than combines, so we would lose the coloring. With Operator it is both 
" colored (purple) and bold.
hi Operator gui=bold
