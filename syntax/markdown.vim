" markdownItalic and Bold are linked to html which was linked to nothing.
" some things with injected code resets highlights (maybe a bug), but this 
" defualts call works.
hi! default link htmlItalic Italic
hi! default link htmlBold Bold

" When using just treesitter, a lot of things are not colored by default but 
" are in fact captured by the treesitter parser.
hi @text.title gui=bold
" Different colors for different levels.
" The rest will just be bold white if I ever use them.
hi! default link @text.title1 @define
hi! default link @text.title2 @constructor
hi! default link @text.title3 @function

hi @text.emphasis gui=italic
hi @text.strong gui=bold

