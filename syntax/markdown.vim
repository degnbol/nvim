" Instead of linking them to htmlItalic and htmlBold
highlight link markdownItalic Italic
highlight link markdownBold Bold

" When using just treesitter, a lot of things are not colored by default but 
" are in fact captured by the treesitter parser.
hi @text.title gui=bold
" Different colors for different levels.
" The rest will just be bold white if I ever use them.
hi link @text.title1 @define
hi link @text.title2 @constructor
hi link @text.title3 @function

hi @text.emphasis gui=italic
hi @text.strong gui=bold

" temp solution for inline code
hi link @code_span String

" TODO have conceallevel=1 hide ` with space so it is invisible but doesn't 
" move text
set conceallevel=0
