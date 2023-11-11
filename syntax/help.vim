" linked to @text by default which is not highlighted
hi! link @text.literal.vimdoc @string
hi! link @text.literal.block.vimdoc @string
" linked to @string by default, which means we can't distinguish it from the literals above.
hi! link @string.special.vimdoc @special
" at begining of each section
hi! link @label.vimdoc @text.title

