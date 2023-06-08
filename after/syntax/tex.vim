" I prefer to conceal this
syn match superscript /\\textsuperscript/ conceal cchar=^
syn match subscript /\\textsubscript/ conceal cchar=_
" ... which removes their color since we are giving them a different syntax 
" group. We can color them as function instead the default Statement, since it 
" matches the conceal color.
hi link subscript Function
hi link superscript Function
" underline was also un-highlighted for some reason
hi link texTypeStyle Function

