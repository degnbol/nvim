" I prefer to conceal this
syn match superscript /\\textsuperscript/ conceal cchar=^
syn match subscript /\\textsubscript/ conceal cchar=_
" ... which removes their color since we are giving them a different syntax 
" group.
" We can color them as Function to match conceal color. Default is Statement. 
" We can also go with texCmdStyle which is given to \texttt which links to 
" Type.
hi link subscript texTypeStyle
hi link superscript texTypeStyle
" underline was also un-highlighted for some reason
hi link texTypeStyle texCmdStyle

