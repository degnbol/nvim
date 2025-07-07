syn keyword Keyword in
syn keyword @function.builtin mkdir sed tr gzip gunzip rm cd cat mv

syn match Operator /\$/
syn match Operator /*/

" flags with either - or + prefix
syn match @flag / \zs-[A-Za-z_-]\+/
syn match @flag / \zs+[A-Za-z_+]\+/
hi def link @flag Function

" paths
syn match @path '[A-Za-z0-9.*_-]*[/.][/A-Za-z0-9.*_-]*' containedin=@function.path contains=Operator
hi def link @path @text.underline

" First word on a line.
syn match Function /^\s*\zs[A-Za-z_-]\+/
" I wanted to have a single syn match call, but the underline doesn't seem to 
" combine with the Function color.
syn match @function.path '^\s*\zs[A-Za-z._-]*/[A-Za-z._-]\+'


