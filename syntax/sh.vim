syn keyword Keyword in
syn keyword @function.builtin mkdir sed tr gzip gunzip rm cd cat mv

syn match Delimiter /\$/
syn match Wildcard /*/
" Hi defined in lua/highlights.lua

" flags with either - or + prefix
syn match @flag / \zs-[A-Za-z_-]\+/
syn match @flag / \zs+[A-Za-z_+]\+/
hi def link @flag Function

" paths.
syn match @path '[A-Za-z0-9.*_-]*[/.][/A-Za-z0-9.*_-]*' contains=FunctionPath,Wildcard
hi def link @path @text.underline

" First word on a line.
syn match Function /^\s*\zs[A-Za-z_-]\+/
" When executing using a path, highlight the executable/script and not the full path as Function.
" I wanted to have a single syn match call, but the underline doesn't seem to 
" combine with the Function color, hence new hi group FunctionPath.
" Regex: \ze[ \n] means there should be a space after the word or end-of-line, 
" and \ze to not highlight those.
" We need to have the path as first word on a line to assume it's a called 
" function so we check for ^\s* before match, however we need to use
" \@<= to look backwards at start of checking, since we are conditioned on 
" @path which doesn't match region from start of line.
" AND to cover cases such as $0:h/path/to/file we also have to look backwards 
" allowing for $0:h, hence \S* added at the end of the atom.
" AND to support e.g. $(git root)/path/to/file which has a space I added the 
" \($(.*)\)\?.
" There's also backticks e.g. `git root` but it's legacy and we prefer $(...) 
" so maybe let's not complicate this regex further.
syn match FunctionPath '\(^\s*\S*\($(.*)\)\?\)\@<=[/A-Za-z0-9._-]*/\zs[A-Za-z._-]\+\ze[ \n]' contained containedin=@path


