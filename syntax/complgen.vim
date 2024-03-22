syntax match Operator /::=/
syntax match Operator /\.\.\./
syntax match Delimiter /[;|()\[\]]/
" Main command name that we are writing completion for.
" Only making it italic leaves it the same colour as subcmds (uncoloured)
syntax match Italic /^[A-Za-z0-9-]\+/
syntax match Flag /--\?[A-Za-z0-9-]\+/
syntax match Flag /++\?[A-Za-z0-9-]\+/
" E.g. +<cmd> (vim flag)
syntax match Flag /--\?\ze[<(\[]/
syntax match Flag /++\?\ze[<(\[]/
" String placed after flag, so it takes precidence over flags written in a 
" string. We could also just say it can't be contained and all that.
" skip anything that is escaped, it might be an escaped quotation mark.
syntax region String start=/\v"/ skip=/\v\\./ end=/\v"/
" placed after flags, again to take precidence.
syntax match Identifier /<[A-Za-z0-9@_-]\+>/
" code verbatim. TODO: injection
syntax region Code matchgroup=Delimiter start=/{{{/ skip=/\v\\./ end=/}}}/
hi def link Code Special
syntax match Comment /^#.*/
