syntax match Operator /::=/
syntax match Delimiter /[;|()]/
syntax match Statement /^[A-Za-z0-9-]\+/
syntax match Flag /--\?[A-Za-z0-9-]\+/
" String placed after flag, so it takes precidence over flags written in a 
" string. We could also just say it can't be contained and all that.
" skip anything that is escaped, it might be an escaped quotation mark.
syntax region String start=/\v"/ skip=/\v\\./ end=/\v"/
" placed after flags, again to take precidence.
syntax match Identifier /<[A-Za-z0-9_-]\+>/
" TODO: maybe add more but this is enough to get started.
