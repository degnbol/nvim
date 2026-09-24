" Appended "." to syntax specific iskeyword so that e.g. `git config alias.root ...`
" isn't understood to contain the shell builtin "alias"
syn iskeyword @,48-57,_,192-255,#,-,.
