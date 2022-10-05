" append $PATH to vim option &path (using , instead of : delim)
exec "setlocal path+=" . substitute($PATH, ":", ",", "g")
" prepend git ROOT (doesn't cause issues when not in git repo)
let ROOT = finddir('.git/..', expand('%:p:h').';')
exec "setlocal path^=" . ROOT
" append $ROOT/src
exec "setlocal path+=" . ROOT . '/src'
" append $ROOT/src/*
exec "setlocal path+=" . ROOT . '/src/*'

" For edit-in-kitty on remotes:
let $PATH .= ':/opt/homebrew/bin:~/.local/bin:~/miniconda3/bin'
