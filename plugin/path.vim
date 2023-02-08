" append $PATH to vim option &path (using , instead of : delim)
exec "setlocal path+=" . substitute($PATH, ":", ",", "g")
" prepend git ROOT (doesn't cause issues when not in git repo)
let ROOT = finddir('.git/..', expand('%:p:h').';')
exec "setlocal path^=" . ROOT
" append $ROOT/src
exec "setlocal path+=" . ROOT . '/src'
" append $ROOT/src/*
exec "setlocal path+=" . ROOT . '/src/*'

" edit-in-kitty on remotes doesn't copy the env variables that are normally 
" present locally so we have to set the necessary ones. These are for julia 
" and python lsp to find the right binaries. These were found by copying all 
" of `env` and checking which fixed the problem by adding them to 
" `edit-in-kitty --env ...` which can be given multiple times like the kitty @ 
" launch command.
" Hardcoding miniforge base env path. It would be cooler to use activated env 
" somehow.
" ~/.local/bin is location of jedi-language-server
let $PATH = $HOME . '/.local/bin/:/opt/homebrew/Caskroom/miniforge/base/bin:/opt/homebrew/bin:' . $PATH
let $CONDA_PREFIX = $HOME . '/miniconda3'

