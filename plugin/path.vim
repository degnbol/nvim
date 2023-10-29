" append $PATH to vim option &path (using , instead of : delim)
exec "setlocal path+=" . substitute($PATH, ":", ",", "g")
" prepend git ROOT (doesn't cause issues when not in git repo)
" Note the '$'. It makes ROOT and env var which means it is available for e.g. 
" path completion so that $ROOT/ will have completion, e.g. in julia or zsh.
" Note that for a nested git repo, we will focus on the immediate root, but 
" also add a reference to the top level root.
let $ROOTTOP = finddir('.git/..', expand('%:p:h').';')
let $ROOT = trim(system('git root'))
exec "setlocal path^=" . $ROOT
exec "setlocal path+=" . $ROOTTOP
exec "setlocal path+=" . $ROOT . '/src'
exec "setlocal path+=" . $ROOT . '/src/*'

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

