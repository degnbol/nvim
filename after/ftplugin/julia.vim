" julia-vim sets it to #=%s=# by default https://github.com/JuliaEditorSupport/julia-vim/blob/master/ftplugin/julia.vim
" Has to be in after/ since we want to override the call loading the file 
" linked above.
setlocal commentstring=#%s

