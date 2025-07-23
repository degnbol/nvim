" bash shows error for any use of keywords "elif", "else", and "then".
" We use "then" with miller calls, e.g. `mlr cat FILE then cat`
" shCondError
" https://github.com/vim/vim/blob/master/runtime/syntax/sh.vim
" Remove "then":
syn clear shCondError | syn keyword shCondError elif else

