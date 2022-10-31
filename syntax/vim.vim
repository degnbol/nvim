#!/usr/bin/env vimscript

" to remind myself that vimscript in config files doesn't support the vim9 
" syntax with #. @preproc since this is what the shebang is colored with in 
" e.g. zsh.
hi link vim9Comment @preproc
hi link vim9LineComment @preproc

