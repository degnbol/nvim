#!/usr/bin/env vimscript
syn match Shebang '^#!.*'
hi link Shebang @preproc

hi def link vimIsCommand @attribute.builtin
