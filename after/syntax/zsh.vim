#!/usr/bin/env vimscript
runtime after/syntax/sh.lua

hi! link zshOperator Operator

" zshString is both "" and '', we need to distinguish to avoid highlighting 
" e.g. : in literal. See syntax/zsh.vim
syn region zshLiteral start="'" end="'"
hi def link zshLiteral zshString
