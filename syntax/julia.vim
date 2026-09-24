" Set flag to show syntax is loaded so no other syntax file will be sourced
let b:current_syntax = "julia"
" Define custom regex syntax highlights missing from treesitter highlights
syn match @punctuation.delimiter ';'
syn keyword @variable.builtin stdin stdout stderr
