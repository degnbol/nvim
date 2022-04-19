" define the syntax identifier texComment with regex for lines starting with %
syntax match texComment /%.*/
" color tex comment as whatever color we have set for comments in general
highlight link texComment Comment
