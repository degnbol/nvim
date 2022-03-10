" insert real tab by default
setlocal noexpandtab
" there was some annoying issue with indents inserted in regular text
set nosmartindent
" t=use textwidth for formatting. NO a=auto format since it is too strong an
" assumtion for general text files, e.g. requirements.txt should respect the
" explicit newline. 
set formatoptions+=t
