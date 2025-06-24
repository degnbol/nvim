" Incorrect indent for incomplete expression.
" I don't always write code from top to bottom, so while editing within code 
" it will guess indent incorrectly, and then trigger auto indent with keys in 
" indentkeys, e.g. ] resulting in line back moving while editing it.
setlocal indentkeys=o,O
