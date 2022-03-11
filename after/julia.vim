" when treesitter runs for julia the indent is wrong after inline for loop,
" but setting the indentexpr to the one from treesitter it works, but then it
" breaks indents after end of function.
" The problem with indent happening when writing a bracket is actually due to
" indentkeys which are a set of keys that when written will trigger autoindent
" of a line. Setting the indent of a line is ALWAYS done by calling
" indentexpr which is GetJuliaIndent() by default. 
" Setting the treesitter config enable flag sets indentexpr=nvim_treesitter#indent()
" For now I compromise by enabling treesitter for julia but disabling
" highlight.
" OK now it somehow doesn't highlight ever? So this is not necessary now for
" some reason.
" TSBufDisable highlight

