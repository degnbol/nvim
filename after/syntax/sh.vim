#!/usr/bin/env vimscript
" recolor "in" from a Statement color to a repeat color
hi @conditional gui=italic
hi link shLoop Repeat

syn keyword Repeat in

au BufEnter *.zsh,*.sh,*.bash hi @parameter guifg=NONE
au BufLeave *.zsh,*.sh,*.bash hi link @parameter Identifier

" color treesitter variable which replaces inconsistent coloring of zshDeref 
" (de-reference, i.e. expanding $). This works better for mixing variables 
" into paths, etc. see testfiles/syn.zsh
" NOTE: since :hi is a global modifier and @variable is cleared by default and 
" used by other filetypes we must change it each time we change focus
" hi link @variable Identifier
au BufEnter *.zsh,*.sh,*.bash hi link @variable Identifier
au BufLeave *.zsh,*.sh,*.bash hi clear @variable

syn match Arg / [-+]\{1,2}[A-Za-z0-9-]\+/
hi link Arg Special

" waste of limited set of colors to destinguish between different function 
" calls based on how core they are
hi link @function.builtin @function.call

syn match Path '\~\?[a-zA-Z0-9_/:$@.*-]*[/.@*][a-zA-Z0-9_/:$@.*-]*'
" same as above but in double quotes
syn match Path '"[a-zA-Z0-9_/:$@.*-]*[/.@*][a-zA-Z0-9_/:$@.*-]*"'
hi Path gui=underline

" source with '.' is captured as a path, so we remove underline:
syn match @function.builtin '^. '
" it is still colored according to @function.call since treesitter takes 
" precidence over regex, however it would be too messy to clear @function.call 
" so we accept this.

" treesitter has @comment but I had to make a regex version that covers the 
" whole line since e.g. a single ' inside a comment would break the highlight 
" for the rest of the file.
syn match Comment '\s*#.*'

" clear treesitter string coloring, allowing coloring of variables inside 
" strings etc.
hi @string guifg=NONE

" A terminal \ means no newline
syn match Operator '\\$'

" estra goodies for specific tools such as miller
syn match MillerOperator '+'
hi link MillerOperator Operator


