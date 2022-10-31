#!/usr/bin/env vimscript
" recolor "in" from a Statement color to a repeat color
hi @conditional gui=italic
hi link shLoop Repeat

syn keyword Repeat in


fun! LocalHi()
    if &ft =~ 'sh\|zsh'
        hi @parameter guifg=NONE
        " color treesitter variable which replaces inconsistent coloring of zshDeref 
        " (de-reference, i.e. expanding $). This works better for mixing variables 
        " into paths, etc. see testfiles/syn.zsh
        " NOTE: since :hi is a global modifier and @variable is cleared by default and 
        " used by other filetypes we must change it each time we change focus
        " hi link @variable Identifier
        hi link @variable Identifier
    elseif &ft =~ 'NvimTree\|TelescopePrompt'
        " Don't do anything when just swithing to an unaffected buffer, e.g. 
        " treeviewer NvimTree or buffer without filetype (empty string)
    elseif &ft == ''
        " exactly nothing didn't work in the regex, it matches anything then
    else
        " reset for other languages
        hi link @parameter Identifier
        hi clear @variable
    endif
endfun
au BufEnter * call LocalHi()


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


