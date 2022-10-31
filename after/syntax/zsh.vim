#!/usr/bin/env vimscript
hi @conditional gui=italic
syn keyword Repeat in
hi @parameter guifg=NONE

syn match Arg / [-+]\{1,2}[A-Za-z-]\+/
hi link Arg Special

" waste of limited set of colors to destinguish between different function 
" calls based on how core they are
hi link @function.builtin @function.call

syn match Path '\~\?[a-zA-Z0-9_/:$@.*-]*[/.@*][a-zA-Z0-9_/:$@.*-]*'
" same as above but in double quotes
syn match Path '"[a-zA-Z0-9_/:$@.*-]*[/.@*][a-zA-Z0-9_/:$@.*-]*"'
hi Path gui=underline

" add sudo to precommands and remove -
" https://github.com/vim/vim/blob/master/runtime/syntax/zsh.vim
syn clear zshPrecommand
syn keyword zshPrecommand noglob nocorrect exec command builtin time sudo nextgroup=zshAfterPrecommand
" things after sudo or time are also @function.call
syn match zshAfterPrecommand ' [a-zA-Z0-9]\+' contained
" zshPrecommand is linked to Special in zsh.vim in vim core but it always 
" colored according to @function.call in the examples where I have found it so 
" far.
" hi link zshPrecommand @function.call
hi def link zshAfterPrecommand @function.call

" source with '.' is captured as a path, so we remove underline:
syn match @function.builtin '^. '
" it is still colored according to @function.call since treesitter takes 
" precidence over regex, however it would be too messy to clear @function.call 
" so we accept this.

" treesitter has @comment but I had to make a regex version that covers the 
" whole line since e.g. a single ' inside a comment would break the highlight 
" for the rest of the file.
syn match Comment '\s*#.*'

syn match @function.call ':[rth]' containedin=zshSubst,zshString
" clear treesitter string coloring, allowing coloring of variables inside 
" strings etc.
hi @string guifg=NONE
hi link zshString String

" color treesitter variable which replaces inconsistent coloring of zshDeref 
" (de-reference, i.e. expanding $). This works better for mixing variables 
" into paths, etc. see testfiles/syn.zsh
hi link @variable Identifier

" slightly more subtle color for things in `` that aren't recognized as 
" function calls, flag arguments, etc.
hi link zshOldSubst String

" in vim core there are keywords such as "clone" that will be recognized in 
" all contexts and highlighted as a keyword, which we color differently from 
" other function calls and italize. Color them according to functions instead.
hi link zshCommands @function.call
