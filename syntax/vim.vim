" captures comments starting with # which isn't a valid comment in neovim 
" vimscript.
syntax match VimGroupName /@[a-z.]*/

" Could be considered a type or a variable.
" By default linked to Normal maybe through variable.
" Identifier currently allows it to be highlighted.
hi def link VimGroupName Identifier
hi def link VimGroupList Identifier
hi def link VimGroup Identifier
" This is e.g. keyword `ALL`
hi def link VimGroupSpecial @identifier.builtin
" The word(s) defined as keyword in `syn keyword <VimGroupName> <VimSynKeyRegion>...`
hi def link VimSynKeyRegion String
hi def link VimSynContains @parameter
hi def link VimSynMtchGrp @parameter
hi def link VimSynReg @parameter

" Function definition
hi def link vimFunction @function
" Was also missing
hi def link vimSubscriptBracket Delimiter
