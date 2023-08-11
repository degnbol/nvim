" extra capitalization
Abolish TH{ere,en,e,is} Th{}
" lazy apostrophe
Abolish {can,is,are,do,does,has,have,was,were,would,should,could}nt {}n't
Abolish {let,there,here,who}s {}'s
iabbrev THeres There's
iabbrev im I'm
iabbrev Im I'm
iabbrev IM I'm
Abolish theyre they're
" flipped letters
Abolish ahve have
Abolish sohw show
Abolish blaance balance
" spelling mistakes
Abolish {despa,sepe}rat{e,es,ed,ing,ely,ion,ions,or} {despe,sepa}rat{}
Abolish eucledian Euclidean
" short forms
Abolish algo{,s} algorithm{}
Abolish tho though
Abolish altho although
Abolish thru through
Abolish probs probably
Abolish ppl people
" This one might be annoying for writing dialogue?
Abolish dialog dialogue
Abolish kinda kind of
Abolish combo{,s} combination{}
Abolish hi{,e}ra{,r}ch{y,ical} hi{e}{r}arch{}
Abolish noone no one
" i..e isn't valid as a keyword, so we have a snippet in luasnippets/all.lua 
" for that kind of typo 
Abolish ie i.e.

" language specific, see lua/keymap
function s:ToggleDanskAbbrev() abort
    if &iminsert
        " Dansk
        iabbrev feks f.eks.
        iunabbrev eg
        iunabbrev Eg
        iunabbrev ti
    else
        " English
        " danish unabbrevs has to be silent! since they might not have been 
        " set yet so will give an error message
        silent! iunabbrev feks
        iabbrev eg e.g.
        iabbrev Eg e.g.
        iabbrev ti it
    endif
endfunction
call s:ToggleDanskAbbrev()
augroup ToggleDanskAbbrev
  autocmd User ToggleDansk :call s:ToggleDanskAbbrev()
augroup END

" project specific
iabbrev hic Hi-C
iabbrev HiC Hi-C
