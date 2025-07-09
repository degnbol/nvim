iabbrev fucntion function
" extra capitalization
Abolish TH{ere,en,e,is} Th{}
" lazy apostrophe
Abolish {ca,is,are,do,does,did,has,have,had,was,were,would,should,could,wo}nt {}n't
Abolish {let,there,here,who}s {}'s
Abolish itøs it's
iabbrev THeres There's
" I by itself is useful in regular text but i can be a variable, and also 
" mentioned in comments, so i should only get abb to I in prose, if possible. 
" iabbrev im I'm " TODO: OMG!!! stop changing this in lua code when writing vim. ...
iabbrev Im I'm
iabbrev IM I'm
" Abolish ive I've " keeps annoying me when writing five
iabbrev Ill I'll
iabbrev youll you'll
" id is a word but I would never spell it Id by itself
iabbrev Id I'd
Abolish {they,you}re {}'re
Abolish yall y'all
" flipped letters
Abolish ahve have
Abolish sohw show
Abolish blaance balance
Abolish sohuld should
Abolish tihnk think
Abolish shoudl should
Abolish udnerstand understand
Abolish palce place
" spelling mistakes
Abolish {despa,sepe}rat{e,es,ed,ing,ely,ion,ions,or} {despe,sepa}rat{}
Abolish flourescent{,ly} fluorescent{}
Abolish eucledian Euclidean
Abolish lifes lives
" short forms
Abolish algo{,s} algorithm{}
Abolish tho though
Abolish altho although
Abolish eventho even though
Abolish inspite in spite
Abolish defacto de facto
Abolish thru through
Abolish passthru passthrough
Abolish probs probably
Abolish ppl people
Abolish dialog dialogue
Abolish avail available
Abolish bc because
Abolish melb Melbourne
" TODO: only do some abbrevs in regular text
" Abolish intro introduction
Abolish prio priority
Abolish prios priorities
" auto-capitalization
iabbrev english English
iabbrev danish Danish
" This one might be annoying for writing dialogue?
Abolish combo{,s} combination{}
Abolish hi{,e}ra{,r}ch{y,ical} hi{e}{r}arch{}
Abolish noone no one
" i..e isn't valid as a keyword, so we have a snippet in luasnippets/all.lua 
" for that kind of typo 
Abolish ie i.e.
Abolish unqiue unique
Abolish simplicies simplices
Abolish pertruding protruding
Abolish effecient efficient
Abolish persuit pursuit
Abolish oc{,c}uring occurring
Abolish feasab{ility,le} feasib{ility,le}
Abolish preceed{,ed,ing} preced{,ed,ing}
Abolish embarass{,ing,ingly} embarrass{,ing,ingly}
Abolish corespond{,s,ing} correspond{,s,ing}
Abolish discernable discernible
Abolish wildtype wild type

" language specific, see lua/keymap
" TODO: make this a bit more convenient, for autocorrecting language specific 
" spelling errors in prose only
function s:ToggleDanskAbbrev() abort
    if &iminsert
        " echom 'Dansk abb'
        " Dansk
        abbrev feks f.eks.
        silent! unabbrev eg
        silent! unabbrev Eg
        " only works with <buffer>
        silent! iunabbrev <buffer> ti
        silent! iunabbrev <buffer> i
    else
        " echom 'English abb'
        " English
        " danish unabbrevs has to be silent! since they might not have been 
        " set yet so will give an error message
        silent! iunabbrev <buffer> feks
        iabbrev eg e.g.
        iabbrev Eg E.g.
        if &ft == 'asciidoc'
            " for regular text where we wouldn't be talking about a variable i 
            " or in Danish where i is a word.
            iabbrev <buffer> i I
            iabbrev <buffer> ti it
        endif
    endif
endfunction
call s:ToggleDanskAbbrev()
augroup ToggleDanskAbbrev
  autocmd User ToggleDansk :call s:ToggleDanskAbbrev()
  autocmd FileType * :call s:ToggleDanskAbbrev()
augroup END

