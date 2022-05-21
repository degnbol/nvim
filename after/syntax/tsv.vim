" conceal end of long (40 chars+) entry with '...'
syn match LargeCell    /\v[^\t]{40}\zs[^\t]+/ conceal cchar=… nextgroup=LargeCellTab
" The tab following such a conceal doesn't from the conceal above.
" We have to also use the nextgroup syntax that looks for a match right after 
" LargeCell. We match any tab but since we are 'contained' we will only be 
" looked for if another group where to contain LargeCellTab so it shouldn't be 
" inefficient.
syn match LargeCellTab /\t/ conceal contained
hi clear Conceal
hi link Conceal Comment
