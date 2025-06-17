" have gf (goto file) work when writing the common $ROOT/PATH pattern.
setlocal includeexpr=substitute(v:fname,'\\$ROOT/','','')

" quick macros for toggling between inline and not inline functions.
" Maybe write full lua functions so we can also support begin/end notation.
" Explanation:
" $[f go to beginning of function we are inside. $ so we don't go to previous 
" function in case we are on first char.
let @i="$[fdwA = \<Esc>JJD"
" Explanation:
" f(%f= we want to go to the = that defines the function but there may be = 
" inside the function args for default values and the contents of the function 
" could be e.g. arg == something so we find it by finding open paren, jumping 
" to matching paren to jump over args, then first = should be it.
let @f="Ifunction \<Esc>f(%f=caw\<BS>\<CR>\<ESC>oend\<ESC>[f=af"
