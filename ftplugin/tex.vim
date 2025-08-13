" Most of the time _ is for math and does not belong within a keyword.
" Other times it appears in e.g. links it has to be escaped, so wouldn't work 
" as a keyword anyways.
setlocal iskeyword-=_

" insert real tab by default
setlocal noexpandtab
" there was some annoying issue with indents inserted in regular text
setlocal nosmartindent
" lists doesn't start with number on first column but uses \item
setlocal fo-=n
" we will often be writing prose that wraps at textwidth columns but not code.
" Let's try with sidescrolloff so the window doesn't scroll when we get close 
" to the edge.
setlocal sidescrolloff=0

" vimtex has a lot of nice default conceals, e.g. greek in maths, \textbf, etc. but it moves the text too much and hides \vspace etc.
setlocal conceallevel=0

setlocal wrap

" in_mathzone and others call stack under the hood:
" https://github.com/lervag/vimtex/blob/c2f38c25375e6fb06654c3de945995c925b286e6/autoload/vimtex/syntax.vim
" It is kinda vimtex's version of treesitter I think.
" Empty when we are not in any environment. Also empty for text in Itemize 
" environment.
function IsText() abort
    if vimtex#syntax#in_mathzone()
        return 0
    endif
    let l:env = vimtex#delim#get_surrounding('env_tex')[1]
    if empty(l:env)
        return 1
    endif
    let l:name = l:env['name']
    if l:name == 'document'
        return 1
    endif
    let l:cmd = vimtex#cmd#get_current()
    if empty(l:cmd)
        return 0
    endif
    return l:cmd['name'] == "\\caption"
    " return len(vimtex#syntax#stack()) == 0
endfunction
" augroup foHack
"     autocmd InsertCharPre <buffer> if !&wrap && IsText() | set fo+=a | else | set fo-=a | endif
" augroup end
" maybe just wrap instead of this hack, although we don't wanna wrap tables

nmap <buffer> <plug>Latex2Unicode v<plug>(vimtex-a$):!$XDG_CONFIG_HOME/nvim/tex/unicode/latex2unicode.sh<CR>
nmap <buffer> <plug>Unicode2Latex v<plug>(vimtex-a$):!$XDG_CONFIG_HOME/nvim/tex/unicode/unicode2latex.sh<CR>
xmap <buffer> <plug>Latex2Unicode_visual :!$XDG_CONFIG_HOME/nvim/tex/unicode/latex2unicode.sh<CR>
xmap <buffer> <plug>Unicode2Latex_visual :!$XDG_CONFIG_HOME/nvim/tex/unicode/unicode2latex.sh<CR>
" doesn't seem to work
" silent! call repeat#set("\<Plug>Latex2Unicode", v:count)
" silent! call repeat#set("\<Plug>Unicode2Latex", v:count)
" convenient macros stored to registers u and l that goes to next math then 
" converts. Can be repeated with a count, and recalled quickly with @@ to 
" quickly convert each math env in a file.
" ']4' and '<leader>lu' are defined in whichkey.lua
let @u=']4 lu'
let @l=']4 lU'

" We aren't using em-dashes much at all in latex due to style guide suggesting 
" spaced en-dash
" https://www.stylemanual.gov.au/grammar-punctuation-and-conventions/punctuation/dashes
" abbrev correction is also very nice in that I have to add a space afterwards 
" so if I really want --- then I can still get it.
" The unicode em-dash can then be mapped with the newunicodechar package
iabbrev --- ⎯

" local to window. Some window I'm switching to sometimes must be setting it 
" so I disable it here.
setlocal signcolumn=no

