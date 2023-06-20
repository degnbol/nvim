" insert real tab by default
setlocal noexpandtab
" there was some annoying issue with indents inserted in regular text
set nosmartindent
" lists doesn't start with number on first column but uses \item
set fo-=n
" we will often be writing prose that wraps at textwidth columns but not code.
" Let's try with sidescrolloff so the window doesn't scroll when we get close 
" to the edge.
set sidescrolloff=0

" vimtex has a lot of nice default conceals, e.g. greek in maths, \textbf, etc.
set conceallevel=2

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
augroup foHack
    autocmd InsertCharPre *.tex if IsText() | set fo+=a | else | set fo-=a | endif
augroup end

nmap <plug>Latex2Unicode v<plug>(vimtex-a$):!$XDG_CONFIG_HOME/nvim/tex/unicode/latex2unicode.sh<CR>
nmap <plug>Unicode2Latex v<plug>(vimtex-a$):!$XDG_CONFIG_HOME/nvim/tex/unicode/unicode2latex.sh<CR>
xmap <plug>Latex2Unicode_visual :!$XDG_CONFIG_HOME/nvim/tex/unicode/latex2unicode.sh<CR>
xmap <plug>Unicode2Latex_visual :!$XDG_CONFIG_HOME/nvim/tex/unicode/unicode2latex.sh<CR>
" doesn't seem to work
" silent! call repeat#set("\<Plug>Latex2Unicode", v:count)
" silent! call repeat#set("\<Plug>Unicode2Latex", v:count)
" convenient macros stored to registers u and l that goes to next math then 
" converts. Can be repeated with a count, and recalled quickly with @@ to 
" quickly convert each math env in a file.
" ']4' and '<leader>lu' are defined in whichkey.lua
let @u=']4 lu'
let @l=']4 lU'

