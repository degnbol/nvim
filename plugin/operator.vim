" operators are the function one can apply to motions, e.g. d is a delete
" operator that can be applied to motions such as ap.
" This allows for making custom operators. New motions can be made with
" operator-append mapping

function Operator(func)
    let &operatorfunc = a:func
    return 'g@'
endfunction

" In lua:
" function DoSomething(type, ...)
"     " Do something with &type that is e.g. line or char
"     " motion is defined by '[ and ']
" end

" function ReplOperator(type, ...)
"     exe 'vaf<Plug>(iron-visual-send)`>j'
" endfunction

" nnoremap <expr> gs Operator('v:lua.DoSomething')

