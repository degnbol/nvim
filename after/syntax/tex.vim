" I prefer to conceal this
syn match superscript /\\textsuperscript/ conceal cchar=^
syn match subscript /\\textsubscript/ conceal cchar=_
" ... which removes their color since we are giving them a different syntax 
" group.
" We can color them as Function to match conceal color. Default is Statement. 
" We can also go with texCmdStyle which is given to \texttt which links to 
" Type.
hi link subscript texTypeStyle
hi link superscript texTypeStyle
" underline was also un-highlighted for some reason
hi link texTypeStyle texCmdStyle

syn match textgreek /\\textalpha/ conceal cchar=α
syn match textgreek /\\textbeta/ conceal cchar=β
syn match textgreek /\\textgamma/ conceal cchar=γ
syn match textgreek /\\textrho/ conceal cchar=ρ
syn match textgreek /\\textpi/ conceal cchar=π
syn match textgreek /\\textpsi/ conceal cchar=ψ
syn match textgreek /\\texttau/ conceal cchar=τ
" and so on...
" same as e.g. \AA angstrom
hi link textgreek SpecialChar

hi texItalStyle gui=italic

" redefine to add conceal (and thinrule)
syntax match texCmdBooktabs "\\\%(top\|mid\|bottom\|thin\)rule\>" conceal cchar=━

" conceal table cell
syn match texCell "\\makecell\[\a\+\]" conceal cchar=☐
hi def link texCell texCmd

" acro is matched on the whole thing, matchgroup on the start and end. 
" Nextgroup combined with the contained for the next group means we only look 
" for the autoinsert start and end patterns right after acro.
syntax region acro matchgroup=texCmdAcro start='\\[Aa]c[sl]\?p\?{' end='}' nextgroup=acroinsert
syntax region acroinsert matchgroup=Delimiter start=/\[/ end=/\]/ contained
" Not linking 'acro', and 'acroinsert', so they go uncolored, since they are 
" basically just text. texCmdAcro is colored in ftplugin/tex.lua

" hl cite key in custom function \citea, \citepdb etc. the same way as in \cite and \citeauthor
syntax match citeCustom /\\cite[A-Za-z]\+\ze[{\[]/ containedin=texCmd nextgroup=citeOpt,citeArg
syntax region citeOpt matchgroup=texDelim start='\[' end='\]' contained nextgroup=citeArg
syntax region citeArg matchgroup=texDelim start='{' end='}' contained
hi def link citeCustom texCmd " not using texCmdRef to make it clear it is custom cmd that adds more text.
hi def link citeArg texRefConcealedArg

" hl args of \import as filepaths
syntax match cmdImport /\\import\ze{/ containedin=texCmd nextgroup=argImport1
syntax region argImport1 matchgroup=texDelim start='{' end='}' contained nextgroup=argImport2
syntax region argImport2 matchgroup=texDelim start='{' end='}' contained
hi def link cmdImport texCmd
hi def link argImport1 texFileArg
hi def link argImport2 texFileArg


