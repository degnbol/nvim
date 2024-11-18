" https://pymolwiki.org/index.php/Selection_Algebra

" Even though treesitter has string covered, we need a regex string match for 
" contained to work.
syntax region String start='"' skip='\\"' end='"'
syntax region String start="'" skip="\\'" end="'"

syntax keyword PymolStringKeyword
    \ all
    \ none
    \ enabled
    \ not
    \ and
    \ or
    \ first
    \ last
    \ model
    \ chain
    \ segi
    \ resn
    \ name
    \ resi
    \ alt
    \ index
    \ id
    \ rank
    \ pepseq
    \ label
    \ in
    \ like
    \ byobject
    \ bysegi
    \ bychain
    \ byres
    \ bycalpha
    \ bymolecule
    \ byfragment
    \ byring
    \ bycell
    \ bound_to
    \ neighbor
    \ extend
    \ within
    \ around
    \ expand
    \ gap
    \ near_to
    \ beyond
    \ partial_charge
    \ formal_charge
    \ b
    \ q
    \ ss
    \ elem
    \ bonded
    \ protected
    \ fixed
    \ restrained
    \ masked
    \ flag
    \ organic
    \ inorganic
    \ solvent
    \ guide
    \ hetatm
    \ hydrogens
    \ backbone
    \ sidechain
    \ metals
    \ donors
    \ acceptors
    \ visible
    \ rep
    \ color
    \ cartoon_color
    \ ribbon_color
    \ center
    \ origin
    \ state
    \ present
    \ x
    \ y
    \ z
    \ numeric_type containedin=String contained

syntax match PymolStringKeyword /\<m\./ containedin=String contained
syntax match PymolStringKeyword /\<c\./ containedin=String contained
syntax match PymolStringKeyword /\<s\./ containedin=String contained
syntax match PymolStringKeyword /\<r\./ containedin=String contained
syntax match PymolStringKeyword /\<i\./ containedin=String contained
syntax match PymolStringKeyword /\<n\./ containedin=String contained
syntax match PymolStringKeyword /\<idx\./ containedin=String contained
syntax match PymolStringKeyword /\<ps\./ containedin=String contained
syntax match PymolStringKeyword /\<bs\./ containedin=String contained
syntax match PymolStringKeyword /\<bc\./ containedin=String contained
syntax match PymolStringKeyword /\<br\./ containedin=String contained
syntax match PymolStringKeyword /\<bca\./ containedin=String contained
syntax match PymolStringKeyword /\<bm\./ containedin=String contained
syntax match PymolStringKeyword /\<bf\./ containedin=String contained
syntax match PymolStringKeyword /\<bto\./ containedin=String contained
syntax match PymolStringKeyword /\<nbr\./ containedin=String contained
syntax match PymolStringKeyword /\<xt\./ containedin=String contained
syntax match PymolStringKeyword /\<w\./ containedin=String contained
syntax match PymolStringKeyword /\<a\./ containedin=String contained
syntax match PymolStringKeyword /\<x\./ containedin=String contained
syntax match PymolStringKeyword /\<nto\./ containedin=String contained
syntax match PymolStringKeyword /\<be\./ containedin=String contained
syntax match PymolStringKeyword /\<pc\./ containedin=String contained
syntax match PymolStringKeyword /\<fc\./ containedin=String contained
syntax match PymolStringKeyword /\<e\./ containedin=String contained
syntax match PymolStringKeyword /\<fxd\./ containedin=String contained
syntax match PymolStringKeyword /\<rst\./ containedin=String contained
syntax match PymolStringKeyword /\<msk\./ containedin=String contained
syntax match PymolStringKeyword /\<f\./ containedin=String contained
syntax match PymolStringKeyword /\<org\./ containedin=String contained
syntax match PymolStringKeyword /\<ino\./ containedin=String contained
syntax match PymolStringKeyword /\<sol\./ containedin=String contained
syntax match PymolStringKeyword /\<pol\./ containedin=String contained
" Polymer matched instead of keyword to avoid it having priority over 
" polymer.protein and polymer.nucleic
syntax match PymolStringKeyword /\<polymer\>/ containedin=String contained
syntax match PymolStringKeyword /\<polymer\.protein\>/ containedin=String contained
syntax match PymolStringKeyword /\<polymer\.nucleic\>/ containedin=String contained
syntax match PymolStringKeyword /\<h\./ containedin=String contained
syntax match PymolStringKeyword /\<bb\./ containedin=String contained
syntax match PymolStringKeyword /\<sc\./ containedin=String contained
syntax match PymolStringKeyword /\<don\./ containedin=String contained
syntax match PymolStringKeyword /\<acc\./ containedin=String contained
syntax match PymolStringKeyword /\<v\./ containedin=String contained
syntax match PymolStringKeyword /\<pr\./ containedin=String contained
syntax match PymolStringKeyword /\<nt\./ containedin=String contained

syntax match PymolStringOperator "/" containedin=String contained
syntax match PymolStringOperator "*" containedin=String contained
syntax match PymolStringOperator "!" containedin=String contained
syntax match PymolStringOperator "&" containedin=String contained
syntax match PymolStringOperator "|" containedin=String contained
" prefix for object/selection. Difference from no prefix is that this is explicitly an object/selection, so can be named the same as any operator.
syntax match PymolStringOperator "%" containedin=String contained
" prefix for object/selection where empty selection is used if object/selection doesn't exist.
syntax match PymolStringOperator "?" containedin=String contained
syntax match PymolStringDelimiter "(" containedin=String contained
syntax match PymolStringDelimiter ")" containedin=String contained

" We can give them colour that is inbetween e.g.
" string with #7AA4A1
" keyword with #AD5C7C
" Resulting in #94808F
" However, treesitter hl takes priority as usual so we would need to disable 
" that, then only use regex hi for string and also make sure that priority 
" works. Otherwise, bold and italics are great by themselves.
hi def link PymolStringKeyword Italic
hi def link PymolStringOperator Bold
hi def link PymolStringDelimiter Bold
