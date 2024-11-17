" Even though treesitter has string covered, we need a regex string match for 
" contained to work.
syntax region String start='"' skip='\"' end='"'
syntax region String start="'" skip="\'" end="'"

syntax keyword Italic
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
    \ polymer
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

syntax match Italic /m\./ containedin=String contained
syntax match Italic /c\./ containedin=String contained
syntax match Italic /s\./ containedin=String contained
syntax match Italic /r\./ containedin=String contained
syntax match Italic /i\./ containedin=String contained
syntax match Italic /n\./ containedin=String contained
syntax match Italic /idx\./ containedin=String contained
syntax match Italic /ps\./ containedin=String contained
syntax match Italic /bs\./ containedin=String contained
syntax match Italic /bc\./ containedin=String contained
syntax match Italic /br\./ containedin=String contained
syntax match Italic /bca\./ containedin=String contained
syntax match Italic /bm\./ containedin=String contained
syntax match Italic /bf\./ containedin=String contained
syntax match Italic /bto\./ containedin=String contained
syntax match Italic /nbr\./ containedin=String contained
syntax match Italic /xt\./ containedin=String contained
syntax match Italic /w\./ containedin=String contained
syntax match Italic /a\./ containedin=String contained
syntax match Italic /x\./ containedin=String contained
syntax match Italic /nto\./ containedin=String contained
syntax match Italic /be\./ containedin=String contained
syntax match Italic /pc\./ containedin=String contained
syntax match Italic /fc\./ containedin=String contained
syntax match Italic /e\./ containedin=String contained
syntax match Italic /fxd\./ containedin=String contained
syntax match Italic /rst\./ containedin=String contained
syntax match Italic /msk\./ containedin=String contained
syntax match Italic /f\./ containedin=String contained
syntax match Italic /org\./ containedin=String contained
syntax match Italic /ino\./ containedin=String contained
syntax match Italic /sol\./ containedin=String contained
syntax match Italic /pol\./ containedin=String contained
syntax match Italic /polymer\.protein/ containedin=String contained
syntax match Italic /polymer\.nucleic/ containedin=String contained
syntax match Italic /h\./ containedin=String contained
syntax match Italic /bb\./ containedin=String contained
syntax match Italic /sc\./ containedin=String contained
syntax match Italic /don\./ containedin=String contained
syntax match Italic /acc\./ containedin=String contained
syntax match Italic /v\./ containedin=String contained
syntax match Italic /pr\./ containedin=String contained
syntax match Italic /nt\./ containedin=String contained

syntax match Bold "/" containedin=String contained
syntax match Bold "*" containedin=String contained
syntax match Bold "!" containedin=String contained
syntax match Bold "&" containedin=String contained
syntax match Bold "|" containedin=String contained
