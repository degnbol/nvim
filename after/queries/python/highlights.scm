;extends
"in" @repeat

; I convert jupyter notebooks to python when opening them, and the cell 
; boundaries are converted to a comment # %% so I want to highlight those 
; boundaries.
(((comment) @cell)
  (#lua-match? @cell "# %%%%")) ; %%%% = two percent signs

