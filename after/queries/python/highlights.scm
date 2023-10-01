;extends
"in" @repeat

; I convert jupyter notebooks to python when opening them, and the block 
; boundaries are converted to a comment # %% so I want to highlight those 
; boundaries.
((module . (comment) @block)
  (#lua-match? @block "^# %%"))

