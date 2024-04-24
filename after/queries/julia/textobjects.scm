;extends
; The following are for a nicer parameter textobj experience, for e.g.
; [a ]a <leader>a<left> <leader>a<right>
; I had some captures in the past that seem to have been fixed now. 
; Look in the history of this file if things start acting up again.
; We purposefully mix arguments and parameters together, 
; we don't want twice as many keybindings.
(argument_list (_) @parameter.inner)
; fix of capture e.g. 1 in [1, 2, 3]
(vector_expression (_) @parameter.inner)
; add capture e.g. 1 in [1; 2; 3] or [1 2 3]
(matrix_row (_) @parameter.inner)

