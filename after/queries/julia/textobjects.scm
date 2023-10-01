;extends
(argument_list (_) @parameter.inner)

; parameters are really what goes in e.g. function definition and arguments are 
; when they are called but we are mixing it all together as parameters for a 
; simple single textobj to call with e.g. ]a, <leader>a<left>, via etc.
; parameter_list may have the child keyword_parameters which we shouldn't 
; capture, so I'm writing out all other options. i.e. I'm trying to do all 
; minus keyword_parameters.
(parameter_list (identifier) (typed_parameter) (optional_parameter) @parameter.inner)
(keyword_parameters (_) @parameter.inner)
; capture e.g. 1 in [1, 2, 3]
(vector_expression (_) @parameter.inner)
; capture e.g. 1 in [1; 2; 3] or [1 2 3]
(matrix_row (_) @parameter.inner)
