; extends

; The leading ':' of an ex command gets prompt colour — mirroring prompts
; elsewhere. Applies in the ui2 cmdline and in injected vim code blocks. Not
; the ',' inside a range (a different node) nor the ':' in a dict literal
; (nested deeper). The ERROR case covers half-typed commands (':', ':lua',
; ':set'), where the colon parses as the first child of an ERROR node.
(script_file (":") @punctuation.special.prompt)
(ERROR . (":") @punctuation.special.prompt)
