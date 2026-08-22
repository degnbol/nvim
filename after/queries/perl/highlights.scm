;extends

; An interpolation inside a pattern is perl, not regex: `$foo` in `s/$foo/x/`
; would otherwise read as the regex end-anchor plus three literals, because
; after/queries/perl/injections.scm covers the whole pattern and the injected
; layer wins at equal priority. Re-assert the base @variable at 101 (> injected
; default 100) over the interpolation's own bytes only, so the rest of the
; pattern keeps its regex colours.
;
; Element access (`$h{k}`, `$a[0]`, `$re->{k}`) is flattened to one @variable
; rather than mirroring the base query's @variable.member / @punctuation.special
; split, which would take a pattern per access shape.
((regexp_content
   [
     (scalar)
     (array)
     (hash)
     (glob)
     (hash_element_expression)
     (array_element_expression)
   ] @variable)
  (#set! priority 101))
