;extends

; ... is captured by @lsp.type.escape.typst that links to comment.
; I don't want ... colored, so I capture it with a higher priority custom
; capture that gets linked to Normal.
((shorthand) @punct.ellipsis
  (#eq? @punct.ellipsis "...")
  (#set! "priority" 200)); higher priority than 125

; Heading marker: the base query captures "=" as @operator (bold, blue)
; because it appears in a broad operator list. Override at higher priority
; so heading markers get @markup.heading.N (dimmed, no bold) instead.
(heading ("=" @markup.heading.1 (#set! "priority" 200)))
(heading ("==" @markup.heading.2 (#set! "priority" 200)))
(heading ("===" @markup.heading.3 (#set! "priority" 200)))
(heading ("====" @markup.heading.4 (#set! "priority" 200)))
(heading ("=====" @markup.heading.5 (#set! "priority" 200)))
(heading ("======" @markup.heading.6 (#set! "priority" 200)))

; Heading text: match markdown's @text.title styling (fg + bold + underdouble).
; The base query captures the whole (heading) as @markup.heading.N (dimmed).
; This captures just the text content at higher priority.
(heading "=" (text) @text.title1 (#set! "priority" 200))
(heading "==" (text) @text.title2 (#set! "priority" 200))
(heading "===" (text) @text.title3 (#set! "priority" 200))
(heading "====" (text) @text.title4 (#set! "priority" 200))
(heading "=====" (text) @text.title5 (#set! "priority" 200))
(heading "======" (text) @text.title6 (#set! "priority" 200))

