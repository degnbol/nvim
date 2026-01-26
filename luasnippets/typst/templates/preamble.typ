// New in 0.14.0
// Character-level justification
// https://typst.app/docs/reference/model/par/#parameters-justification-limits
// The "spacing" controls spacing between words, and "tracking" controls spacing between characters.
// "spacing" is set to defaults and "tracking" to their suggested values, however default is zero.
#set par(justify: true, justification-limits: (
    spacing: (min: 100% * 2 / 3, max: 150%),
    tracking: (min: -0.01em, max: 0.02em),
))

// Show links written verbatum ("https://..." instead of "click here") with monospace font.
#import "@preview/linkify:0.1.1"
#import linkify.display: url-as-raw
#show: url-as-raw
// Underline links (either case)
#show link: underline

#set page(numbering: "1")

#set heading(numbering: "1.1")

// Add padding (larger page margins) on either side of caption, left-adjust, and change "Table" and "Figure" to smallcaps.
#show figure.caption: it => pad(x: 2em,
align(left)[                                                       
    #smallcaps[#it.supplement #it.counter.display(it.numbering)]#it.separator#it.body
]
)

// Maybe set this to have horizontal division when inline and vertical in block mode.
// That would replicate the default of latex, but not sure if needed.
// https://typst.app/docs/reference/math/frac/#parameters-style
#show math.equation.where(block: false): set math.frac(style: "horizontal")

// Glossary/acronyms/abbreviations/nomenclature
#import "@preview/glossy:0.9.0": *
// #show: init-glossary.with(yaml("glossary.yaml"))

// Chemistry
// Reclassify operators to "normal" to remove math-mode spacing around bonds
// #import "@preview/chemformula:0.1.1": ch as _ch
// #let ch(formula) = {
//   show sym.minus: math.class("normal", sym.minus)
//   show "=": math.class("normal", "=")
//   show ",": math.class("normal", ",")
//   _ch(formula)
// }
// #import "@preview/alchemist:0.1.8": *

#set document(
    title: [TITLE]
)

#title()

= Section
