// Glossary/acronyms/abbreviations/nomenclature
#import "@preview/glossy:0.9.0": *
// #show: init-glossary.with(yaml("glossary.yaml"))

// New in 0.14.0
// Character-level justification
// https://typst.app/docs/reference/model/par/#parameters-justification-limits
// The "spacing" controls spacing between words, and "tracking" controls spacing between characters.
// "spacing" is set to defaults and "tracking" to their suggested values, however default is zero.
#set par(justify: true, justification-limits: (
    spacing: (min: 100% * 2 / 3, max: 150%),
    tracking: (min: -0.01em, max: 0.02em),
))

// Maybe set this to have horizontal division when inline and vertical in block mode.
// That would replicate the default of latex, but not sure if needed.
// https://typst.app/docs/reference/math/frac/#parameters-style
#show math.equation.where(block: false): set math.frac(style: "horizontal")

#set document(
    title: [TITLE]
)

#title()

= Section
