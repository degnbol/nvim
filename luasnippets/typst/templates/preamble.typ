// New in 0.14.0
// Character-level justification
// https://typst.app/docs/reference/model/par/#parameters-justification-limits
#set par(justify: true, justification-limits: (
    spacing: (min: 100% * 2 / 3, max: 150%),
    tracking: (min: -0.01em, max: 0.02em),
))
