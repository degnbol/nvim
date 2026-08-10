-- The DocModel: a parsed docstring as ordered, language-agnostic sections.
-- Styles (`style/<style>.lua`) parse prose into this shape and render it back;
-- reconcile transforms it against a signature. Anything a style does not model
-- is preserved verbatim (`other`) so no hand-written content is ever lost.

local M = {}

--- One documented parameter. `type` is nil when the docstring omits it (this
--- config's Google style is typeless — the annotated signature carries types).
--- `description` may be multiline, stored newline-joined and dedented.
--- @class docstring.Param
--- @field name string
--- @field type? string
--- @field description string

--- The single, nameless return slot. `type` is refreshed from the return
--- annotation by reconcile; nil when absent.
--- @class docstring.Return
--- @field type? string
--- @field description string

--- One documented exception. `type` is the exception class; there is no
--- signature source for these, so reconcile passes them through verbatim.
--- @class docstring.Raise
--- @field type string
--- @field description string

--- A section the style does not model (e.g. `Example:`, `Note:`), kept verbatim
--- in place. `header` is the text before the colon; `body` is the raw indented
--- lines beneath it, dedented to the section-body column.
--- @class docstring.Section
--- @field header string
--- @field body string

--- A parsed docstring.
--- @class docstring.Model
--- @field summary string  free text before the first section (may be multiline)
--- @field params docstring.Param[]  documented signature params, in doc order
--- @field returns? docstring.Return
--- @field raises docstring.Raise[]
--- @field unmatched docstring.Param[]  entries with no matching signature param
--- @field other docstring.Section[]  unmodeled sections, in encounter order

--- An empty model — the starting point for fresh generation.
--- @return docstring.Model model
function M.empty()
    --- @type docstring.Model
    local model = {
        summary = "",
        params = {},
        returns = nil,
        raises = {},
        unmatched = {},
        other = {},
    }
    return model
end

return M
