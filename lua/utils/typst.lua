local util = require "utils/init"

local M = {}

--- Typst rules that set the fonts of `look`, for a document's start. As in
--- LaTeX with unicode-math, operator names and quoted text in math use the text
--- font, and letters, symbols and numbers the math font. Typst selects a face
--- by weight and style, not by name, so each face is chosen by its weight.
--- @param look MathLook
--- @return string typ
function M.fonts_typ(look)
    local text = look.text
    -- Math letters are `symbol` elements, and numbers and quoted text, also
    -- inside `op`, are `text`. The math font is the fallback for the text font,
    -- as `bold("…")` maps letters to math alphanumerics. The equation show-set
    -- rule sets the weight, as typst equations default to weight 450.
    return util.fill([[
#set text(font: "{{text}}", weight: {{regular}})
#set strong(delta: 0)
#show emph: set text(weight: {{italic}})
#show strong: set text(weight: {{bold}})
#show emph: it => { show strong: set text(weight: {{bold_italic}}); it }
#show strong: it => { show emph: set text(weight: {{bold_italic}}); it }
#show math.equation: set text(font: "{{math}}", weight: {{regular}})
#show math.equation: it => {
  show text: t => if t.text.match(regex("^[\d.,]+$")) == none { set text(font: ("{{text}}", "{{math}}")); t } else { t }
  it
}]], {
        text = text.family,
        math = look.math.family,
        regular = text.regular.weight,
        italic = text.italic.weight,
        bold = text.bold.weight,
        bold_italic = text.bold_italic.weight,
    })
end

--- Typst rules that fit math to a terminal cell `box_in` tall, for after the
--- fonts are set:
--- - The text size is set so that the math x-height, cap-height and descender
---   depth best match `t`. Each target i gives a cell length C_i = m_i / t_i,
---   where t_i is `t.x`, `t.cap` or `t.desc`, and m_i is the ink height of `x`
---   or `H` above the baseline or the ink depth of `p`, in the math font at the
---   current size. The fit scales C = sqrt(min C_i · max C_i) to `box_in`,
---   which makes the largest over-size and under-size errors equal and
---   opposite in log.
--- - Each equation starts with an invisible strut `box_in` tall,
---   `t.below · box_in` of it below the baseline, so it extends at least that
---   far. Inline equations extend to the union of the strut and their ink only
---   with snacks' `top-edge: "bounds", bottom-edge: "bounds"` text edges.
--- - The page margin is 0, so an inline equation's page is its cell box.
--- @param t {x: number, cap: number, desc: number, below: number} fractions of the cell height: x-height, cap-height, descender depth, baseline to cell bottom
--- @param box_in number cell height in inches
--- @return string typ
function M.cell_fit_typ(t, box_in)
    -- The measurements come before the strut rule, or they would see the strut.
    -- The strut's class "opening" adds no spacing, like the start of an
    -- equation. The label stops the rebuilt equation from matching its own show
    -- rule again.
    return util.fill([[
#set page(margin: 0pt)
#show: body => context {
  let ink(eq, bottom-edge) = measure({ show math.equation: set text(top-edge: "bounds", bottom-edge: bottom-edge); eq }).height
  let cells = (
    ink($x$, "baseline") / {{x}},
    ink($H$, "baseline") / {{cap}},
    (ink($p$, "bounds") - ink($p$, "baseline")) / {{desc}},
  )
  let cell = calc.sqrt(calc.min(..cells).pt() * calc.max(..cells).pt()) * 1pt
  set text(size: text.size * ({{box}}in / cell))
  show math.equation: it => if it.at("label", default: none) == <math-cell-strut> { it } else {
    [#math.equation(block: it.block, math.class("opening", box(height: {{box}}in, baseline: {{below}}in)) + it.body) <math-cell-strut>]
  }
  body
}]], { x = t.x, cap = t.cap, desc = t.desc, box = box_in, below = t.below * box_in })
end

--- Checks whether typst finds a font family, without blocking.
--- @param family string
--- @param on_result fun(found: boolean, stderr: string) called in the main loop
function M.has_font(family, on_result)
    vim.system({ "typst", "fonts" }, { text = true }, vim.schedule_wrap(function(res)
        on_result(vim.list_contains(vim.split(res.stdout or "", "\n"), family), res.stderr or "")
    end))
end

return M
