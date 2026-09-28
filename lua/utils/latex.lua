local M = {}

--- Whether a delimited LaTeX math source is inline (`$…$`, `$`…`$`, `\(…\)`)
--- rather than display (`$$…$$`, `\[…\]`).
--- @param src string math source including its delimiters, surrounding whitespace allowed
--- @return boolean
function M.is_inline(src)
    local s = vim.trim(src)
    return s:find("^%$[^$]") ~= nil or s:find("^\\%(") ~= nil
end

--- The body of a delimited LaTeX math source (`$…$`, `$`…`$`, `$$…$$`,
--- `\(…\)`, `\[…\]`). Whitespace outside the delimiters is trimmed, whitespace
--- inside is kept.
--- @param src string math source including its delimiters
--- @return string body
function M.math_body(src)
    local body = vim.trim(src):gsub("^%$+`?", ""):gsub("`?%$+$", "")
        :gsub("^\\[%[%(]", ""):gsub("\\[%]%)]$", "")
    return body
end

--- TeX that makes `\setmathfont` also record each call, for `cell_fit_tex` to
--- run again. Place before the first `\setmathfont`.
M.record_setmathfont_tex = [[
\ExplSyntaxOn
\tl_new:N \g__mathcellfit_setmathfont_tl
\NewCommandCopy \__mathcellfit_setmathfont:w \setmathfont
\RenewDocumentCommand \setmathfont { O{} m O{} }
  {
    \tl_gput_right:Nn \g__mathcellfit_setmathfont_tl { \__mathcellfit_setmathfont:w [#1] {#2} [#3] }
    \__mathcellfit_setmathfont:w [#1] {#2} [#3]
  }
\ExplSyntaxOff]]

--- Preamble TeX that sets the fonts of `look` with fontspec and unicode-math.
--- @param look MathLook
--- @return string tex
function M.fonts_tex(look)
    local text = look.text
    return ("\\setmainfont{%s}[UprightFont=* %s, ItalicFont=* %s, BoldFont=* %s, BoldItalicFont=* %s]\n\\setmathfont{%s}")
        :format(text.family, text.regular.name, text.italic.name, text.bold.name, text.bold_italic.name, look.math.file)
end

--- TeX that fits math to a terminal cell `box_in` tall, for a document body. Defines:
--- - `\MathCellFit`: sets the font size so that the math x-height, cap-height and
---   descender depth best match `t`. Each target i gives a cell length
---   C_i = m_i / t_i, where t_i is `t.x`, `t.cap` or `t.desc`, and m_i is the
---   height of `x`, the height of `H` or the depth of `p`, measured in the math
---   font at the current size. The fit scales C = sqrt(min C_i · max C_i) to
---   `box_in`, which makes the largest over-size and under-size errors equal
---   and opposite in log.
--- - `\MathCellStrut`: a strut `box_in` tall, `t.below · box_in` of it below the baseline.
--- Then runs the `\setmathfont` calls recorded by `record_setmathfont_tex` again
--- at the fitted size. unicode-math selects script-style glyphs by point size:
--- below the mean of the size at `\setmathfont` and its script size. It also
--- declares script sizes for that size only. Without this, text-style math at a
--- smaller fitted size would be drawn with script glyphs, and its scripts at
--- LaTeX's default ratios.
--- Place after the fonts are set, and once per document.
--- @param t {x: number, cap: number, desc: number, below: number} fractions of the cell height: x-height, cap-height, descender depth, baseline to cell bottom
--- @param box_in number cell height in inches
--- @return string tex
function M.cell_fit_tex(t, box_in)
    return ([[
\ExplSyntaxOn
\hbox_set:Nn \l_tmpa_box { $x$ }
\fp_const:Nn \c__mathcellfit_x_fp { \dim_to_fp:n { \box_ht:N \l_tmpa_box } / %.9f }
\hbox_set:Nn \l_tmpa_box { $H$ }
\fp_const:Nn \c__mathcellfit_cap_fp { \dim_to_fp:n { \box_ht:N \l_tmpa_box } / %.9f }
\hbox_set:Nn \l_tmpa_box { $p$ }
\fp_const:Nn \c__mathcellfit_desc_fp { \dim_to_fp:n { \box_dp:N \l_tmpa_box } / %.9f }
\fp_const:Nn \c__mathcellfit_size_fp
  {
    \use:c { f@size } * \dim_to_fp:n { %.9fin } / sqrt
      (
        min ( \c__mathcellfit_x_fp , \c__mathcellfit_cap_fp , \c__mathcellfit_desc_fp ) *
        max ( \c__mathcellfit_x_fp , \c__mathcellfit_cap_fp , \c__mathcellfit_desc_fp )
      )
  }
\cs_new_protected:Npx \MathCellFit
  {
    \exp_not:N \fontsize
      { \fp_use:N \c__mathcellfit_size_fp pt } { \fp_use:N \c__mathcellfit_size_fp pt }
    \exp_not:N \selectfont
  }
\cs_new_protected:Npn \MathCellStrut { \rule [ -%.9fin ] { 0pt } { %.9fin } }
\group_begin: \MathCellFit \g__mathcellfit_setmathfont_tl \group_end:
\ExplSyntaxOff]]):format(t.x, t.cap, t.desc, box_in, t.below * box_in, box_in)
end

return M
