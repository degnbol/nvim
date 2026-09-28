--- A text face of a font family.
--- @class MathLook.Face
--- @field name string fontspec face suffix after the family, e.g. "Book Italic"
--- @field weight integer OS/2 weight class, for typst

--- Fonts for rendered math.
--- @class MathLook
--- @field math {family: string, file: string} math font: family name (typst) and file name (fontspec, from tectonic's bundle)
--- @field text {family: string, regular: MathLook.Face, italic: MathLook.Face, bold: MathLook.Face, bold_italic: MathLook.Face} text font, for `\text`, operator names and quoted text

-- TeX Gyre DejaVu Math: of the full math fonts in tectonic's bundle, its
-- x-height, cap-height and descender best fit Operator Mono's. The text font is
-- the terminal font (kitty.conf `font_family`, `italic_font`, …). Typst needs
-- the math font installed (~/dotfiles/fonts/tex-gyre-math/install.sh).
--- @type MathLook
return {
    math = { family = "TeX Gyre DejaVu Math", file = "texgyredejavu-math.otf" },
    text = {
        family = "Operator Mono SSm Lig",
        regular = { name = "Medium", weight = 350 },
        italic = { name = "Book Italic", weight = 325 },
        bold = { name = "Bold", weight = 400 },
        bold_italic = { name = "Bold Italic", weight = 400 },
    },
}
