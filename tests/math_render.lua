-- Test helpers that compile and rasterise math documents.
local latex = require "utils/latex"

local M = {}

--- Kitty cell size in px for Operator Mono 12 pt.
M.cell = { width = 15, height = 35 }
--- Operator Mono 12 pt in `M.cell`, as kitty.cell_fractions reports it.
M.fractions = { x = 13.2 / 35, cap = 16.5 / 35, desc = 0.135, below = 6 / 35 }

--- Standalone LaTeX document with inline math `body`, like snacks' latex template.
--- @param fonts string preamble TeX that sets the fonts
--- @param body string math mode content
--- @param box_in number|nil cell height in inches to fit the math to, or nil for no fit
--- @return string tex
function M.tex_math(fonts, body, box_in)
    return table.concat({
        "\\documentclass[preview,border=0pt,varwidth,12pt]{standalone}",
        "\\usepackage{amsmath, unicode-math}",
        latex.record_setmathfont_tex,
        fonts,
        "\\begin{document}",
        box_in and table.concat({
            latex.cell_fit_tex(M.fractions, box_in),
            "{\\MathCellFit\\begin{math}\\MathCellStrut " .. body .. "\\end{math}}",
        }, "\n") or ("$" .. body .. "$"),
        "\\end{document}",
    }, "\n")
end

--- Compile a LaTeX document with tectonic. Fails the test on a compile error.
--- @param tex string document source
--- @return string pdf path
function M.compile_tex(tex)
    local dir = vim.fn.tempname()
    vim.fn.mkdir(dir, "p")
    vim.fn.writefile(vim.split(tex, "\n"), dir .. "/m.tex")
    local res = vim.system({ "tectonic", "-Z", "continue-on-errors", "--keep-logs", "--outdir", dir, dir .. "/m.tex" }):wait()
    assert.are.equal(0, res.code, res.stderr)
    local errors = vim.tbl_filter(function(l) return l:find("^!") end, vim.fn.readfile(dir .. "/m.log"))
    assert.are.same({}, errors)
    return dir .. "/m.pdf"
end

--- Compile a typst document. Fails the test on an error or a warning.
--- @param typ string document source
--- @return string pdf path
function M.compile_typ(typ)
    local dir = vim.fn.tempname()
    vim.fn.mkdir(dir, "p")
    vim.fn.writefile(vim.split(typ, "\n"), dir .. "/m.typ")
    local res = vim.system({ "typst", "compile", dir .. "/m.typ", dir .. "/m.pdf" }, { text = true }):wait()
    assert.are.equal(0, res.code, res.stderr)
    assert.are.equal("", res.stderr)
    return dir .. "/m.pdf"
end

--- Rasterise page 1 of `pdf` as snacks' math convert does.
--- @param pdf string
--- @param dpi number
--- @return table png `info` as snacks' identify step stores it, and the `ink` bounding box in px: `top`, `bottom`, `left`, `right` (0 = top/left edge, bottom/right exclusive)
function M.rasterise(pdf, dpi)
    local png = pdf:gsub("%.pdf$", dpi .. ".png")
    local res = vim.system({ "magick", "-density", tostring(dpi), pdf .. "[0]", png }):wait()
    assert.are.equal(0, res.code, res.stderr)
    local out = assert(vim.system({ "magick", png, "-format", "%w %h %x %y %@", "info:" }, { text = true }):wait().stdout)
    local w, h, x, y, ink_w, ink_h, ink_x, ink_y = out:match("^(%d+) (%d+) ([%d.]+) ([%d.]+) (%d+)x(%d+)%+(%d+)%+(%d+)")
    return {
        info = {
            size = { width = tonumber(w), height = tonumber(h) },
            dpi = { width = tonumber(x), height = tonumber(y) },
        },
        ink = {
            top = tonumber(ink_y),
            bottom = tonumber(ink_y) + tonumber(ink_h),
            left = tonumber(ink_x),
            right = tonumber(ink_x) + tonumber(ink_w),
        },
    }
end

--- Fonts embedded in `pdf`.
--- @param pdf string
--- @return string[] sorted PostScript names without the subset prefix
function M.pdf_fonts(pdf)
    local res = vim.system({ "pdffonts", pdf }, { text = true }):wait()
    assert.are.equal(0, res.code, res.stderr)
    local fonts = {}
    for i, line in ipairs(vim.split(vim.trim(res.stdout), "\n")) do
        if i > 2 then
            fonts[#fonts + 1] = line:match("^%S+"):gsub("^%u+%+", ""):gsub("%-Identity%-H$", "")
        end
    end
    table.sort(fonts)
    return fonts
end

return M
