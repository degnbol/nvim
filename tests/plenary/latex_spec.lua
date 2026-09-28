---@diagnostic disable: undefined-global
local latex = require "utils/latex"
local image_placement = require "utils/image_placement"
local math_look = require "math_look"
local render = dofile("tests/math_render.lua")

describe("is_inline", function()
    for _, src in ipairs { "$x$", "$`x`$", "\\(x\\)", " $x$ " } do
        it("is true for " .. src, function()
            assert.is_true(latex.is_inline(src))
        end)
    end

    for _, src in ipairs { "$$x$$", "\\[x\\]" } do
        it("is false for " .. src, function()
            assert.is_false(latex.is_inline(src))
        end)
    end
end)

describe("math_body", function()
    it("strips $", function()
        assert.are.equal("x", latex.math_body("$x$"))
    end)

    it("keeps whitespace inside $$", function()
        assert.are.equal(" x ", latex.math_body("$$ x $$"))
    end)

    it("strips \\[ \\]", function()
        assert.are.equal("x", latex.math_body("\\[x\\]"))
    end)

    it("keeps the newline before an environment", function()
        local body = latex.math_body("$$\n\\begin{aligned}x\\end{aligned}\n$$")
        assert.is_truthy(body:find("^\n\\begin{aligned}"))
    end)

    it("starts with \\begin when the environment follows the delimiter", function()
        local body = latex.math_body("$$\\begin{align}x\\end{align}$$")
        assert.is_truthy(body:find("^\\begin"))
    end)

    it("strips a backtick inside $", function()
        assert.are.equal("x", latex.math_body("$`x`$"))
    end)
end)

describe("fonts_tex", function()
    if vim.fn.executable("tectonic") == 0 or vim.fn.executable("pdffonts") == 0 then
        pending("tectonic or pdffonts not installed")
        return
    end

    it("selects each face by name", function()
        local pdf = render.compile_tex(table.concat({
            "\\documentclass{standalone}",
            "\\usepackage{amsmath, unicode-math}",
            latex.fonts_tex(math_look),
            "\\begin{document}",
            "$x \\text{a \\textit{b} \\textbf{c} \\textbf{\\textit{d}}}$",
            "\\end{document}",
        }, "\n"))
        assert.are.same({
            "OperatorMonoSSmLig-Bold", "OperatorMonoSSmLig-BoldItalic", "OperatorMonoSSmLig-BookItalic",
            "OperatorMonoSSmLig-Medium", "TeXGyreDejaVuMath-Regular",
        }, render.pdf_fonts(pdf))
    end)
end)

describe("cell_fit_tex", function()
    if vim.fn.executable("tectonic") == 0 or vim.fn.executable("magick") == 0 then
        pending("tectonic or magick not installed")
        return
    end
    vim.cmd.packadd "snacks.nvim"
    local snacks = require "snacks"
    local rasterise = render.rasterise

    local operator = render.fractions
    local cell_w, cell_h = render.cell.width, render.cell.height
    local density = 192
    local box_in = image_placement.row_height_in(cell_h / cell_w, density)
    snacks.image.terminal.size = function()
        return { cell_width = cell_w, cell_height = cell_h, scale = cell_w / 8 }
    end

    --- Compile inline math `body`, fit to the cell unless `fit` is false.
    --- @param fonts string preamble TeX that sets the fonts
    --- @param body string math mode content
    --- @param fit boolean|nil
    --- @return string pdf path
    local function compile(fonts, body, fit)
        return render.compile_tex(render.tex_math(fonts, body, fit ~= false and box_in or nil))
    end

    for name, fonts in pairs {
        ["Latin Modern"] = "\\setmathfont{latinmodern-math.otf}",
        ["Fira Math"] = "\\setmathfont{FiraMath-Regular.otf}",
        ["math_look"] = latex.fonts_tex(math_look),
        -- As a document header's font after the template's.
        ["math_look, then Fira Math"] = latex.fonts_tex(math_look) .. "\n\\setmathfont{FiraMath-Regular.otf}",
    } do
        describe(name, function()
            -- H has a flat base, so its ink bottom is the baseline. Italic x overshoots it.
            local pdfs = { x = compile(fonts, "x"), H = compile(fonts, "H"), p = compile(fonts, "p") }
            local H = rasterise(pdfs.H, density)

            it("shows $H$ one cell tall", function()
                local fit = snacks.image.util.fit("", { width = 99, height = 99 }, { info = H.info })
                assert.are.equal(1, fit.height)
            end)

            it("makes the PNG the snapped box", function()
                assert.are.equal(math.floor(density * box_in), H.info.size.height)
            end)

            it("sits $H$ on the cell baseline", function()
                local h = H.info.size.height
                assert.is_true(math.abs(h - H.ink.bottom - operator.below * h) <= 1)
            end)

            it("balances x-height, cap-height and descender errors", function()
                -- At 10× density, so px rounding is below the tolerance.
                local hi = density * 10
                local glyphs = vim.tbl_map(function(pdf) return rasterise(pdf, hi) end, pdfs)
                local h = glyphs.H.info.size.height
                local baseline = glyphs.H.ink.bottom
                -- Rendered size over target size, per metric.
                local ratios = {
                    (baseline - glyphs.x.ink.top) / (operator.x * h),
                    (baseline - glyphs.H.ink.top) / (operator.cap * h),
                    (glyphs.p.ink.bottom - baseline) / (operator.desc * h),
                }
                assert.is_true(math.abs(math.max(unpack(ratios)) * math.min(unpack(ratios)) - 1) < 0.02,
                    vim.inspect(ratios))
            end)

            it("collapses $(y)x_{p}^{2}$ to one row", function()
                local png = rasterise(compile(fonts, "(y)x_{p}^{2}"), density)
                assert.is_true(snacks.image.util.fit("", { width = 99, height = 99 }, { info = png.info }).height <= 2)
            end)

            it("grows the PNG for $x_i$", function()
                local png = rasterise(compile(fonts, "x_i"), density)
                assert.is_true(png.info.size.height > H.info.size.height)
            end)

            it("keeps the text-style $H$ at the fitted size", function()
                -- A script-style glyph has another aspect ratio.
                local hi = density * 10
                local fitted = rasterise(pdfs.H, hi).ink
                local normal = rasterise(compile(fonts, "H", false), hi).ink
                local aspect = function(ink) return (ink.right - ink.left) / (ink.bottom - ink.top) end
                assert.is_true(math.abs(aspect(fitted) / aspect(normal) - 1) < 0.01,
                    vim.inspect { fitted = aspect(fitted), normal = aspect(normal) })
            end)
        end)
    end

    it("keeps the last math font", function()
        local fonts = latex.fonts_tex(math_look) .. "\n\\setmathfont{FiraMath-Regular.otf}"
        assert.are.same({ "FiraMath-Regular" }, render.pdf_fonts(compile(fonts, "x")))
    end)
end)
