---@diagnostic disable: undefined-global
local util = require "utils/init"
local latex = require "utils/latex"
local typst = require "utils/typst"
local image_placement = require "utils/image_placement"
local look = require "math_look"

local render = dofile("tests/math_render.lua")

for _, exe in ipairs { "typst", "tectonic", "magick", "pdffonts" } do
    if vim.fn.executable(exe) == 0 then
        pending(exe .. " not installed")
        return
    end
end
local has_font ---@type boolean|nil
typst.has_font(look.math.family, function(found) has_font = found end)
vim.wait(10000, function() return has_font ~= nil end)
if not has_font then
    pending(look.math.family .. " not installed")
    return
end
vim.cmd.packadd "snacks.nvim"
local snacks = require "snacks"
local tpl = require("snacks.picker.util").tpl

local operator = render.fractions
local cell_w, cell_h = render.cell.width, render.cell.height
local density = 192
local box_in = image_placement.row_height_in(cell_h / cell_w, density)
snacks.image.terminal.size = function()
    return { cell_width = cell_w, cell_height = cell_h, scale = cell_w / 8 }
end

--- Compile typst math `content` in snacks' typst template, with the fonts and
--- the cell fit.
--- @param content string math with its `$` delimiters
--- @return string pdf path
local function compile(content)
    local doc = util.insert_around(snacks.image.config.math.typst.tpl, "${header}",
        typst.fonts_typ(look) .. "\n", "\n" .. typst.cell_fit_typ(operator, box_in))
    return render.compile_typ(tpl(doc, { color = "#000000", header = "", content = content }, { indent = true, prefix = "$" }))
end

describe("fonts_typ", function()
    for content, fonts in pairs {
        ["$sin$"] = { "OperatorMonoSSmLig-Medium" },
        ['$"a 1"$'] = { "OperatorMonoSSmLig-Medium" },
        ["$123$"] = { "TeXGyreDejaVuMath-Regular" },
        ["$x$"] = { "TeXGyreDejaVuMath-Regular" },
        ['$bold("x")$'] = { "TeXGyreDejaVuMath-Regular" },
        -- The spaces are in the regular face.
        ["_a_ *b* *_c_* _*d*_"] = {
            "OperatorMonoSSmLig-Bold", "OperatorMonoSSmLig-BoldItalic", "OperatorMonoSSmLig-BookItalic", "OperatorMonoSSmLig-Medium",
        },
    } do
        it("draws " .. content .. " in " .. table.concat(fonts, ", "), function()
            assert.are.same(fonts, render.pdf_fonts(compile(content)))
        end)
    end
end)

describe("cell_fit_typ", function()
    -- H has a flat base, so its ink bottom is the baseline.
    local pdfs = { x = compile("$x$"), H = compile("$H$"), p = compile("$p$") }
    local H = render.rasterise(pdfs.H, density)
    local box_px = math.floor(density * box_in)

    it("shows $H$ one cell tall", function()
        local fit = snacks.image.util.fit("", { width = 99, height = 99 }, { info = H.info })
        assert.are.equal(1, fit.height)
    end)

    it("makes the PNG the snapped box", function()
        assert.are.equal(box_px, H.info.size.height)
    end)

    it("sits $H$ on the cell baseline", function()
        local h = H.info.size.height
        assert.is_true(math.abs(h - H.ink.bottom - operator.below * h) <= 1)
    end)

    it("balances x-height, cap-height and descender errors", function()
        local hi = density * 10
        local glyphs = vim.tbl_map(function(pdf) return render.rasterise(pdf, hi) end, pdfs)
        local h = glyphs.H.info.size.height
        local baseline = glyphs.H.ink.bottom
        local ratios = {
            (baseline - glyphs.x.ink.top) / (operator.x * h),
            (baseline - glyphs.H.ink.top) / (operator.cap * h),
            (glyphs.p.ink.bottom - baseline) / (operator.desc * h),
        }
        assert.is_true(math.abs(math.max(unpack(ratios)) * math.min(unpack(ratios)) - 1) < 0.02,
            vim.inspect(ratios))
    end)

    for _, content in ipairs { "$x_i$", "$p_q$" } do
        it("grows the PNG for " .. content, function()
            local png = render.rasterise(compile(content), density)
            assert.is_true(png.info.size.height > box_px)
        end)
    end

    it("grows the PNG above the cell for $H^(2^(2^2))$", function()
        -- Only the superscripts extend beyond the cell, so the ink bottom stays the baseline.
        local png = render.rasterise(compile("$H^(2^(2^2))$"), density)
        local h = png.info.size.height
        assert.is_true(h > box_px)
        assert.is_true(math.abs(h - png.ink.bottom - operator.below * box_px) <= 1)
    end)

    it("floors display math to one cell", function()
        local png = render.rasterise(compile("$ x $"), density)
        assert.are.equal(box_px, png.info.size.height)
    end)

    it("matches LaTeX's $x$ and $H$ ink", function()
        local hi = density * 10
        for _, body in ipairs { "x", "H" } do
            local tex = render.rasterise(render.compile_tex(render.tex_math(latex.fonts_tex(look), body, box_in)), hi).ink
            local typ = render.rasterise(pdfs[body], hi).ink
            assert.is_true(math.abs((tex.right - tex.left) - (typ.right - typ.left)) <= 1, body .. " width")
            assert.is_true(math.abs((tex.bottom - tex.top) - (typ.bottom - typ.top)) <= 1, body .. " height")
        end
    end)
end)

describe("str_literal", function()
    it("round-trips \\, \", a newline and a tab", function()
        local s = 'a\\b "c"\nd\te'
        local res = vim.system({ "typst", "eval", typst.str_literal(s), "--format", "json" }, { text = true }):wait()
        assert.are.equal(0, res.code, res.stderr)
        assert.are.equal(s, vim.json.decode(res.stdout))
    end)
end)

describe("mitex_typ", function()
    local density_hi = density * 10
    for latex_body, typ in pairs {
        x = "$x$",
        H = "$H$",
        ["\\frac{k_{cat}}{K_M}"] = "$(k_(c a t))/K_M$",
    } do
        it("renders inline " .. latex_body .. " as " .. typ, function()
            local mitex = render.rasterise(compile(typst.mitex_typ(latex_body, true)), density_hi)
            local native = render.rasterise(compile(typ), density_hi)
            assert.are.same(native.info.size, mitex.info.size)
            assert.are.same(native.ink, mitex.ink)
        end)
    end

    it("renders display aligned taller than one cell", function()
        local body = "\\begin{aligned} a &= b \\\\ c &= d \\end{aligned}"
        local png = render.rasterise(compile(typst.mitex_typ(body, false)), density)
        assert.is_true(png.info.size.height > math.floor(density * box_in))
    end)
end)
