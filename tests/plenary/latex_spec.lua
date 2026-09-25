---@diagnostic disable: undefined-global
local latex = require "utils/latex"
local image_placement = require "utils/image_placement"

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

describe("cell_fit_tex", function()
    if vim.fn.executable("tectonic") == 0 or vim.fn.executable("magick") == 0 then
        pending("tectonic or magick not installed")
        return
    end
    vim.cmd.packadd "snacks.nvim"
    local snacks = require "snacks"

    -- Operator Mono 12 pt, as kitty.cell_fractions reports it, in 15×35 px cells.
    local cell_w, cell_h = 15, 35
    local operator = { x = 13.2 / 35, cap = 16.5 / 35, desc = 0.135, below = 6 / 35 }
    local density = 192
    local box_in = image_placement.row_height_in(cell_h / cell_w, density)
    snacks.image.terminal.size = function()
        return { cell_width = cell_w, cell_height = cell_h, scale = cell_w / 8 }
    end

    --- Compile math `body` in a standalone document, as snacks does for math.
    --- @param header string TeX after \begin{document}
    --- @param body string math mode content
    --- @return string pdf path
    local function compile(header, body)
        local dir = vim.fn.tempname()
        vim.fn.mkdir(dir, "p")
        local tex = table.concat({
            "\\documentclass[preview,border=0pt,varwidth,12pt]{standalone}",
            "\\usepackage{amsmath, unicode-math}",
            "\\begin{document}",
            header,
            latex.cell_fit_tex(operator, box_in),
            "{\\MathCellFit\\begin{math}\\MathCellStrut " .. body .. "\\end{math}}",
            "\\end{document}",
        }, "\n")
        vim.fn.writefile(vim.split(tex, "\n"), dir .. "/m.tex")
        local res = vim.system({ "tectonic", "-Z", "continue-on-errors", "--keep-logs", "--outdir", dir, dir .. "/m.tex" }):wait()
        assert.are.equal(0, res.code, res.stderr)
        local errors = vim.tbl_filter(function(l) return l:find("^!") end, vim.fn.readfile(dir .. "/m.log"))
        assert.are.same({}, errors)
        return dir .. "/m.pdf"
    end

    --- Rasterise page 1 of `pdf` as snacks' math convert does.
    --- @param pdf string
    --- @param dpi number
    --- @return table png `info` as snacks' identify step stores it, and `ink` rows `{top, bottom}` (0 = top edge, bottom exclusive)
    local function rasterise(pdf, dpi)
        local png = pdf:gsub("%.pdf$", dpi .. ".png")
        local res = vim.system({ "magick", "-density", tostring(dpi), pdf .. "[0]", png }):wait()
        assert.are.equal(0, res.code, res.stderr)
        local out = vim.system({ "magick", png, "-format", "%w %h %x %y %@", "info:" }, { text = true }):wait().stdout
        local w, h, x, y, ink_h, ink_y = out:match("^(%d+) (%d+) ([%d.]+) ([%d.]+) %d+x(%d+)%+%d+%+(%d+)")
        return {
            info = {
                size = { width = tonumber(w), height = tonumber(h) },
                dpi = { width = tonumber(x), height = tonumber(y) },
            },
            ink = { top = tonumber(ink_y), bottom = tonumber(ink_y) + tonumber(ink_h) },
        }
    end

    for name, header in pairs { ["Latin Modern"] = "", ["Fira Math"] = "\\setmathfont{FiraMath-Regular.otf}" } do
        describe(name, function()
            -- H has a flat base, so its ink bottom is the baseline. Italic x overshoots it.
            local pdfs = { x = compile(header, "x"), H = compile(header, "H"), p = compile(header, "p") }
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
                local png = rasterise(compile(header, "(y)x_{p}^{2}"), density)
                assert.is_true(snacks.image.util.fit("", { width = 99, height = 99 }, { info = png.info }).height <= 2)
            end)
        end)
    end
end)
