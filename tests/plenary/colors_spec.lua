---@diagnostic disable: undefined-global
local colors = require("utils/colors")
local RGB = colors.RGB
local HSL = colors.HSL
local Lab = colors.Lab

-- Helper for approximate float comparisons.
local function approx_equal(a, b, tolerance)
    tolerance = tolerance or 0.001
    return math.abs(a - b) < tolerance
end

-- Assert two tables are approximately equal element-wise.
local function assert_approx_equal(actual, expected, tolerance, msg)
    tolerance = tolerance or 0.001
    for i = 1, #expected do
        assert(
            approx_equal(actual[i], expected[i], tolerance),
            string.format(
                "%s: index %d: expected %.6f, got %.6f (tolerance %.6f)",
                msg or "mismatch",
                i,
                expected[i],
                actual[i],
                tolerance
            )
        )
    end
end

describe("colors", function()
    describe("hex_to_sRGB / sRGB_to_hex round-trip", function()
        local test_cases = {
            { hex = "#000000", rgb = { 0, 0, 0 }, name = "black" },
            { hex = "#FFFFFF", rgb = { 255, 255, 255 }, name = "white" },
            { hex = "#FF0000", rgb = { 255, 0, 0 }, name = "red" },
            { hex = "#00FF00", rgb = { 0, 255, 0 }, name = "green" },
            { hex = "#0000FF", rgb = { 0, 0, 255 }, name = "blue" },
        }

        for _, tc in ipairs(test_cases) do
            it("converts " .. tc.name .. " hex to sRGB", function()
                local result = colors.hex_to_sRGB(tc.hex)
                assert.are.equal(tc.rgb[1], result.r)
                assert.are.equal(tc.rgb[2], result.g)
                assert.are.equal(tc.rgb[3], result.b)
            end)

            it("converts " .. tc.name .. " sRGB to hex", function()
                local result = colors.sRGB_to_hex(RGB(tc.rgb))
                assert.are.equal(tc.hex, result)
            end)

            it("round-trips " .. tc.name .. " through hex -> sRGB -> hex", function()
                local rgb = colors.hex_to_sRGB(tc.hex)
                local hex = colors.sRGB_to_hex(rgb)
                assert.are.equal(tc.hex, hex)
            end)
        end
    end)

    describe("HSL_to_sRGB / sRGB_to_HSL round-trip", function()
        it("converts black HSL(0,0,0) to RGB(0,0,0)", function()
            local result = colors.HSL_to_sRGB(HSL { 0, 0, 0 })
            assert.are.equal(0, result.r)
            assert.are.equal(0, result.g)
            assert.are.equal(0, result.b)
        end)

        it("converts white HSL(0,0,1) to RGB(255,255,255)", function()
            local result = colors.HSL_to_sRGB(HSL { 0, 0, 1 })
            assert.are.equal(255, result.r)
            assert.are.equal(255, result.g)
            assert.are.equal(255, result.b)
        end)

        it("converts red HSL(0,1,0.5) to RGB(255,0,0)", function()
            local result = colors.HSL_to_sRGB(HSL { 0, 1, 0.5 })
            assert.are.equal(255, result.r)
            assert.are.equal(0, result.g)
            assert.are.equal(0, result.b)
        end)

        it("achromatic greys have s=0", function()
            -- Test several grey values.
            local greys = { 0, 64, 128, 192, 255 }
            for _, v in ipairs(greys) do
                local hsl = colors.sRGB_to_HSL(RGB { v, v, v })
                assert.are.equal(0, hsl.s)
            end
        end)

        it("round-trips red through sRGB -> HSL -> sRGB", function()
            local original = RGB { 255, 0, 0 }
            local hsl = colors.sRGB_to_HSL(original)
            local result = colors.HSL_to_sRGB(hsl)
            assert.are.equal(original.r, result.r)
            assert.are.equal(original.g, result.g)
            assert.are.equal(original.b, result.b)
        end)

        it("round-trips green through sRGB -> HSL -> sRGB", function()
            local original = RGB { 0, 255, 0 }
            local hsl = colors.sRGB_to_HSL(original)
            local result = colors.HSL_to_sRGB(hsl)
            assert.are.equal(original.r, result.r)
            assert.are.equal(original.g, result.g)
            assert.are.equal(original.b, result.b)
        end)

        it("round-trips blue through sRGB -> HSL -> sRGB", function()
            local original = RGB { 0, 0, 255 }
            local hsl = colors.sRGB_to_HSL(original)
            local result = colors.HSL_to_sRGB(hsl)
            assert.are.equal(original.r, result.r)
            assert.are.equal(original.g, result.g)
            assert.are.equal(original.b, result.b)
        end)

        it("round-trips mid-grey through sRGB -> HSL -> sRGB", function()
            local original = RGB { 128, 128, 128 }
            local hsl = colors.sRGB_to_HSL(original)
            local result = colors.HSL_to_sRGB(hsl)
            assert.are.equal(original.r, result.r)
            assert.are.equal(original.g, result.g)
            assert.are.equal(original.b, result.b)
        end)
    end)

    describe("sRGB_to_lRGB / lRGB_to_sRGB round-trip", function()
        it("round-trips black through sRGB -> lRGB -> sRGB", function()
            local original = RGB { 0, 0, 0 }
            local linear = colors.sRGB_to_lRGB(original)
            local result = colors.lRGB_to_sRGB(linear)
            assert.are.equal(original.r, result.r)
            assert.are.equal(original.g, result.g)
            assert.are.equal(original.b, result.b)
        end)

        it("round-trips white through sRGB -> lRGB -> sRGB", function()
            local original = RGB { 255, 255, 255 }
            local linear = colors.sRGB_to_lRGB(original)
            local result = colors.lRGB_to_sRGB(linear)
            assert.are.equal(original.r, result.r)
            assert.are.equal(original.g, result.g)
            assert.are.equal(original.b, result.b)
        end)

        it("round-trips mid-grey through sRGB -> lRGB -> sRGB", function()
            local original = RGB { 128, 128, 128 }
            local linear = colors.sRGB_to_lRGB(original)
            local result = colors.lRGB_to_sRGB(linear)
            assert.are.equal(original.r, result.r)
            assert.are.equal(original.g, result.g)
            assert.are.equal(original.b, result.b)
        end)

        it("black sRGB converts to black linear RGB", function()
            local linear = colors.sRGB_to_lRGB(RGB { 0, 0, 0 })
            assert_approx_equal(linear, { 0, 0, 0 }, 0.001, "black linear")
        end)

        it("white sRGB converts to white linear RGB", function()
            local linear = colors.sRGB_to_lRGB(RGB { 255, 255, 255 })
            assert_approx_equal(linear, { 1, 1, 1 }, 0.001, "white linear")
        end)
    end)

    describe("Color type constructors", function()
        it("RGB accepts positional arguments", function()
            local c = RGB { 255, 128, 0 }
            assert.are.equal(255, c.r)
            assert.are.equal(128, c.g)
            assert.are.equal(0, c.b)
        end)

        it("RGB accepts named arguments", function()
            local c = RGB { r = 255, g = 128, b = 0 }
            assert.are.equal(255, c.r)
            assert.are.equal(128, c.g)
            assert.are.equal(0, c.b)
        end)

        it("HSL accepts positional arguments", function()
            local c = HSL { 0.5, 0.8, 0.3 }
            assert.are.equal(0.5, c.h)
            assert.are.equal(0.8, c.s)
            assert.are.equal(0.3, c.l)
        end)

        it("HSL accepts named arguments", function()
            local c = HSL { h = 0.5, s = 0.8, l = 0.3 }
            assert.are.equal(0.5, c.h)
            assert.are.equal(0.8, c.s)
            assert.are.equal(0.3, c.l)
        end)

        it("Lab accepts positional arguments", function()
            local c = Lab { 50, 20, -30 }
            assert.are.equal(50, c.L)
            assert.are.equal(20, c.a)
            assert.are.equal(-30, c.b)
        end)

        it("Lab accepts named arguments", function()
            local c = Lab { L = 50, a = 20, b = -30 }
            assert.are.equal(50, c.L)
            assert.are.equal(20, c.a)
            assert.are.equal(-30, c.b)
        end)

        it("RGB numeric indexing works", function()
            local c = RGB { 255, 128, 0 }
            assert.are.equal(255, c[1])
            assert.are.equal(128, c[2])
            assert.are.equal(0, c[3])
        end)
    end)

    describe("mix2_lab", function()
        local lab1 = Lab { 20, 10, -5 }
        local lab2 = Lab { 80, -20, 30 }

        it("strength=0 returns first colour", function()
            local result = colors.mix2_lab(lab1, lab2, 0)
            assert_approx_equal(result, { lab1.L, lab1.a, lab1.b }, 0.001, "strength=0")
        end)

        it("strength=1 returns second colour", function()
            local result = colors.mix2_lab(lab1, lab2, 1)
            assert_approx_equal(result, { lab2.L, lab2.a, lab2.b }, 0.001, "strength=1")
        end)

        it("strength=0.5 returns midpoint", function()
            local result = colors.mix2_lab(lab1, lab2, 0.5)
            local expected_L = (lab1.L + lab2.L) / 2
            local expected_a = (lab1.a + lab2.a) / 2
            local expected_b = (lab1.b + lab2.b) / 2
            assert_approx_equal(result, { expected_L, expected_a, expected_b }, 0.001, "strength=0.5")
        end)

        it("default strength is 0.5", function()
            local result = colors.mix2_lab(lab1, lab2)
            local expected_L = (lab1.L + lab2.L) / 2
            local expected_a = (lab1.a + lab2.a) / 2
            local expected_b = (lab1.b + lab2.b) / 2
            assert_approx_equal(result, { expected_L, expected_a, expected_b }, 0.001, "default strength")
        end)

        it("strength=0.25 is closer to first colour", function()
            local result = colors.mix2_lab(lab1, lab2, 0.25)
            local expected_L = lab1.L * 0.75 + lab2.L * 0.25
            local expected_a = lab1.a * 0.75 + lab2.a * 0.25
            local expected_b = lab1.b * 0.75 + lab2.b * 0.25
            assert_approx_equal(result, { expected_L, expected_a, expected_b }, 0.001, "strength=0.25")
        end)
    end)
end)
