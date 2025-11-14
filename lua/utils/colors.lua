-- Taken from
-- https://github.com/rktjmp/lush.nvim/blob/main/lua/lush/vivid/rgb/convert.lua
-- https://github.com/rktjmp/lush.nvim/blob/main/lua/lush/vivid/hsl/convert.lua
-- Which in turn adapted the code from https://github.com/EmmanuelOga/columns/blob/master/utils/color.lua
-- Implemented:
-- - Conversions between HEX, RGB, and HSL.
-- - Mixing two colors in HSL space.

local round = require "utils.math".round
local pi2 = math.pi * 2

local M = {}

---Convert from table with r, g, and b to hex string.
---Example: {r=255, g=255, b=255} -> "#FFFFFF"
---@param rgb table
---@return string
function M.rgb_to_hex(rgb)
    return string.format("#%02X%02X%02X", rgb.r, rgb.g, rgb.b)
end

---Convert from hex color string to table with r, g, and b.
---@param hex_str string
---@return table
function M.hex_to_rgb(hex_str)
    -- normalise
    local hex = "[abcdef0-9][abcdef0-9]"
    local pat = "^#(" .. hex .. ")(" .. hex .. ")(" .. hex .. ")$"
    hex_str = string.lower(hex_str)

    -- smoke test
    assert(string.find(hex_str, pat) ~= nil,
        "hex_to_rgb: invalid hex_str: " .. tostring(hex_str))

    -- convert
    local r, g, b = string.match(hex_str, pat)
    r, g, b = tonumber(r, 16), tonumber(g, 16), tonumber(b, 16)

    return { r = r, g = g, b = b }
end

---Convert HSL to RGB.
---RGB values in range [0;255] and HSL values in [0;1].
---Example {h=0, s=0, l=0} -> {r=0, g=0, b=0}
---@param hsl table
---@return table
function M.hsl_to_rgb(hsl)
    local r, g, b
    local h, s, l = hsl.h, hsl.s, hsl.l

    if s == 0 then
        r, g, b = l, l, l -- achromatic
    else
        local function hue2rgb(p, q, t)
            if t < 0 then t = t + 1 end
            if t > 1 then t = t - 1 end
            if t < 1 / 6 then return p + (q - p) * 6 * t end
            if t < 1 / 2 then return q end
            if t < 2 / 3 then return p + (q - p) * (2 / 3 - t) * 6 end
            return p
        end

        local q
        if l < 0.5 then q = l * (1 + s) else q = l + s - l * s end
        local p = 2 * l - q

        r = hue2rgb(p, q, h + 1 / 3)
        g = hue2rgb(p, q, h)
        b = hue2rgb(p, q, h - 1 / 3)
    end

    return {
        r = round(r * 255),
        g = round(g * 255),
        b = round(b * 255)
    }
end

---Convert table with r,g,b to table with h,s,l.
---RGB values in range [0;255] and HSL values in [0;1].
---@param rgb table
---@return table
function M.rgb_to_hsl(rgb)
    local r, g, b = rgb.r / 255, rgb.g / 255, rgb.b / 255

    local max, min = math.max(r, g, b), math.min(r, g, b)
    local h, s, l

    l = (max + min) / 2

    if max == min then
        h, s = 0, 0 -- achromatic
    else
        local d = max - min
        if l > 0.5 then s = d / (2 - max - min) else s = d / (max + min) end
        if max == r then
            h = (g - b) / d
            if g < b then h = h + 6 end
        elseif max == g then
            h = (b - r) / d + 2
        elseif max == b then
            h = (r - g) / d + 4
        end
        h = h / 6
    end

    return { h = h, s = s, l = l }
end

function M.hex_to_hsl(hex)
    return M.rgb_to_hsl(M.hex_to_rgb(hex))
end

function M.hsl_to_hex(hsl)
    return M.rgb_to_hex(M.hsl_to_rgb(hsl))
end

-- Mix two colors in HSL space.
-- https://stackoverflow.com/questions/35816179/calculation-algorithm-to-mix-3-hsl-colors
---@param color table {h=, s=, l=}
---@param target table {h=, s=, l=}
---@param strength number [0;1] where 0 means return `color` and 1 means `target`.
---@return table {h=, s=, l=}
function M.mix(color, target, strength)
    -- convert colors to vector
    local cv = {
        x = math.cos(color.h * pi2) * color.s,
        y = math.sin(color.h * pi2) * color.s,
        z = color.l
    }
    local tv = {
        x = math.cos(target.h * pi2) * target.s,
        y = math.sin(target.h * pi2) * target.s,
        z = target.l
    }
    -- combine
    local rv = {
        x = (cv.x * (1 - strength)) + (tv.x * strength),
        y = (cv.y * (1 - strength)) + (tv.y * strength),
        z = (cv.z * (1 - strength)) + (tv.z * strength),
    }
    -- back to color
    return {
        h = math.atan2(rv.y, rv.x) / pi2,
        s = math.sqrt(rv.x * rv.x + rv.y * rv.y),
        l = rv.z
    }
end


return M
