-- Some parts inspired by
-- https://github.com/rktjmp/lush.nvim/blob/main/lua/lush/vivid/rgb/convert.lua
-- https://github.com/rktjmp/lush.nvim/blob/main/lua/lush/vivid/hsl/convert.lua
-- Which in turn adapted the code from https://github.com/EmmanuelOga/columns/blob/master/utils/color.lua
-- Conversions between HEX and various color spaces.
-- Other conversions based on spectal.js and
-- https://www.easyrgb.com/en/math.php
-- Mixing colors in LAB color space and also realistic color mixing from
-- https://github.com/rvanwijnen/spectral.js/blob/3.0.0/spectral.js

local round = require "utils.math".round
local mul_mat_vec = require "utils.math".mul_mat_vec
local pi2 = math.pi * 2
local C = require "utils.color_constants"

local M = {}

-- Construct new color type given keys, e.g. keys={r=1, g=2, b=3}
local function new_color_type(keys)
    local color_type = {}
    -- Allow for using rgb.r etc. notation while still being able to use rgb[0].
    color_type.__index = function(self, key)
        return self[keys[key]]
    end
    color_type.__newindex = function(self, key, value)
        self[keys[key]] = value
    end
    setmetatable(color_type, {
        -- Instantiate by calling e.g. RGB {10, 20, 30} or RGB {r=10, g=20, b=30}
        __call = function(cls, table)
            -- Check if given as integer or string indexing.
            if table[1] == nil then
                local ret = {}
                for key, i in pairs(keys) do
                    ret[i] = table[key]
                end
                return setmetatable(ret, cls)
            end
            return setmetatable(table, cls)
        end,
    })
    return color_type
end
M.RGB = new_color_type { r = 1, g = 2, b = 3 }
M.HSL = new_color_type { h = 1, s = 2, l = 3 }
M.XYZ = new_color_type { x = 1, y = 2, z = 3 }
M.Lab = new_color_type { L = 1, a = 2, b = 3 }

---Convert from table with r, g, and b to hex string.
---Example: {r=255, g=255, b=255} -> "#FFFFFF"
---@param sRGB table
---@return string
function M.sRGB_to_hex(sRGB)
    sRGB = M.RGB(sRGB)
    return string.format("#%02X%02X%02X", sRGB.r, sRGB.g, sRGB.b)
end

---Convert from hex color string to table with r, g, and b.
---@param hex_str string
---@return table
function M.hex_to_sRGB(hex_str)
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

    return M.RGB { r, g, b }
end

---Convert HSL to RGB.
---RGB values in range [0;255] and HSL values in [0;1].
---Example {h=0, s=0, l=0} -> {r=0, g=0, b=0}
---@param HSL table
---@return table
function M.HSL_to_sRGB(HSL)
    HSL = M.HSL(HSL)
    local h, s, l = HSL.h, HSL.s, HSL.l
    local r, g, b

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

    return M.RGB {
        round(r * 255),
        round(g * 255),
        round(b * 255)
    }
end

---Convert table with r,g,b to table with h,s,l.
---RGB values in range [0;255] and HSL values in [0;1].
---@param RGB table
---@return table
function M.sRGB_to_HSL(RGB)
    RGB = M.RGB(RGB)
    local r, g, b = RGB.r / 255, RGB.g / 255, RGB.b / 255

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

    return M.HSL { h, s, l }
end

---Convert from sRGB with r,g,b in range [0;255] to linear RGB in range [0;1]
---Undoes gamma-correction from an RGB-encoded color.
---https://en.wikipedia.org/wiki/SRGB#Specification_of_the_transformation
---https://stackoverflow.com/questions/596216/formula-to-determine-brightness-of-rgb-color
---@param sRGB table
---@return table lRGB
function M.sRGB_to_lRGB(sRGB)
    sRGB = M.RGB(sRGB)
    -- Normalise.
    local r = sRGB.r / 255
    local g = sRGB.g / 255
    local b = sRGB.b / 255

    -- Undo gamma-correction.
    local function uncompand(t)
        if t > 0.04045 then
            return ((t + 0.055) / 1.055) ^ C.GAMMA
        else
            return t / 12.92
        end
    end

    return M.RGB {
        uncompand(r),
        uncompand(g),
        uncompand(b),
    }
end

---Inverse of rgb_to_linear.
---https://www.easyrgb.com/en/math.php
---@param lRGB table
---@return table
function M.lRGB_to_sRGB(lRGB)
    lRGB = M.RGB(lRGB)
    local function compand(t)
        if t > 0.0031308 then
            return 1.055 * (t ^ (1 / C.GAMMA)) - 0.055
        else
            return t * 12.92
        end
    end

    return M.RGB {
        round(compand(lRGB.r) * 255),
        round(compand(lRGB.g) * 255),
        round(compand(lRGB.b) * 255),
    }
end

---Convert from linear RGB to XYZ.
---@param lRGB table
---@return table XYZ
function M.lRGB_to_XYZ(lRGB)
    return mul_mat_vec(C.CONVERSION.RGB_XYZ, M.RGB(lRGB))
end

---Inverse of lRGB_to_XYZ.
---@param XYZ table
---@return table lRGB
function M.XYZ_to_lRGB(XYZ)
    return mul_mat_vec(C.CONVERSION.XYZ_RGB, M.XYZ(XYZ))
end

---Convert from CIE XYZ to CIE Lab color space.
---Using default D65 illuminant.
---@param XYZ table {x=, y=, z=}
---@return table LAB {L=, a=, b=}
function M.XYZ_to_CIELab(XYZ)
    XYZ = M.XYZ(XYZ)

    -- Reference white point for D65 illuminant.
    local x_n = 0.95047
    local y_n = 1.00000
    local z_n = 1.08883

    local function f(t)
        if t > 0.008856 then
            return t ^ (1 / 3)
        else
            return 7.787 * t + 16 / 116
        end
    end

    local x = f(XYZ.x / x_n)
    local y = f(XYZ.y / y_n)
    local z = f(XYZ.z / z_n)

    return M.Lab {
        116 * y - 16,
        500 * (x - y),
        200 * (y - z),
    }
end

---Inverse of XYZ_to_CIELab.
---@param CIELab table
---@return table XYZ
function M.CIELab_to_XYZ(CIELab)
    CIELab = M.Lab(CIELab)

    -- Reference white point for D65 illuminant.
    local x_n = 0.95047
    local y_n = 1.00000
    local z_n = 1.08883

    local y = (CIELab.L + 16) / 116
    local x = CIELab.a / 500 + y
    local z = y - CIELab.b / 200

    local function f(t)
        if t ^ 3 > 0.00856 then
            return t ^ 3
        else
            return (t - 16 / 116) / 7.787
        end
    end

    return M.XYZ {
        f(x) * x_n,
        f(y) * y_n,
        f(z) * z_n,
    }
end

---Converts XYZ values to OKLab color space.
---@param XYZ table
---@return table OKLab
function M.XYZ_to_OKLab(XYZ)
    XYZ = M.XYZ(XYZ)
    local lms = mul_mat_vec(C.CONVERSION.XYZ_LMS, XYZ)
    -- Cube root.
    for i, v in ipairs(lms) do
        lms[i] = v ^ (1/3)
    end
    local OKLab = mul_mat_vec(C.CONVERSION.LMS_LAB, lms)
    return M.Lab(OKLab)
end


---Converts OKLab color space values to XYZ.
---@param OKLab table
---@return table XYZ
function M.OKLab_to_XYZ(OKLab)
    OKLab = M.Lab(OKLab)
    local lms = mul_mat_vec(C.CONVERSION.LAB_LMS, OKLab)
    -- Cubed.
    for i, v in ipairs(lms) do
        lms[i] = v ^ 3
    end
    local XYZ = mul_mat_vec(C.CONVERSION.LMS_XYZ, lms)
    return M.XYZ(XYZ)
end

-- Mix two colors in Lab space.
---@param Lab1 table {L=, a=, b=}
---@param Lab2 table {L=, a=, b=}
---@param strength? number [0;1] where 0 means return `color` and 1 means `target`. Default=0.5.
---@return table Lab {L=, a=, b=}
function M.mix2_lab(Lab1, Lab2, strength)
    Lab1 = M.Lab(Lab1)
    Lab2 = M.Lab(Lab2)
    -- Default equal.
    strength = strength or 0.5
    -- Lab color space is supposed be more or less perceptually uniform which I
    -- take to mean linear interpolation is fine.
    return M.Lab {
        (Lab1.L * (1 - strength)) + (Lab2.L * strength),
        (Lab1.a * (1 - strength)) + (Lab2.a * strength),
        (Lab1.b * (1 - strength)) + (Lab2.b * strength),
    }
end

-- Mix multiple colors in Lab space.
-- https://stackoverflow.com/questions/35816179/calculation-algorithm-to-mix-3-hsl-colors
---@param Labs table Each key is a LAB color table, each value is the weight. {{h=, s=, l=}=weight, ...}
---@return table {L=, a=, b=}
function M.mix_lab(Labs)
    local L = 0
    local a = 0
    local b = 0

    -- Allow for weights not summing to 1.
    local total_weight = 0
    for _, weight in pairs(Labs) do
        total_weight = total_weight + weight
    end

    for Lab, weight in pairs(Labs) do
        Lab = M.Lab(Lab)
        L = L + Lab.L * weight / total_weight
        a = a + Lab.a * weight / total_weight
        b = b + Lab.b * weight / total_weight
    end

    return M.Lab { L, a, b }
end

---Converts spectral reflectance values to XYZ using the CIE color matching functions.
---@param R table array of reflectance values.
---@return table XYZ CIE XYZ values.
function M.R_to_XYZ(R)
    assert(#R == C.R_SAMPLES, "Wrong # R samples: " .. #R)
    local XYZ = mul_mat_vec(C.CIE.CMF, R)
    return M.XYZ(XYZ)
end

---Converts linear RGB values to spectral reflectance values.
---Uses pre-calculated reflectances.
---@param lRGB table linear RGB values {r=, g=, b=} in [0;1]
---@return table R reflectance values
function M.lRGB_to_R(lRGB)
    lRGB = M.RGB(lRGB)
    local w = math.min(lRGB.r, lRGB.g, lRGB.b)

    local rw = lRGB.r - w
    local gw = lRGB.g - w
    local bw = lRGB.b - w

    local c = math.min(gw, bw);
    local m = math.min(rw, bw);
    local y = math.min(rw, gw);
    local r = math.max(0, math.min(rw - bw, rw - gw));
    local g = math.max(0, math.min(gw - bw, gw - rw));
    local b = math.max(0, math.min(bw - gw, bw - rw));

    local R = {}

    for i = 1, C.R_SAMPLES do
        R[i] = math.max(
            C.EPSILON,
            w * C.BASE_SPECTRA.W[i] + c * C.BASE_SPECTRA.C[i] + m * C.BASE_SPECTRA.M[i] + y * C.BASE_SPECTRA.Y[i] +
            r * C.BASE_SPECTRA.R[i] + g * C.BASE_SPECTRA.G[i] + b * C.BASE_SPECTRA.B[i]
        )
    end

    return R
end

---Computes the Kubelka–Munk absorption/scattering parameter KS for a given spectral reflectance R.
---
---In Kubelka–Munk theory, the KS function reflects the ratio that controls the conversion from spectral
---reflectance to an equivalent absorption/scattering coefficient. The formulation
---`(1 - R)² / (2 * R)` is a common approximation that assumes a diffusely scattering medium.
---@param R number The spectral reflectance value.
---@return number KS The computed KS value.
function M.KS(R)
    return (1 - R)^2 / (2*R)
end

---Computes the Kubelka–Munk mixing coefficient KM from a given KS value.
---
---The KM function transforms the KS parameter into a measure that can be linearly mixed.
---This conversion is essential because the Kubelka–Munk model assumes that when pigments are
---mixed, the resulting reflectance is a function of the weighted combination of the pigment
---absorption and scattering properties. The formula used here:
---
---KM(KS) = 1 + KS - √(KS² + 2KS)
---
---provides the appropriate transformation for blending multiple pigment spectra.
---@param KS number The KS value (absorption/scattering parameter).
---@return number KM The computed KM mixing coefficient.
function M.KM(KS)
    return 1 + KS - math.sqrt(KS^2 + 2*KS)
end


---Mixes multiple colors using a model based on the Kubelka–Munk theory.
---
---This function implements a mixing algorithm that is inspired by the Kubelka–Munk theory,
---which models how light interacts with diffusely scattering and absorbing layers (such as pigments or paints).
---The approach is as follows:
---
---- For each wavelength band (with R_SAMPLES samples), compute a weighted average of the KS values.
---- Weights are determined by the square of a factor that considers both the square-root of the color's
---  luminance and its tinting strength multiplied by a user-specified weight.
---- The resulting weighted KS average is then converted back using the KM function to obtain the
---  mixed spectral reflectance.
---
---In effect, this method blends pigments based on their optical absorption and scattering properties,
---providing a physically motivated approximation for pigment mixing as described by Kubelka–Munk.
---
---@param colors table vector of colors, each color having fields XYZ and R
---@param weights table vector of weights aka. factor for each color. 
---@return table R vector of reflectances.
function M._mix(colors, weights)
    local R_ret = {}
    for i_R = 1, C.R_SAMPLES do
        local ksMix = 0
        local totalConcentration = 0
        for i_color, color in ipairs(colors) do
            local R = color.R[i_R]
            local weight = weights[i_color]
            -- tintingStrength can be customised for a given colour but is just 1 by default.
            -- local concentration = weight ^ 2 * color.tintingStrength ^ 2 * color.luminance
            local luminance = math.max(C.EPSILON, color.XYZ.y)
            -- Luminance can also be modified in spectral.js but defaults to max(EPSILON, XYZ[1]) which is XYZ[2] in lua 1-indexing.
            local concentration = weight ^ 2 * luminance
            totalConcentration = totalConcentration + concentration
            ksMix = ksMix + M.KS(R) * concentration
        end
        R_ret[i_R] = M.KM(ksMix / totalConcentration)
    end
    return R_ret
end

-- Convenient chaining functions.

---RGB to CIE Lab using default D65 illuminant.
---https://kaizoudou.com/from-rgb-to-lab-color-space/
---http://brucelindbloom.com/index.html
---@param sRGB table {r=, g=, b=}
---@return table CIELab
function M.sRGB_to_CIELab(sRGB)
    local lRGB = M.sRGB_to_lRGB(sRGB)
    local XYZ = M.lRGB_to_XYZ(lRGB)
    return M.XYZ_to_CIELab(XYZ)
end

---CIE Lab to sRGB conversion.
---Inverse of sRGB_to_CIELab.
---@param CIELab table
---@return table sRGB
function M.CIELab_to_sRGB(CIELab)
    local XYZ = M.CIELab_to_XYZ(CIELab)
    local lRGB = M.XYZ_to_lRGB(XYZ)
    return M.lRGB_to_sRGB(lRGB)
end

function M.hex_to_HSL(hex)
    return M.sRGB_to_HSL(M.hex_to_sRGB(hex))
end

function M.HSL_to_hex(hsl)
    return M.sRGB_to_hex(M.HSL_to_sRGB(hsl))
end

function M.hex_to_CIELab(hex)
    return M.sRGB_to_CIELab(M.hex_to_sRGB(hex))
end

function M.CIELab_to_hex(cielab)
    return M.sRGB_to_hex(M.CIELab_to_sRGB(cielab))
end

---Converts linear RGB values directly to OKLab by first converting to XYZ.
---@param lRGB table Linear RGB values.
---@return table OKLab
function M.lRGB_to_OKLab(lRGB)
    return M.XYZ_to_OKLab(M.lRGB_to_XYZ(lRGB))
end

function M.hex_to_R(hex_str)
    local sRGB = M.hex_to_sRGB(hex_str)
    local lRGB = M.sRGB_to_lRGB(sRGB)
    return M.lRGB_to_R(lRGB)
end

function M.R_to_hex(R)
    local XYZ = M.R_to_XYZ(R)
    local lRGB = M.XYZ_to_lRGB(XYZ)
    local sRGB = M.lRGB_to_sRGB(lRGB)
    return M.sRGB_to_hex(sRGB)
end

---Mix colors given their reflectances.
---@param Rs table vector of reflectance vectors.
---@param weights? table weight aka. factor for each color. If omitted then equal weights are given.
---@return table R Vector of reflectances.
function M.mix(Rs, weights)
    local colors = {}
    weights = weights or {}
    for i, R in ipairs(Rs) do
        colors[i] = {R=R, XYZ=M.R_to_XYZ(R)}
        weights[i] = weights[i] or (1/#Rs)
    end
    return M._mix(colors, weights)
end


return M
