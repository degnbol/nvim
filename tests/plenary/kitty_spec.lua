---@diagnostic disable: undefined-global
local kitty = require "utils/kitty"
local colors = require "utils/colors"

describe("kitty.text_alpha_fx", function()
    if vim.fn.executable("magick") == 0 then
        pending("magick not on PATH")
        return
    end

    ---Alphas, 0–1, of a 256-step 0-to-1 alpha ramp.
    ---@param fx string|nil `-fx` expression to apply to the alpha channel first
    ---@return number[] alphas
    local function ramp_alphas(fx)
        local args = { "magick", "-size", "1x256", "gradient:", "-alpha", "copy" }
        if fx then vim.list_extend(args, { "-channel", "A", "-fx", fx, "+channel" }) end
        vim.list_extend(args, { "-alpha", "extract", "-depth", "8", "gray:-" })
        local res = vim.system(args):wait()
        assert.are.equal(0, res.code, res.stderr)
        return vim.tbl_map(function(byte) return byte / 255 end, { res.stdout:byte(1, -1) })
    end

    local composition = { gamma = 1.7, contrast = 30 }
    for _, pair in ipairs { { "#e0e2ea", "#14161b" }, { "#575279", "#faf4ed" } } do
        it(("matches kitty's shader for %s on %s"):format(pair[1], pair[2]), function()
            local L_fg, L_bg = colors.hex_to_luminance(pair[1]), colors.hex_to_luminance(pair[2])
            local w = (1 - L_fg + L_bg) * 0.5
            local alphas = ramp_alphas()
            local adjusted = ramp_alphas(kitty.text_alpha_fx(composition, L_fg, L_bg))
            assert.are.equal(256, #adjusted)
            for i, a in ipairs(alphas) do
                local mixed = a + (a ^ (1 / composition.gamma) - a) * w
                local expected = math.min(math.max(mixed * (1 + composition.contrast / 100), 0), 1)
                assert.near(expected, adjusted[i], 1 / 255)
            end
        end)
    end
end)

describe("kitty text_composition script", function()
    if vim.fn.executable("kitty") == 0 then
        pending("kitty not on PATH")
        return
    end

    ---Runs the script against a kitty.conf with only a `text_composition_strategy` line.
    ---@param strategy string option value
    ---@return vim.SystemCompleted
    local function run(strategy)
        local dir = vim.fn.tempname()
        vim.fn.mkdir(dir)
        vim.fn.writefile({ "text_composition_strategy " .. strategy }, dir .. "/kitty.conf")
        local script = vim.fn.stdpath("config") .. "/scripts/kitty_text_composition.py"
        return vim.system({ "kitty", "+launch", script },
            { text = true, env = { KITTY_CONFIG_DIRECTORY = dir } }):wait()
    end

    ---The script's decoded output for `strategy`. Asserts that it exits 0.
    ---@param strategy string option value
    ---@return kitty.TextComposition
    local function parse(strategy)
        local res = run(strategy)
        assert.are.equal(0, res.code, res.stderr)
        return vim.json.decode(res.stdout)
    end

    if vim.fn.has("mac") == 1 then
        it("maps platform to 1.7 30 on macOS", function()
            assert.are.same({ gamma = 1.7, contrast = 30 }, parse("platform"))
        end)
    end

    it("reads gamma and contrast", function()
        assert.are.same({ gamma = 2, contrast = 50 }, parse("2 50"))
    end)

    it("clamps gamma and contrast", function()
        assert.are.same({ gamma = 0.01, contrast = 100 }, parse("0 200"))
    end)

    it("rejects legacy", function()
        assert.are_not.equal(0, run("legacy").code)
    end)
end)

describe("kitty.cell_fractions", function()
    if vim.fn.executable("kitty") == 0 then
        pending("kitty not on PATH")
        return
    end

    local fr = assert(kitty.cell_fractions())

    it("names the face", function()
        assert.is_string(fr.postscript_name)
    end)

    it("orders the vertical layout", function()
        assert.is_true(0 < fr.desc)
        assert.is_true(fr.desc < fr.x)
        assert.is_true(fr.x < fr.cap)
        assert.is_true(fr.cap < 1 - fr.below)
        assert.is_true(1 - fr.below < 1)
    end)

    it("rejects a fixed-size modify_font offset", function()
        local dir = vim.fn.tempname()
        vim.fn.mkdir(dir)
        vim.fn.writefile({ "modify_font baseline 3" }, dir .. "/kitty.conf")
        local script = vim.fn.stdpath("config") .. "/scripts/kitty_cell_fractions.py"
        local res = vim.system({ "kitty", "+launch", script },
            { text = true, env = { KITTY_CONFIG_DIRECTORY = dir } }):wait()
        assert.are_not.equal(0, res.code)
        assert.matches("baseline 3", res.stderr, 1, true)
    end)

    -- Snacks reads the cell size with an ioctl on fd 1 and silently falls back to 9×18 when fd 1 is not a tty.
    if not (kitty.term() and vim.uv.guess_handle(1) == "tty") then
        pending("matches the live cell aspect (stdout is not a kitty tty)")
        return
    end

    it("matches the live cell aspect", function()
        vim.cmd.packadd("snacks.nvim")
        local cell = require("snacks").image.terminal.size()
        -- Kitty rounds both sides to whole pixels, which moves height / width by up to (0.5 + 0.5·aspect) / width.
        -- The 10× measurement rounds too, so allow twice that.
        assert.near(cell.cell_height / cell.cell_width, fr.aspect, (1 + fr.aspect) / cell.cell_width)
    end)
end)
