---@diagnostic disable: undefined-global
local kitty = require "utils/kitty"

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
