---@diagnostic disable: undefined-global
vim.cmd.packadd "snacks.nvim"
local snacks = require "snacks"
local convert = require "snacks.image.convert"

describe("snacks_image_cache", function()
    local dir = vim.fn.tempname()
    vim.fn.mkdir(dir, "p")
    snacks.image.config.cache = dir .. "/cache"
    snacks.image.terminal.size = function() return { cell_width = 8, cell_height = 16, scale = 1 } end
    require("snacks_image_cache").setup()
    local src = dir .. "/a.png"

    --- @return string[] output file of each conversion step of `src`
    local function cache_files()
        return vim.tbl_map(function(step) return step.file end, convert.convert({ src = src }).steps)
    end

    --- Write `content` to `src` and set its access and modification times to `mtime`.
    --- @param content string
    --- @param mtime number seconds since the epoch
    local function write(content, mtime)
        vim.fn.writefile({ content }, src, "b")
        vim.uv.fs_utime(src, mtime, mtime)
    end

    it("keeps the cache files of an unchanged source", function()
        write("x", 1000)
        assert.are.same(cache_files(), cache_files())
    end)

    it("changes the cache files when the source mtime changes", function()
        write("x", 1000)
        local before = cache_files()
        write("x", 2000)
        assert.are_not.same(before, cache_files())
    end)

    it("changes the cache files when the source size changes", function()
        write("x", 1000)
        local before = cache_files()
        write("xy", 1000)
        assert.are_not.same(before, cache_files())
    end)

    it("keeps the .<format>.png tail that selects the magick args", function()
        local math_src = dir .. "/a.math.tex"
        vim.fn.writefile({ "x" }, math_src)
        local pngs = vim.tbl_filter(function(file) return file:match("%.png$") end,
            vim.tbl_map(function(step) return step.file end, convert.convert({ src = math_src }).steps))
        assert.is_truthy(vim.fs.basename(pngs[1]):match("%.math%.png$"))
    end)

    describe("after a conversion", function()
        local cache = snacks.image.config.cache
        local jpg = dir .. "/b.jpg"

        --- Conversion of `jpg` with the files of its current version already in the cache,
        --- so `run` finishes without running a command.
        --- @return snacks.image.Convert
        --- @return string[] names of the files of the current version
        --- @return string[] names of the files of other versions of `jpg`, and of another source
        local function converted()
            vim.fn.writefile({ "x" }, jpg)
            vim.uv.fs_utime(jpg, 1000, 1000)
            vim.fs.rm(cache, { recursive = true, force = true })
            vim.fn.mkdir(cache, "p")
            local c = convert.convert({ src = jpg })
            local current = vim.tbl_map(function(step) return vim.fs.basename(step.file) end, c.steps)
            -- The name snacks gives a converted file of `jpg`, without our stamp.
            local unstamped = current[1]:gsub("^%d+%-%d+%-", ""):match("^(.-)%.")
            local others = { "1-1-" .. unstamped .. ".png", unstamped .. ".png", "x" .. unstamped .. ".png" }
            for _, name in ipairs(vim.list_extend(vim.list_extend({}, current), others)) do
                vim.fn.writefile({}, cache .. "/" .. name)
            end
            return c, current, others
        end

        --- @param name string
        --- @return boolean
        local function cached(name) return vim.uv.fs_stat(cache .. "/" .. name) ~= nil end

        it("removes the files of other versions of the source only", function()
            local c, current, others = converted()
            c:run()
            assert.are.same({ true, true }, vim.tbl_map(cached, current))
            assert.are.same({ false, false, true }, vim.tbl_map(cached, others))
        end)

        it("removes nothing when the source changed since the conversion started", function()
            local c, current, others = converted()
            vim.uv.fs_utime(jpg, 2000, 2000)
            c:run()
            assert.are.same({ true, true }, vim.tbl_map(cached, current))
            assert.are.same({ true, true, true }, vim.tbl_map(cached, others))
        end)
    end)
end)
