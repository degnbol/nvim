---@diagnostic disable: undefined-global
-- Tests for utils.file: a bounded window of lines out of a file, gzip
-- included. The cap is the point — reaching a deep line must never mean
-- reading a whole file to get there.

local file = require "utils.file"
local util = require "utils/init"

describe("read_lines", function()
    local dir
    local max_bytes = file.MAX_BYTES

    before_each(function()
        dir = vim.fn.tempname()
        vim.fn.mkdir(dir, "p")
    end)

    after_each(function()
        file.MAX_BYTES = max_bytes
        vim.fs.rm(dir, { recursive = true, force = true })
    end)

    --- Write `lines` to a fixture named `name`, gzipping it when `name` ends in
    --- .gz. `flags` go to writefile — "b" for no trailing newline.
    local function fixture(name, lines, flags)
        local path = dir .. "/" .. name
        local plain = (path:gsub("%.gz$", ""))
        assert.are.equal(0, vim.fn.writefile(lines, plain, flags or ""))
        if plain ~= path then
            assert.are.equal(0, vim.system({ "gzip", "-f", plain }):wait().code)
        end
        return path
    end

    local function numbered(n)
        local lines = {}
        for i = 1, n do
            lines[i] = "line " .. i
        end
        return lines
    end

    it("reads a window out of the middle", function()
        local path = fixture("a.txt", numbered(300))
        assert.are.same({ "line 170", "line 171", "line 172" }, file.read_lines(path, 170, 3))
    end)

    it("stops at the end of the file", function()
        local path = fixture("a.txt", numbered(300))
        assert.are.same({ "line 299", "line 300" }, file.read_lines(path, 299, 10))
    end)

    it("returns nothing past the end of the file", function()
        local path = fixture("a.txt", numbered(300))
        assert.are.same({}, file.read_lines(path, 301, 10))
    end)

    it("returns nothing for an empty file", function()
        local path = fixture("empty.txt", {})
        assert.are.same({}, file.read_lines(path, 1, 10))
    end)

    it("keeps a last line that has no trailing newline", function()
        local path = fixture("bare.txt", numbered(3), "b")
        assert.are.same({ "line 1", "line 2", "line 3" }, file.read_lines(path, 1, 10))
    end)

    it("reads a window out of a gzipped file", function()
        local path = fixture("a.tsv.gz", numbered(300))
        assert.are.same({ "line 170", "line 171", "line 172" }, file.read_lines(path, 170, 3))
    end)

    it("reports a corrupt archive", function()
        local path = dir .. "/bad.gz" -- not via fixture(), which would gzip it
        assert.are.equal(0, vim.fn.writefile({ "not gzip at all" }, path))
        assert.has_error(function() file.read_lines(path, 1, 10) end)
    end)

    -- 20000 numbered lines are ~215 kB, so a cap of 128 kB cuts the file in
    -- half. It has to stay above the chunk size to bind at all: the cap is
    -- weighed once per chunk read.
    local CAP = 128 * 2 ^ 10

    for _, ext in ipairs({ "txt", "txt.gz" }) do
        it("reads a shallow line out of a capped " .. ext, function()
            local path = fixture("big." .. ext, numbered(20000))
            file.MAX_BYTES = CAP
            assert.are.same({ "line 3" }, file.read_lines(path, 3, 1))
        end)

        it("comes up empty when the cap stops short of the window in a " .. ext, function()
            local path = fixture("big." .. ext, numbered(20000))
            file.MAX_BYTES = CAP
            assert.are.same({}, file.read_lines(path, 19000, 5))
        end)

        it("returns whole lines up to where the cap stopped a " .. ext, function()
            local path = fixture("big." .. ext, numbered(20000))
            file.MAX_BYTES = CAP
            local lines = file.read_lines(path, 1, 20000)
            assert.is_true(#lines > 0 and #lines < 20000)
            assert.are.equal("line " .. #lines, lines[#lines])
        end)
    end

    it("never returns a line the cap cut in half", function()
        local path = fixture("wide.txt", { string.rep("x", 200000), "line 2" })
        file.MAX_BYTES = CAP
        assert.are.same({}, file.read_lines(path, 1, 10))
    end)

    it("stops inflating as soon as the window is full", function()
        -- A gzip stream is a concatenation of members, so one small member
        -- repeated stands in for a huge archive: ~1 GB inflated, MBs on disk.
        -- Reading it whole takes seconds; reading the first two lines must not.
        local member = assert(util.readtext(fixture("member.txt.gz", vim.fn["repeat"]({ "x" }, 100000))))
        local path = dir .. "/huge.txt.gz"
        local out = assert(io.open(path, "wb"))
        out:write(string.rep(member, 5000))
        out:close()

        local start = vim.uv.hrtime()
        assert.are.same({ "x", "x" }, file.read_lines(path, 1, 2))
        local window_ns = vim.uv.hrtime() - start

        -- Measured rather than a fixed bound: how long the whole archive takes
        -- is the machine's business, being a fraction of it is the code's.
        start = vim.uv.hrtime()
        assert.are.equal(0, vim.system({ "gzip", "-cd", "--", path }, { stdout = false }):wait().code)
        local whole_ns = vim.uv.hrtime() - start
        assert.is_true(window_ns * 2 < whole_ns,
            ("%.0f ms for two lines of a %.0f ms archive"):format(window_ns / 1e6, whole_ns / 1e6))
    end)
end)
