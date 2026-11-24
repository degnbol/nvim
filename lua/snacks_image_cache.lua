--- Key snacks.image's conversion cache on the source file's mtime and size, not only its path.
--- snacks names each converted file after the source path and reuses any that
--- exists (snacks/image/convert.lua `Convert.new`, `Convert:_resolve`), so an
--- image changed on disk keeps showing its old conversion, also after a restart
--- (folke/snacks.nvim#2949, closed without a fix). Remove once fixed upstream.
--- Limits: the key covers the source file only, not files it includes
--- (`\input`, `#include`, linked images).
local util = require "utils/init"

local M = {}

--- Stamp of a file's version.
--- @param path string
--- @return string|nil stamp mtime in ns, `-`, size in bytes; nil when the file does not exist
local function stamp(path)
    local stat = vim.uv.fs_stat(path)
    return stat and ("%d%09d-%d"):format(stat.mtime.sec, stat.mtime.nsec, stat.size)
end

--- Remove the converted files of a source, stamped or not, except those of one version.
--- @param dir string cache directory
--- @param unstamped string file name prefix snacks gives the source
--- @param kept string file name prefix of the version to keep
local function remove_stale(dir, unstamped, kept)
    for name in vim.fs.dir(dir) do
        if vim.startswith(name:gsub("^%d+%-%d+%-", ""), unstamped .. ".")
            and not vim.startswith(name, kept .. ".") then
            vim.fs.rm(dir .. "/" .. name)
        end
    end
end

--- Add the source file's version stamp to the name of every converted file.
--- A changed source then gets new converted files, which snacks displays as a
--- new image. After a successful conversion of the source's current version,
--- the files of other versions are removed, as snacks never removes cache files.
--- The stamp goes first: snacks reads the magick args from the `.<format>.png`
--- end of the name (`commands.convert` in convert.lua).
--- Patches the private `Convert` class of `snacks.image.convert` in place. Call once.
function M.setup()
    local Convert = util.upvalue(require("snacks.image.convert").convert, "Convert")

    local resolve = Convert.resolve
    function Convert:resolve()
        self._stamp = stamp(self.src)
        if self._stamp then
            self._unstamped = self.prefix
            self.prefix = self._stamp .. "-" .. self.prefix
        end
        return resolve(self)
    end

    local on_done = Convert.on_done
    function Convert:on_done()
        on_done(self)
        if self._stamp and self:ready() and stamp(self.src) == self._stamp then
            remove_stale(require("snacks").image.config.cache, self._unstamped, self.prefix)
        end
    end
end

return M
