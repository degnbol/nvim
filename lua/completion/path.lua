local lib = require "blink.cmp.sources.path.lib"
local paths = require "utils.paths"

local M = {}

---@type blink.cmp.PathOpts
---@diagnostic disable-next-line: missing-fields
M.opts = {
    get_cwd = function(context) return paths.buf_dir(context.bufnr) end,
}

--- The directory whose entries complete the keyword, as blink's path source
--- resolves it, after the `$(git root)` and $NAME forms in the word before the
--- keyword are expanded against the buffer.
---@param opts blink.cmp.PathOpts read for `get_cwd` and `ignore_root_slash`
---@param context blink.cmp.Context
---@return string|nil dir nil when the directory does not exist
function M.expanded_dirname(opts, context)
    local start_col = context.bounds.start_col
    local before = context.line:sub(1, start_col - 1)
    -- `$(git root)/…` contains a space, so try it before a plain word.
    local i_word = before:match("()%$%b()%S*$") or before:match("()%S*$")
    local word = before:sub(i_word)
    -- regex.PATH cannot get past `//`, which a value ending in `/` leaves.
    local expanded = (paths.expand_dollars(word, context.bufnr):gsub("/+", "/"))
    -- A fresh bounds: lib.candidates and other providers read the real context.
    local view = setmetatable({
        line = before:sub(1, i_word - 1) .. expanded .. context.line:sub(start_col),
        bounds = { start_col = start_col + #expanded - #word, length = context.bounds.length },
    }, { __index = context })
    local dir = lib.dirname(opts, view)
    if dir and paths.is_dir(dir) then return dir end
    return nil
end

return M
