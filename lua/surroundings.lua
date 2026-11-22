-- mini.surround surrounding specifications, by surrounding identifier.
local MiniSurround = require "mini.surround"
local util = require "utils/init"

---A mini.surround region (`:h MiniSurround-glossary`): 1-indexed lines and
---byte columns, `to` inclusive. A column past the end of a line is its newline.
---@alias SurroundRegion { from: { line: integer, col: integer }, to?: { line: integer, col: integer } }
---@alias SurroundRegionPair { left: SurroundRegion, right: SurroundRegion }

---A region pair of the current buffer widened to whole lines. `left` becomes
---the line where `left` starts, with its newline. `right` becomes the line
---where `right` starts, with the newline before it.
---@param pair SurroundRegionPair
---@return SurroundRegionPair
local function whole_lines(pair)
    local first, last = pair.left.from.line, pair.right.from.line
    local function n_cols(line) return #util.get_line(line - 1) end
    return {
        left = { from = { line = first, col = 1 }, to = { line = first, col = n_cols(first) + 1 } },
        right = { from = { line = last - 1, col = n_cols(last - 1) + 1 }, to = { line = last, col = n_cols(last) } },
    }
end

local fence = MiniSurround.gen_spec.input.treesitter { outer = "@fence.outer", inner = "@fence.inner" }

return {
    -- [d]elete, e.g. only as cstd to delete surround """
    ['d'] = { output = { left = '', right = '' } },
    -- normally used for tags such as <div> </div> which is cool but I don't use them so t for triple
    ['t'] = {
        -- figured out these patterns from the [[ ]] example in :h MiniSurround.config
        input = { { '%"%"%"().-()%"%"%"', "'''().-()'''" } },
        output = { left = '"""', right = '"""' }
    },
    -- Taken from `:h MiniSurround.config` for use in lua.
    -- I never really use the open versions of brackets that adds space so we override '['.
    ['['] = {
        input = { '%[%[().-()%]%]' },
        output = { left = '[[', right = ']]' },
    },
    -- for latex ``...'', ~ is chosen since it is shift+`
    ['~'] = {
        input = { "``().-()''" },
        output = { left = '``', right = "''" },
    },
    F = {
        input = function() return vim.tbl_map(whole_lines, fence()) end,
        output = { left = '```\n', right = '\n```' },
    },
}
