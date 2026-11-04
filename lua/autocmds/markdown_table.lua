local util = require "utils/init"
local md_table = require "utils.markdown_table"

--- Cap on column width: the count if given, else the default.
--- Read before leaving visual mode, which resets the count.
--- @return integer
local function max_width()
    return vim.v.count > 0 and vim.v.count or md_table.MAX_WIDTH
end

--- Align the markdown tables intersecting rows [r1, r2] of the current buffer.
--- @param r1 integer 0-indexed, inclusive
--- @param r2 integer 0-indexed, inclusive
--- @param width integer cap on column width
local function align_rows(r1, r2, width)
    if md_table.align(0, r1, r2, { max_width = width, align_cells = true }) == 0 then
        vim.notify("not on a markdown table", vim.log.levels.INFO)
    end
end

vim.api.nvim_create_autocmd("FileType", {
    group = vim.api.nvim_create_augroup("markdown_table", { clear = true }),
    pattern = md_table.FILETYPES,
    callback = function(ev)
        vim.keymap.set("n", "<localleader>a", function()
            local row = util.get_cursor()
            align_rows(row, row, max_width())
        end, { buf = ev.buf, desc = "Align markdown table" })
        vim.keymap.set("x", "<localleader>a", function()
            local width = max_width()
            local r1, _, r2 = util.get_visual_range()
            align_rows(r1 - 1, r2 - 1, width)
        end, { buf = ev.buf, desc = "Align markdown tables" })
    end,
})
