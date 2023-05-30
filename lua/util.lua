#!/usr/bin/env lua
M = {}

function M.end_visual()
    -- magic from
    -- https://github.com/neovim/neovim/issues/19770
    vim.api.nvim_feedkeys('\027', 'xt', false)
end

function M.get_visual_range()
    local r1, c1 = unpack(vim.api.nvim_buf_get_mark(0, "<"))
    local r2, c2 = unpack(vim.api.nvim_buf_get_mark(0, ">"))
    -- handle edge-case where final char is unicode or other multibyte char
    local char = vim.api.nvim_buf_get_text(0, r2-1, c2, r2-1, c2+1, {})[1]
    -- if multibyte then char is only half the symbol and won't match a broad pattern like:
    if not char:match('[%w%p%s]') then c2=c2+1 end
    return r1, c1, r2, c2
end

-- return whether r1, c1 is before r2, c2
function M.before(r1, c1, r2, c2)
    return r1 < r2 or (r1 == r2 and c1 < c2)
end

return M
