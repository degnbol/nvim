#!/usr/bin/env lua
-- this script sets vartabstop on file save based on longest cell in each column.

-- Run automatically sometimes. Maybe also run after pressing tab. 
vim.api.nvim_create_autocmd({"BufEnter", "BufWritePost"},
{buffer=0, callback=function()
    local vartabstop = {}
    
    -- check at most 100 lines.
    local lines = vim.api.nvim_buf_get_lines(0, 0, 100, false)
    for _, line in ipairs(lines) do
        local fields = vim.split(line, '\t', true)
        -- when the line is shorter than or equal in length to a line seen so far
        for i = 1, math.min(#vartabstop, #fields) do
            vartabstop[i] = math.max(vartabstop[i], fields[i]:len())
        end
        -- when the line is longer than any line seen so far
        for i = #vartabstop+1, #fields do
            table.insert(vartabstop, fields[i]:len())
        end
    end

    -- tab represented a minimum 2 spaces.
    for i, v in ipairs(vartabstop) do
        vartabstop[i] = v+2
    end
    
    -- apply
    vim.opt_local.vartabstop = vartabstop
end})

