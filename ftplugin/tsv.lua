#!/usr/bin/env lua
api = vim.api
opt_local = vim.opt_local
require 'split'

function Vartabstop()
    local vartabstop = {}
    
    -- first take values from header, so we don't have to check for nil for every 
    -- line.
    local line = api.nvim_buf_get_lines(0, 0, 1, true)[1]
    for field in line:split('\t') do
        table.insert(vartabstop, field:len())
    end
    
    -- iterate each line
    local lines = api.nvim_buf_get_lines(0, 1, -1, true)
    for _, line in ipairs(lines) do
        local i = 1
        for field in line:split('\t') do
            vartabstop[i] = math.max(vartabstop[i], field:len())
            i=i+1
        end
    end
    
    -- tab represented as minimum 2 spaces.
    for i, v in ipairs(vartabstop) do
        vartabstop[i] = v+2
    end
    
    opt_local.vartabstop = vartabstop
end

-- run automatically sometimes. Maybe add after making changes, e.g. leaving 
-- insert mode or TextChanged(I) event(s).
api.nvim_create_autocmd(
    {"BufEnter", "BufWritePost"},
    {pattern="*.tsv", callback=Vartabstop}
)

