#!/usr/bin/env lua

local ns = vim.api.nvim_create_namespace("textcolor")
vim.api.nvim_set_hl_ns(ns)

local colors = {"gray", "red", "green", "blue", "magenta", "cyan", "yellow", "purple", "pink", "orange"}
local fgs = {}
for _, col in ipairs(colors) do
    fgs[col] = col
    -- add yellowtext etc. which are (usually) darker variants made for being able to read the text on a white paper background.
    fgs[col .. "text"] = col
end
-- make hl groups
for col, fg in pairs(fgs) do
    vim.api.nvim_set_hl(ns, col, {fg=fg})
end

local function textcolorRange(buf, linestart, lineend)
    local lines = vim.api.nvim_buf_get_lines(buf, linestart, lineend, false) -- true errors on new empty file
    if table.concat(lines):match('\\textcolor') then
        -- only works if the whole command is on the same line when first completed
        for iLine, line in ipairs(lines) do
            local ln = linestart+iLine-1
            -- may have multiple matches on a single line
            for col, c1, content in line:gmatch('\\textcolor{(%a+)}()(%b{})') do
                if fgs[col] ~= nil then
                    -- -1=zero-index
                    vim.api.nvim_buf_add_highlight(buf, ns, col, ln, c1, c1+#content-2)
                end
            end
        end
    end
end

vim.api.nvim_buf_attach(0, false, {
    on_lines=function (_, buf, changedtick, linenum0, linenum1, linenumLast, preCount, del_codepoints, del_codeunits)
        vim.api.nvim_buf_clear_namespace(buf, ns, linenum0, math.max(linenum1, linenumLast))
        textcolorRange(buf, linenum0, math.max(linenum1, linenumLast)+1)
    end,
    on_reload=function (_, buf)
        textcolorRange(buf, 0, -1)
    end
})

-- color newly opened file
textcolorRange(0, 0, -1)

