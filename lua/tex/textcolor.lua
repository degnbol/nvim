local hi = require "utils/highlights"

local M = {}

local ns = vim.api.nvim_create_namespace("textcolor")

local colors = { "gray", "red", "green", "blue", "magenta", "cyan", "yellow", "purple", "pink", "orange" }
---@type table<string, string> xcolor name -> fg colour name
local fgs = {}
for _, col in ipairs(colors) do
    fgs[col] = col
    -- add yellowtext etc. which are (usually) darker variants made for being able to read the text on a white paper background.
    fgs[col .. "text"] = col
end

---Highlight group for an xcolor name.
---@param name string xcolor name
---@return string group
local function group_of(name)
    return "texTextcolor" .. name:sub(1, 1):upper() .. name:sub(2)
end

---Bodies of the `\textcolor{<name>}{<body>}` commands on a line, for the base
---xcolor names and their `<name>text` variants. A command split over lines is
---not found. Empty bodies are skipped.
---@param line string
---@return { col: integer, end_col: integer, group: string }[] spans 0-indexed byte columns, end exclusive, braces excluded. `group` is the highlight group for the colour.
function M.spans(line)
    local spans = {}
    for name, c1, body in line:gmatch('\\textcolor{(%a+)}()(%b{})') do
        if fgs[name] and #body > 2 then
            -- c1 is the 1-indexed `{`, so the 0-indexed start of the body.
            spans[#spans + 1] = { col = c1, end_col = c1 + #body - 2, group = group_of(name) }
        end
    end
    return spans
end

local function on_range(_, _, buf, brow, _, erow, ecol)
    -- End (lnum, 0) means EOL of the preceding line is the last included point.
    local last = ecol == 0 and erow - 1 or erow
    if last < brow then return end
    for i, line in ipairs(vim.api.nvim_buf_get_lines(buf, brow, last + 1, false)) do
        for _, span in ipairs(M.spans(line)) do
            vim.api.nvim_buf_set_extmark(buf, ns, brow + i - 1, span.col, {
                end_col = span.end_col,
                hl_group = span.group,
                ephemeral = true,
            })
        end
    end
end

---Colour `\textcolor` bodies in tex windows, and define their groups on every
---ColorScheme. Call once.
function M.setup()
    hi.onColorScheme(function()
        for name, fg in pairs(fgs) do
            hi.set(group_of(name), { fg = fg })
        end
    end)
    vim.api.nvim_set_decoration_provider(ns, {
        on_win = function(_, _, buf) return vim.bo[buf].filetype == "tex" end,
        on_range = on_range,
    })
end

return M
