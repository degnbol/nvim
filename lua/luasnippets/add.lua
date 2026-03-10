-- Quick snippet creation from visual selection.
-- Appends a new snippet to luasnippets/<ft>.lua, then opens it for editing.

local util = require("utils/init")

local M = {}

local rtp = vim.opt.runtimepath:get()[1]

---@param ft string
---@return string
local function snippets_file(ft)
    return rtp .. "/luasnippets/" .. ft .. ".lua"
end

---@param str string
---@return string
local function escape_t(str)
    return str:gsub('\\', '\\\\'):gsub('"', '\\"')
end

---Format lines as a multi-line LuaSnip snippet body.
---@param lines string[]
---@return string[]
local function format_body(lines)
    if #lines == 1 then
        return { '{t"' .. escape_t(lines[1]) .. '"}' }
    end
    local out = { "{t{" }
    for _, line in ipairs(lines) do
        out[#out + 1] = '\t"' .. escape_t(line) .. '",'
    end
    out[#out + 1] = "}}"
    return out
end

---Create a snippet entry as lines.
---@param trigger string
---@param desc string
---@param lines string[]
---@return string[] entry_lines
---@return number body_start line offset (0-based) where body begins within entry
local function format_snippet(trigger, desc, lines)
    local body = format_body(lines)
    local entry = {}
    local desc_part = desc ~= "" and ', dscr="' .. escape_t(desc) .. '"' or ""
    entry[#entry + 1] = 's({trig="' .. trigger .. '"' .. desc_part .. "},"
    local body_start = #entry
    for _, l in ipairs(body) do
        entry[#entry + 1] = l
    end
    entry[#entry + 1] = "),"
    return entry, body_start
end

---Ensure the snippets file exists with the standard header.
---@param path string
local function ensure_file(path)
    if vim.fn.filereadable(path) == 1 then return end
    local f = io.open(path, "w")
    if not f then return end
    f:write("---@diagnostic disable: undefined-global\nreturn {\n}\n")
    f:close()
end

---Append snippet lines to the file before the final `}`.
---Returns the line number (1-based) where the body starts in the file.
---@param path string
---@param entry string[]
---@param body_offset number
---@return number? body_line
local function append_snippet(path, entry, body_offset)
    local file_lines = {}
    for line in io.lines(path) do
        file_lines[#file_lines + 1] = line
    end

    -- Find last `}` line (the return table's closing brace)
    local insert_at
    for idx = #file_lines, 1, -1 do
        if file_lines[idx]:match("^}%s*$") then
            insert_at = idx
            break
        end
    end
    if not insert_at then
        vim.notify("Could not find closing brace in " .. path, vim.log.levels.ERROR)
        return nil
    end

    -- Insert blank line + entry before the closing brace
    local new = {}
    for idx = 1, insert_at - 1 do
        new[#new + 1] = file_lines[idx]
    end
    new[#new + 1] = ""
    local entry_start = #new + 1
    for _, l in ipairs(entry) do
        new[#new + 1] = l
    end
    for idx = insert_at, #file_lines do
        new[#new + 1] = file_lines[idx]
    end

    local f = io.open(path, "w")
    if not f then return nil end
    f:write(table.concat(new, "\n") .. "\n")
    f:close()

    return entry_start + body_offset
end

---Add a snippet from visual selection. Prompts for trigger and description,
---then opens the file at the new snippet for further editing (placeholders, etc.).
function M.add_from_visual()
    local r1, c1, r2, c2 = util.get_visual_range()
    local lines = vim.api.nvim_buf_get_text(0, r1 - 1, c1, r2 - 1, c2, {})
    if #lines == 0 then
        vim.notify("No selection", vim.log.levels.WARN)
        return
    end

    local ft = vim.bo.filetype:match("^[^.]+")

    vim.ui.input({ prompt = "Trigger: " }, function(trigger)
        if not trigger or trigger == "" then return end
        vim.ui.input({ prompt = "Description: " }, function(desc)
            desc = desc or ""
            local path = snippets_file(ft)
            ensure_file(path)
            local entry, body_offset = format_snippet(trigger, desc, lines)
            local body_line = append_snippet(path, entry, body_offset)
            -- Open file at the snippet body for editing
            vim.cmd.edit(path)
            if body_line then
                vim.api.nvim_win_set_cursor(0, { body_line, 0 })
            end
        end)
    end)
end

---Open the luasnippets file for the current filetype.
function M.edit()
    local ft = vim.bo.filetype:match("^[^.]+")
    local path = snippets_file(ft)
    if vim.fn.filereadable(path) == 0 then
        ensure_file(path)
    end
    vim.cmd.edit(path)
end

return M
