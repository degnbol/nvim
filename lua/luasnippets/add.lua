-- Quick snippet creation from visual selection.
-- Appends a new snippet to luasnippets/<ft>.lua, then opens it for editing.

local util = require("utils/init")
local ls = require("luasnip")

local M = {}

local rtp = vim.opt.runtimepath:get()[1]

---@param ft string
---@return string
local function snippets_file(ft)
    return rtp .. "/luasnippets/" .. ft .. ".lua"
end

---Escape for Lua double-quoted string context.
---@param str string
---@return string
local function escape_lua(str)
    return str:gsub('\\', '\\\\'):gsub('"', '\\"')
end

---Build a LuaSnip snippet that expands into a new snippet definition.
---Tab stops: 1=trigger name, 2=description, 3=choice between t"" and fmta() body.
---@param lines string[]
---@return table snippet
local function build_snippet(lines)
    -- Build t"" option lines
    local t_opt
    if #lines == 1 then
        t_opt = { '{t"' .. escape_lua(lines[1]) .. '"}),'}
    else
        t_opt = { "{t{" }
        for _, line in ipairs(lines) do
            t_opt[#t_opt + 1] = '\t"' .. escape_lua(line) .. '",'
        end
        t_opt[#t_opt + 1] = "}}),"
    end

    -- Build fmta option lines
    local fmta_opt
    if #lines == 1 then
        fmta_opt = { 'fmta("' .. escape_lua(lines[1]) .. '", {})),'}
    else
        fmta_opt = { "fmta([[" }
        for _, line in ipairs(lines) do
            fmta_opt[#fmta_opt + 1] = line
        end
        fmta_opt[#fmta_opt + 1] = "]], {})),"
    end

    return ls.s("", {
        ls.t('s({trig="'),
        ls.i(1, "name"),
        ls.t('", dscr="'),
        ls.i(2, "desc"),
        ls.t({'"},', ''}),
        ls.c(3, {
            ls.t(t_opt),
            ls.t(fmta_opt),
        }),
        ls.t({'', ''}),
    })
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

---Insert a blank line before the final `}` and return its line number (1-based).
---@param path string
---@return number? blank_line
local function insert_blank_line(path)
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

    local new = {}
    for idx = 1, insert_at - 1 do
        new[#new + 1] = file_lines[idx]
    end
    new[#new + 1] = ""
    local blank_line = #new
    for idx = insert_at, #file_lines do
        new[#new + 1] = file_lines[idx]
    end

    local f = io.open(path, "w")
    if not f then return nil end
    f:write(table.concat(new, "\n") .. "\n")
    f:close()

    return blank_line
end

---Add a snippet from visual selection. Opens the snippets file and expands the
---full entry as a snippet — tab through trigger name, description, and body format.
function M.add_from_visual()
    local r1, c1, r2, c2 = util.get_visual_range()
    local lines = vim.api.nvim_buf_get_text(0, r1 - 1, c1, r2 - 1, c2, {})
    if #lines == 0 then
        vim.notify("No selection", vim.log.levels.WARN)
        return
    end

    local ft = vim.bo.filetype:match("^[^.]+")
    local path = snippets_file(ft)
    ensure_file(path)

    local snippet = build_snippet(lines)
    local blank_line = insert_blank_line(path)

    vim.cmd.edit(path)
    if blank_line then
        -- Force blink.cmp to apply buffer-local keymaps (snippet_forward/backward)
        -- to this new buffer. Blink only applies them on InsertEnter.
        vim.api.nvim_exec_autocmds('InsertEnter', { buffer = 0 })
        ls.snip_expand(snippet, { pos = { blank_line - 1, 0 } })
    end
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
