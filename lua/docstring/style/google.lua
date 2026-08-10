-- Google-style docstring parse/render/locate. Operates on *dedented* docstring
-- text (common body indentation already stripped): the language layer removes it
-- on extraction and re-applies it on write, keeping this module indentation- and
-- buffer-agnostic. Round-trips a DocModel: `parse(render(m)) == m`.
--
-- Entry format is `name: description` or `name (type): description`; this
-- config's convention is typeless (types live in the annotated signature), so a
-- param with no `type` renders without the parenthesised group. Multiline
-- descriptions continue on deeper-indented lines.

local M = {}

--- Section header for orphan entries — documented params absent from the
--- signature. Rendered as its own section (rather than an inline marker) so it
--- round-trips by header recognition and stays maximally visible.
M.ORPHAN_HEADER = "Orphaned Args"

--- Header text (before the colon) → the model section it feeds. Aliases map to
--- the same canonical target; render always emits the canonical header.
local HEADER_KIND = {
    Args = "params",
    Arguments = "params",
    Parameters = "params",
    Returns = "returns",
    Return = "returns",
    Raises = "raises",
    Raise = "raises",
    [M.ORPHAN_HEADER] = "unmatched",
}

local BODY_INDENT = "    "
local CONT_INDENT = "        "

--- The header a line declares (`Args:` → `Args`), or nil if it is not a bare
--- `Header:` line at column zero.
--- @param line string
--- @return string|nil
local function header_of(line)
    return line:match("^([%w][%w ]*):%s*$")
end

--- Match one entry line into its parts, or nil if the line is not an entry (a
--- name followed by `:` or `(type):`). The single source of truth for "is this
--- an entry line", shared by `entry_of` (parse) and `locate_doc_entries` (byte
--- ranges) so the two never diverge. The name allows leading stars (`*args`,
--- `**kwargs`) and dots (qualified exception classes in `Raises:`).
--- @param line string
--- @return string|nil name
--- @return integer|nil scol  0-indexed byte column where the name starts
--- @return string|nil type  nil when the parenthesised group is absent
--- @return string|nil description
local function match_entry(line)
    local indent, name, rest = line:match("^(%s*)([%*%w_%.]+)%s*(.-)$")
    if not name then return end
    local type, desc = rest:match("^%((.-)%)%s*:%s?(.*)$")
    if not type then
        desc = rest:match("^:%s?(.*)$")
        if not desc then return end
    end
    return name, #indent, type, desc
end

--- The `{ name, type, description }` an entry line declares, or nil if the line
--- is not an entry.
--- @param line string
--- @return docstring.Param|nil
local function entry_of(line)
    local name, _, type, desc = match_entry(line)
    if not name then return end
    --- @type docstring.Param
    local param = { name = name, type = type, description = assert(desc) }
    return param
end

--- Indentation width (in characters) of a line, or nil for a blank line.
--- @param line string
--- @return integer|nil
local function indent_of(line)
    if line:match("^%s*$") then return end
    return #line:match("^%s*")
end

--- Whether the `Header:` line at `lines[i]` opens a section: a following body
--- line (after any blanks) must be indented. A bare `Word:` with no indented
--- body is prose, not a section — so a summary ending in `Note:` is not promoted
--- to a section, and a section render always emits an indented body (see
--- `render_section`), keeping the two consistent.
--- @param lines string[]
--- @param i integer
--- @return boolean
local function heads_section(lines, i)
    for j = i + 1, #lines do
        local indent = indent_of(lines[j])
        if indent then return indent > 0 end
    end
    return false
end

--- Split lines into `{ summary, sections }`, where each section is
--- `{ header, lines }` (body lines, trailing blanks trimmed). Summary is every
--- line before the first header, trailing blanks trimmed.
--- @param lines string[]
--- @return string[] summary
--- @return { header: string, lines: string[] }[] sections
local function split_sections(lines)
    local summary = {}
    local sections = {}
    --- @type { header: string, lines: string[] }|nil
    local current = nil
    for i, line in ipairs(lines) do
        local header = header_of(line)
        if header and not heads_section(lines, i) then header = nil end
        if header then
            current = { header = header, lines = {} }
            sections[#sections + 1] = current
        elseif current then
            current.lines[#current.lines + 1] = line
        else
            summary[#summary + 1] = line
        end
    end
    local function trim_trailing(list)
        while #list > 0 and list[#list]:match("^%s*$") do
            list[#list] = nil
        end
    end
    trim_trailing(summary)
    for _, section in ipairs(sections) do
        trim_trailing(section.lines)
    end
    return summary, sections
end

--- Parse a param/raise/orphan section body into entries. A new entry begins at
--- the base-indent entry line; deeper or blank lines continue the current
--- entry's description (each continuation left-stripped, joined by newlines).
--- @param body string[]
--- @return docstring.Param[] entries
local function parse_entries(body)
    local base = nil
    for _, line in ipairs(body) do
        local indent = indent_of(line)
        if indent then
            base = indent
            break
        end
    end
    local entries = {}
    --- @type docstring.Param|nil
    local current = nil
    for _, line in ipairs(body) do
        local indent = indent_of(line)
        local entry = indent == base and entry_of(line) or nil
        if entry then
            entries[#entries + 1] = entry
            current = entry
        elseif current then
            local cont = line:match("^%s*(.-)%s*$")
            current.description = current.description .. "\n" .. cont
        end
    end
    for _, entry in ipairs(entries) do
        entry.description = entry.description:gsub("%s+$", "")
    end
    return entries
end

--- The type prefix of a `type: desc` first line, or nil if the line is a plain
--- description. A type is a dotted name with an optional balanced subscript
--- (`bool`, `typing.Optional`, `dict[str, int]`), never multi-word prose — so a
--- description that merely contains a colon (`The parsed config: see below`) is
--- not mistaken for a type. Residual ambiguity remains for a single word before
--- the colon (`Note: it varies.` reads as `type=Note`): indistinguishable from a
--- typed return without out-of-band knowledge of the docstring's convention.
--- @param line string
--- @return string|nil type
--- @return string|nil description
local function split_return_type(line)
    local type, desc = line:match("^([%w_%.]+%b[])%s*:%s?(.*)$")
    if not type then
        type, desc = line:match("^([%w_%.]+)%s*:%s?(.*)$")
    end
    return type, desc
end

--- Parse a Returns body into a single nameless slot. `type` comes from a leading
--- `type:` prefix on the first line, if present.
--- @param body string[]
--- @return docstring.Return|nil
local function parse_return(body)
    if #body == 0 then return end
    local lines = {}
    for _, line in ipairs(body) do
        lines[#lines + 1] = line:match("^%s*(.-)%s*$")
    end
    local first = table.remove(lines, 1)
    local type, desc = split_return_type(first)
    if not type then
        desc = first
    end
    table.insert(lines, 1, desc)
    --- @type docstring.Return
    local ret = { type = type, description = table.concat(lines, "\n"):gsub("%s+$", "") }
    return ret
end

--- The `other` section body dedented to its minimum indent, so render can
--- re-indent uniformly.
--- @param body string[]
--- @return string
local function dedent_body(body)
    local min = nil
    for _, line in ipairs(body) do
        local indent = indent_of(line)
        if indent and (not min or indent < min) then min = indent end
    end
    min = min or 0
    local lines = {}
    for _, line in ipairs(body) do
        lines[#lines + 1] = line:match("^%s*$") and "" or line:sub(min + 1)
    end
    return table.concat(lines, "\n")
end

--- Parse dedented Google-style docstring text into a DocModel.
--- @param text string
--- @return docstring.Model model
function M.parse(text)
    local lines = vim.split(text, "\n", { plain = true })
    local summary, sections = split_sections(lines)

    local model = require("docstring.model").empty()
    model.summary = table.concat(summary, "\n")
    for _, section in ipairs(sections) do
        local kind = HEADER_KIND[section.header]
        if kind == "params" then
            for _, entry in ipairs(parse_entries(section.lines)) do
                model.params[#model.params + 1] = entry
            end
        elseif kind == "unmatched" then
            for _, entry in ipairs(parse_entries(section.lines)) do
                model.unmatched[#model.unmatched + 1] = entry
            end
        elseif kind == "raises" then
            for _, entry in ipairs(parse_entries(section.lines)) do
                --- @type docstring.Raise
                local raise = { type = entry.name, description = entry.description }
                model.raises[#model.raises + 1] = raise
            end
        elseif kind == "returns" then
            model.returns = parse_return(section.lines)
        else
            model.other[#model.other + 1] =
                { header = section.header, body = dedent_body(section.lines) }
        end
    end
    return model
end

--- Render one param entry: `name (type): first` with continuation lines at
--- CONT_INDENT.
--- @param param docstring.Param
--- @param out string[]
local function render_entry(param, out)
    local head = param.type and (param.name .. " (" .. param.type .. ")") or param.name
    local desc = vim.split(param.description, "\n", { plain = true })
    out[#out + 1] = BODY_INDENT .. (desc[1] ~= "" and (head .. ": " .. desc[1]) or (head .. ":"))
    for i = 2, #desc do
        out[#out + 1] = desc[i] == "" and "" or (CONT_INDENT .. desc[i])
    end
end

--- Append a section (blank separator, header, body) to `out` if `render_body`
--- produces any lines.
--- @param out string[]
--- @param header string
--- @param render_body fun(body: string[])
local function render_section(out, header, render_body)
    local body = {}
    render_body(body)
    if #body == 0 then return end
    if #out > 0 then out[#out + 1] = "" end
    out[#out + 1] = header .. ":"
    for _, line in ipairs(body) do
        out[#out + 1] = line
    end
end

--- Render a DocModel to dedented Google-style docstring text.
--- @param model docstring.Model
--- @return string text
function M.render(model)
    local out = {}
    if model.summary ~= "" then
        for _, line in ipairs(vim.split(model.summary, "\n", { plain = true })) do
            out[#out + 1] = line
        end
    end

    render_section(out, "Args", function(body)
        for _, param in ipairs(model.params) do
            render_entry(param, body)
        end
    end)
    render_section(out, "Returns", function(body)
        local ret = model.returns
        if not ret then return end
        local desc = vim.split(ret.description, "\n", { plain = true })
        local first = ret.type and ret.type ~= "" and (ret.type .. ":") or ""
        if desc[1] ~= "" then
            first = first == "" and desc[1] or (first .. " " .. desc[1])
        end
        if first == "" then return end -- no type and no description → no section
        body[#body + 1] = BODY_INDENT .. first
        for i = 2, #desc do
            body[#body + 1] = desc[i] == "" and "" or (CONT_INDENT .. desc[i])
        end
    end)
    render_section(out, "Raises", function(body)
        for _, raise in ipairs(model.raises) do
            render_entry({ name = raise.type, description = raise.description }, body)
        end
    end)
    for _, section in ipairs(model.other) do
        render_section(out, section.header, function(body)
            for _, line in ipairs(vim.split(section.body, "\n", { plain = true })) do
                body[#body + 1] = line == "" and "" or BODY_INDENT .. line
            end
        end)
    end
    render_section(out, M.ORPHAN_HEADER, function(body)
        for _, param in ipairs(model.unmatched) do
            render_entry(param, body)
        end
    end)

    return table.concat(out, "\n")
end

--- Locate each documented param's name token in dedented Google-style text.
--- Ranges are 0-indexed byte (row, col), relative to `text`; the caller offsets
--- them by the docstring's buffer position. Covers `Args:` and orphan entries,
--- not `Raises:` (exception classes are not signature params).
--- @param text string
--- @return { name: string, srow: integer, scol: integer, erow: integer, ecol: integer }[]
function M.locate_doc_entries(text)
    local lines = vim.split(text, "\n", { plain = true })
    local kind = nil
    local base = nil
    local located = {}
    for row, line in ipairs(lines) do
        local header = header_of(line)
        if header then
            kind = HEADER_KIND[header]
            base = nil
        elseif kind == "params" or kind == "unmatched" then
            local indent = indent_of(line)
            if indent then
                base = base or indent
                if indent == base then
                    local name, scol = match_entry(line)
                    if name then
                        located[#located + 1] = {
                            name = name,
                            srow = row - 1,
                            scol = scol,
                            erow = row - 1,
                            ecol = scol + #name,
                        }
                    end
                end
            end
        end
    end
    return located
end

return M
