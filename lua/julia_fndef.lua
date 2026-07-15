--- Toggle a Julia function definition between short (`f(x) = …`) and long
--- (`function f(x) … end`) form, and conform an overload set to one form.

local M = {}

local ts = require "utils/treesitter"
local gettext = vim.treesitter.get_node_text

--- @class JuliaDef
--- @field range integer[]  {start_row, start_col, end_row, end_col}, 0-indexed
--- @field form "short"|"long"
--- @field sig string   signature text (callee + args + any `::T` / `where`)
--- @field callee string  function name (`f` or qualified `Base.show`)
--- @field body string   short: RHS text; long: the block's sole expression
--- @field single boolean  whether the body is a single expression (inline-able)

--- Descend to the `call_expression` inside a short-form LHS or long-form
--- signature. The node may be the call directly, or a `typed_expression`
--- (`f(x)::T`) / `where_expression` (`f(x) where T`) wrapping it (possibly
--- nested, e.g. `f(x) where T where S`).
--- @param node TSNode|nil
--- @return TSNode|nil call
local function call_of(node)
    if not node then return nil end
    local t = node:type()
    if t == "call_expression" then return node end
    if t == "typed_expression" or t == "where_expression" then
        return call_of(node:named_child(0))
    end
    return nil
end

--- Extract the definition described by a node, or nil if it is not a function
--- definition (e.g. a plain `a = 1` assignment).
--- @param node TSNode
--- @param buf integer
--- @return JuliaDef|nil
local function parse_def(node, buf)
    local t = node:type()
    if t == "function_definition" then
        local sig, block
        for child in node:iter_children() do
            local ct = child:type()
            if ct == "signature" then sig = child
            elseif ct == "block" then block = child end
        end
        local call = sig and call_of(sig:named_child(0))
        local name = call and call:named_child(0)
        if not (sig and name) then return nil end
        local n = block and block:named_child_count() or 0
        local expr = block and block:named_child(0)
        --- @type JuliaDef
        local def = {
            range = { node:range() },
            form = "long",
            sig = gettext(sig, buf),
            callee = gettext(name, buf),
            single = n == 1,
            body = expr and gettext(expr, buf) or "",
        }
        return def
    elseif t == "assignment" then
        local lhs = node:named_child(0)
        local call = call_of(lhs)
        local name = call and call:named_child(0)
        if not (lhs and name) then return nil end
        local op = node:named_child(1)
        if not op or gettext(op, buf) ~= "=" then return nil end
        local rhs = node:named_child(2)
        --- @type JuliaDef
        local def = {
            range = { node:range() },
            form = "short",
            sig = gettext(lhs, buf),
            callee = gettext(name, buf),
            single = true,
            body = rhs and gettext(rhs, buf) or "",
        }
        return def
    end
    return nil
end

--- Ascend from the cursor to the nearest function definition, stopping before
--- any enclosing `macrocall_expression` so a `@inline` prefix is preserved.
--- @return JuliaDef|nil
local function def_under_cursor()
    return ts.ancestor(function(node) return parse_def(node, 0) end)
end

--- One indentation level for the current buffer.
--- @return string
local function indent_unit()
    if vim.bo.expandtab then return string.rep(" ", vim.fn.shiftwidth()) end
    return "\t"
end

--- Replacement lines for a definition rewritten to `target` form. The body is
--- indented relative to `base` (the def line's leading whitespace); the first
--- line carries no indent as it is inserted after the def's existing prefix.
--- @param def JuliaDef
--- @param target "short"|"long"
--- @param base string
--- @param unit string
--- @return string[]
local function rewrite_lines(def, target, base, unit)
    if target == "long" then
        local lines = { "function " .. def.sig }
        for _, l in ipairs(vim.split(def.body, "\n", { plain = true })) do
            lines[#lines + 1] = base .. unit .. l
        end
        lines[#lines + 1] = base .. "end"
        return lines
    end
    return vim.split(def.sig .. " = " .. def.body, "\n", { plain = true })
end

--- Rewrite a definition in place to `target` form (single buffer edit).
--- @param def JuliaDef
--- @param target "short"|"long"
local function apply(def, target)
    local srow, scol, erow, ecol = def.range[1], def.range[2], def.range[3], def.range[4]
    local first = vim.api.nvim_buf_get_lines(0, srow, srow + 1, false)[1] or ""
    local base = first:match("^%s*")
    vim.api.nvim_buf_set_text(0, srow, scol, erow, ecol,
        rewrite_lines(def, target, base, indent_unit()))
end

--- Toggle the function definition under the cursor between short and long form.
function M.toggle()
    local def = def_under_cursor()
    if not def then
        vim.notify("no function definition under cursor", vim.log.levels.WARN)
        return
    end
    local target = def.form == "short" and "long" or "short"
    if target == "short" and not def.single then
        vim.notify("body is not a single expression", vim.log.levels.WARN)
        return
    end
    apply(def, target)
end

--- Conform every same-name method in the buffer to the current form of the
--- definition under the cursor (the unchanged exemplar). Long→short is skipped
--- for multi-statement bodies, which cannot be inlined.
function M.conform()
    local exemplar = def_under_cursor()
    if not exemplar then
        vim.notify("no function definition under cursor", vim.log.levels.WARN)
        return
    end
    local target = exemplar.form
    local root = vim.treesitter.get_parser(0, "julia"):parse()[1]:root()
    local query = vim.treesitter.query.parse("julia",
        "[(function_definition) @d (assignment) @d]")

    local defs = {}
    for _, node in query:iter_captures(root, 0) do
        local def = parse_def(node, 0)
        if def and def.callee == exemplar.callee and def.form ~= target then
            defs[#defs + 1] = def
        end
    end
    -- Edit bottom-to-top so earlier defs' stored ranges stay valid.
    table.sort(defs, function(a, b) return a.range[1] > b.range[1] end)

    local skipped = {}
    local edited = 0
    for _, def in ipairs(defs) do
        if target == "short" and not def.single then
            skipped[#skipped + 1] = def.range[1] + 1
        else
            if edited > 0 then vim.cmd("silent! undojoin") end
            apply(def, target)
            edited = edited + 1
        end
    end

    if #skipped > 0 then
        table.sort(skipped)
        vim.notify(("kept %d multi-statement %s def(s) at line(s) %s"):format(
            #skipped, exemplar.callee, table.concat(skipped, ", ")),
            vim.log.levels.WARN)
    end
end

return M
