-- Python treesitter locators for the docstring library: find the function under
-- the cursor, extract its ordered signature params (with the Google/Python
-- filtering rules), and locate its docstring node + dedented text. The forward
-- (code→structure) half of the pipeline; styles handle the prose half.

local ts = require "utils/treesitter"

local M = {}

--- One signature parameter. `type` is the annotation text, nil when unannotated.
--- @class docstring.SigParam
--- @field name string
--- @field type? string

--- Params documented under a different name (or not at all) in `Args:`.
--- @type table<string, true>
local SKIP_NAMES = { self = true, cls = true }

--- The parameter name a node declares, with splat stars: `identifier` →
--- `"x"`, `list_splat_pattern` → `"*args"`, `dictionary_splat_pattern` →
--- `"**kwargs"`; nil for anything else.
--- @param node TSNode
--- @param src integer|string  buffer number or source text
--- @return string|nil
local function name_of(node, src)
    local kind = node:type()
    if kind == "identifier" then
        return vim.treesitter.get_node_text(node, src)
    elseif kind == "list_splat_pattern" then
        return "*" .. vim.treesitter.get_node_text(assert(node:named_child(0)), src)
    elseif kind == "dictionary_splat_pattern" then
        return "**" .. vim.treesitter.get_node_text(assert(node:named_child(0)), src)
    end
end

--- The `{ name, type }` a `parameters` child declares, or nil for a node that is
--- not a documentable parameter (`/` and `*` separators, comments).
--- @param node TSNode
--- @param src integer|string  buffer number or source text
--- @return docstring.SigParam|nil
local function param_of(node, src)
    local kind = node:type()
    local name, type_node
    if kind == "identifier" or kind == "list_splat_pattern" or kind == "dictionary_splat_pattern" then
        name = name_of(node, src)
    elseif kind == "typed_parameter" then
        name = name_of(assert(node:named_child(0)), src)
        type_node = node:field("type")[1]
    elseif kind == "default_parameter" then
        name = name_of(node:field("name")[1], src)
    elseif kind == "typed_default_parameter" then
        name = name_of(node:field("name")[1], src)
        type_node = node:field("type")[1]
    end
    if not name then return end
    --- @type docstring.SigParam
    local param = { name = name, type = type_node and vim.treesitter.get_node_text(type_node, src) }
    return param
end

--- The `function_definition` enclosing `node` (default: node under the cursor),
--- or nil if the cursor is not in a function.
--- @param node TSNode|nil
--- @return TSNode|nil
function M.locate_target(node)
    return ts.ancestor("function_definition", node)
end

--- The ordered documentable parameters of a function, applying the Google/Python
--- filtering rules: `self`/`cls` dropped, `/` and `*` separators skipped,
--- `*args`/`**kwargs` kept with their stars, defaults ignored for the name.
--- @param target TSNode  a `function_definition` node
--- @param src integer|string  buffer number or source text
--- @return docstring.SigParam[] params
function M.signature_params(target, src)
    local params = {}
    local parameters = target:field("parameters")[1]
    if not parameters then return params end
    for child in parameters:iter_children() do
        local param = param_of(child, src)
        if param and not SKIP_NAMES[param.name] then
            params[#params + 1] = param
        end
    end
    return params
end

--- Strip a docstring's common indentation (inspect.cleandoc semantics): the
--- first line is left-stripped, subsequent lines lose the minimum indent among
--- them, and leading/trailing blank lines are dropped. Yields the dedented text
--- the style parser expects; the buffer writer re-applies the body indent.
--- @param text string  raw string content, between the quotes
--- @return string
local function cleandoc(text)
    local lines = vim.split(text, "\n", { plain = true })
    local margin = nil
    for i = 2, #lines do
        local line = lines[i]
        if not line:match("^%s*$") then
            local indent = #line:match("^%s*")
            if not margin or indent < margin then margin = indent end
        end
    end
    margin = margin or 0
    lines[1] = lines[1]:match("^%s*(.*)$")
    for i = 2, #lines do
        lines[i] = lines[i]:match("^%s*$") and "" or lines[i]:sub(margin + 1)
    end
    while #lines > 0 and lines[1] == "" do table.remove(lines, 1) end
    while #lines > 0 and lines[#lines] == "" do table.remove(lines) end
    return table.concat(lines, "\n")
end

--- The row a fresh docstring should be inserted at (0-indexed): the line of the
--- function body's first statement, which the new docstring precedes. nil unless
--- only whitespace precedes that statement on its line — an inline body
--- (`def f(): return a`, or a `): return a` after a multi-line signature) cannot
--- take a docstring line without reflowing, so it is rejected.
--- @param target TSNode  a `function_definition` node
--- @param src integer|string  buffer number or source text
--- @return integer|nil row
function M.body_row(target, src)
    local body = target:field("body")[1]
    if not body or body:type() ~= "block" then return end
    local first = body:named_child(0)
    if not first then return end
    local row, col = first:range()
    local line
    if type(src) == "number" then
        line = vim.api.nvim_buf_get_lines(src, row, row + 1, true)[1]
    else
        line = vim.split(src --[[@as string]], "\n", { plain = true })[row + 1]
    end
    if line:sub(1, col):match("%S") then return end
    return row
end

--- The docstring string node of a function and its dedented text, or nil for a
--- function whose body does not open with a plain string literal. An f-string
--- (`interpolation` child) is not a docstring — CPython sets `__doc__` to nil for
--- it — so it is rejected. Implicit concatenation (`"a" "b"`, a
--- `concatenated_string` node) is a docstring in CPython but is not handled here;
--- it reads as no docstring, a known limitation.
--- @param target TSNode  a `function_definition` node
--- @param src integer|string  buffer number or source text
--- @return TSNode|nil node
--- @return string|nil text  dedented docstring content (no quotes)
function M.locate_docstring(target, src)
    local body = target:field("body")[1]
    if not body then return end
    local first = body:named_child(0)
    if not first or first:type() ~= "expression_statement" then return end
    local string_node = first:named_child(0)
    if not string_node or string_node:type() ~= "string" then return end
    local content = nil
    for child in string_node:iter_children() do
        if child:type() == "interpolation" then return end
        if child:type() == "string_content" then content = child end
    end
    local text = content and vim.treesitter.get_node_text(content, src) or ""
    return string_node, cleandoc(text)
end

return M
