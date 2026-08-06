-- Matcher for argparse / Tap `add_argument` option strings. Given Python source
-- and a symbol rename (old → new), it finds `add_argument(...)` calls whose
-- argparse-derived dest equals `old` and returns edits over the *specific*
-- string token that produced the dest, so the CLI binding tracks the rename.
-- The dest-derivation rule (first `--long` else first `-x`, strip leading
-- dashes, `-`→`_`, unless `dest=` overrides; a positional's dest is its name) is
-- shared truth with the argparse diagnostic — see
-- notes/PLAN-argparse-diagnostics.md.
--
-- Only plain string literals are inspected: option strings and `dest=` values
-- built from f-strings, concatenation, or other expressions can't be resolved
-- statically, so a call with a non-literal `dest=` is left untouched entirely.
-- Assumes the default `-` prefix char; `prefix_chars=` is not read.

local ts = require "utils/treesitter"

local M = {}

--- @class lsp_rename.Edit  A replacement over a source range, byte columns, 0-indexed.
--- @field srow integer
--- @field scol integer
--- @field erow integer
--- @field ecol integer
--- @field newText string

--- The `string_content` child of a *plain* string literal — the text between
--- the quotes. Returns nil for a non-`string` node (e.g. `concatenated_string`),
--- an f-string (has an `interpolation` child, so no static value), or an empty
--- string (no `string_content` child).
--- @param node TSNode
--- @return TSNode|nil
local function literal_content(node)
    if node:type() ~= "string" then return end
    local content = nil
    for child in node:iter_children() do
        local kind = child:type()
        if kind == "interpolation" then return end
        if kind == "string_content" then content = child end
    end
    return content
end

--- The method/function name of a `call`'s function node: the trailing identifier
--- of an attribute (`p.add_argument`) or a bare identifier (`add_argument`).
--- @param fn TSNode
--- @param src string
--- @return string|nil
local function called_name(fn, src)
    local name = fn:type() == "attribute" and fn:field("attribute")[1] or fn
    if name and name:type() == "identifier" then
        return vim.treesitter.get_node_text(name, src)
    end
end

--- Collect an argument_list's plain-literal option strings and its `dest=` value.
--- @param args TSNode  an `argument_list` node
--- @param src string
--- @return { node: TSNode, text: string }[] options  option strings, quotes stripped
--- @return { node: TSNode, text: string }|nil dest  a literal `dest=` value, quotes stripped
--- @return boolean has_dest  whether a `dest=` kwarg is present, even if non-literal
local function collect_options(args, src)
    --- @type { node: TSNode, text: string }[]
    local options = {}
    local dest = nil
    local has_dest = false
    for child in args:iter_children() do
        local kind = child:type()
        if kind == "string" then
            local content = literal_content(child)
            if content then
                options[#options + 1] =
                    { node = content, text = vim.treesitter.get_node_text(content, src) }
            end
        elseif kind == "keyword_argument" then
            local name = vim.treesitter.get_node_text(child:field("name")[1], src)
            if name == "dest" then
                has_dest = true
                local content = literal_content(child:field("value")[1])
                if content then
                    dest = { node = content, text = vim.treesitter.get_node_text(content, src) }
                end
            end
        end
    end
    return options, dest, has_dest
end

--- Derive the argparse dest and the string token that produced it from collected
--- option strings. `kind` distinguishes how the token must be rewritten on
--- rename: `"positional"`/`"dest"` take the new name verbatim, `"flag"` re-adds
--- the leading dashes and maps `_`→`-`.
--- @param options { node: TSNode, text: string }[]
--- @param dest { node: TSNode, text: string }|nil  literal `dest=` value
--- @param has_dest boolean
--- @return string|nil dest
--- @return TSNode|nil token  the `string_content` node to rewrite
--- @return string|nil kind  "positional" | "flag" | "dest"
local function derive_dest_from(options, dest, has_dest)
    -- A `dest=` kwarg overrides the option strings, decoupling them from the
    -- attribute. Its presence — even non-literal — means the option strings must
    -- not be touched; bail unless we can rewrite the literal dest value itself.
    if has_dest then
        if not dest then return end
        return dest.text, dest.node, "dest"
    end

    local first = options[1]
    if not first then return end
    if first.text:sub(1, 1) ~= "-" then
        -- Positional (`"infiles"`): the name is the dest and its own token.
        return first.text, first.node, "positional"
    end

    -- Flag(s): dest comes from the first `--long`, else the first option string.
    local flag = first
    for _, option in ipairs(options) do
        if option.text:sub(1, 2) == "--" then
            flag = option
            break
        end
    end
    local name = flag.text:gsub("^%-+", ""):gsub("%-", "_")
    return name, flag.node, "flag"
end

--- @param call TSNode  an `add_argument` `call` node
--- @param src string
--- @return string|nil dest
--- @return TSNode|nil token
--- @return string|nil kind
local function derive_dest(call, src)
    local args = call:field("arguments")[1]
    if not args then return end
    return derive_dest_from(collect_options(args, src))
end

--- Whether `node`'s start position lies within `outer`'s range (start-inclusive,
--- end-exclusive). Coordinates are 0-indexed byte (row, col).
--- @param node TSNode
--- @param outer TSNode
--- @return boolean
local function starts_within(node, outer)
    local nr, nc = node:range()
    local sr, sc, er, ec = outer:range()
    local after_start = nr > sr or (nr == sr and nc >= sc)
    local before_end = nr < er or (nr == er and nc < ec)
    return after_start and before_end
end

--- Edits that retarget `add_argument` option strings from `old` to `new`.
--- @param src string  Python source
--- @param old string  current symbol name
--- @param new string  new symbol name
--- @return lsp_rename.Edit[] edits
function M.edits(src, old, new)
    local root = vim.treesitter.get_string_parser(src, "python"):parse()[1]:root()
    local calls = vim.treesitter.query.parse("python", "(call) @call")

    --- @type lsp_rename.Edit[]
    local edits = {}
    for _, call in calls:iter_captures(root, src) do
        local fn = call:field("function")[1]
        if fn and called_name(fn, src) == "add_argument" then
            local dest, token, kind = derive_dest(call, src)
            if dest == old and token then
                local newText = new
                if kind == "flag" then
                    local dashes = vim.treesitter.get_node_text(token, src):match("^%-+")
                    newText = dashes .. new:gsub("_", "-")
                end
                local srow, scol, erow, ecol = token:range()
                edits[#edits + 1] =
                    { srow = srow, scol = scol, erow = erow, ecol = ecol, newText = newText }
            end
        end
    end
    return edits
end

--- The argparse dest configured by the option string (or literal `dest=` value)
--- under `node`, or nil if `node` is not inside such a string. The nearest
--- enclosing `call` must be `add_argument`; a string in any other call bails
--- rather than climbing to an outer `add_argument`.
--- @param node TSNode  the node under the cursor
--- @param src string
--- @return string|nil dest
function M.dest_at(node, src)
    local call = ts.ancestor("call", node)
    if not call then return end
    local fn = call:field("function")[1]
    if not fn or called_name(fn, src) ~= "add_argument" then return end
    local args = call:field("arguments")[1]
    if not args then return end

    local options, dest, has_dest = collect_options(args, src)

    --- @type TSNode[]  the string literals a cursor may sit in to trigger a rename
    local strings = {}
    for _, option in ipairs(options) do
        strings[#strings + 1] = option.node:parent()
    end
    if dest then
        strings[#strings + 1] = dest.node:parent()
    end

    for _, string_node in ipairs(strings) do
        if starts_within(node, string_node) then
            return (derive_dest_from(options, dest, has_dest))
        end
    end
end

--- The dest of the first `add_argument` call in `src`, or nil if there is none.
--- Test seam for the shared dest-derivation fixtures.
--- @param src string
--- @return string|nil
function M.dest_of(src)
    local root = vim.treesitter.get_string_parser(src, "python"):parse()[1]:root()
    local calls = vim.treesitter.query.parse("python", "(call) @call")
    for _, call in calls:iter_captures(root, src) do
        local fn = call:field("function")[1]
        if fn and called_name(fn, src) == "add_argument" then
            return (derive_dest(call, src))
        end
    end
end

return M
