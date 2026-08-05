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

--- Derive the argparse dest of an `add_argument` call and the string token that
--- produced it. `kind` distinguishes how the token must be rewritten on rename:
--- `"positional"`/`"dest"` take the new name verbatim, `"flag"` re-adds the
--- leading dashes and maps `_`→`-`.
--- @param call TSNode  an `add_argument` `call` node
--- @param src string
--- @return string|nil dest
--- @return TSNode|nil token  the `string_content` node to rewrite
--- @return string|nil kind  "positional" | "flag" | "dest"
local function derive_dest(call, src)
    local args = call:field("arguments")[1]
    if not args then return end

    --- @type { node: TSNode, text: string }[]  option strings, quotes stripped
    local options = {}
    local has_dest = false
    local dest_content = nil
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
                local value = child:field("value")[1]
                dest_content = value and literal_content(value)
            end
        end
    end

    -- A `dest=` kwarg overrides the option strings, decoupling them from the
    -- attribute. Its presence — even non-literal — means the option strings must
    -- not be touched; bail unless we can rewrite the literal dest value itself.
    if has_dest then
        if not dest_content then return end
        return vim.treesitter.get_node_text(dest_content, src), dest_content, "dest"
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
    local dest = flag.text:gsub("^%-+", ""):gsub("%-", "_")
    return dest, flag.node, "flag"
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
