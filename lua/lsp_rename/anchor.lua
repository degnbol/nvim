-- Resolve a cursor sitting on an argparse/Tap `add_argument` option string to
-- the class attribute that string configures, so a rename started there has a
-- real symbol to anchor on. The reverse of the handler-wrapper augment: that
-- adds string edits to a symbol rename; this starts the rename from a string.
--
-- Only Tap subclasses are handled: the derived dest maps to a declared class
-- attribute of the same name, which the language server tracks as a real symbol.
-- Plain argparse (`args.<dest>`) has no such declaration and is left untouched.

local argparse = require "lsp_rename.argparse"
local ts = require "utils/treesitter"

local M = {}

--- Whether a `class_definition` lists a base whose name is `Tap`, bare or dotted
--- (`Tap`, `tap.Tap`). An aliased import (`Tap as T`) is not recognised.
--- @param class TSNode
--- @param src string
--- @return boolean
local function subclasses_tap(class, src)
    local bases = class:field("superclasses")[1]
    if not bases then return false end
    for base in bases:iter_children() do
        if base:named() and vim.treesitter.get_node_text(base, src):match("[%w_]+$") == "Tap" then
            return true
        end
    end
    return false
end

--- The declaring identifier of a class attribute named `name`, searching the
--- direct statements of a `class_definition`'s body.
--- @param class TSNode
--- @param name string
--- @param src string
--- @return TSNode|nil
local function attribute(class, name, src)
    local body = class:field("body")[1]
    if not body then return end
    for stmt in body:iter_children() do
        if stmt:type() == "expression_statement" then
            local assign = stmt:named_child(0)
            if assign and assign:type() == "assignment" then
                local left = assign:field("left")[1]
                if left and left:type() == "identifier"
                    and vim.treesitter.get_node_text(left, src) == name then
                    return left
                end
            end
        end
    end
end

--- The position to redirect a rename to when the cursor sits on an `add_argument`
--- option string inside a `Tap` subclass: the declaring identifier of the class
--- attribute the string configures. nil when the cursor is not on such a string,
--- the enclosing class is not a `Tap` subclass, or it declares no matching
--- attribute.
--- @param bufnr integer
--- @param row integer  0-indexed cursor row
--- @param col integer  0-indexed cursor byte column
--- @return integer|nil arow  0-indexed
--- @return integer|nil acol  0-indexed byte column
function M.position(bufnr, row, col)
    local node = vim.treesitter.get_node({ bufnr = bufnr, pos = { row, col } })
    if not node then return end
    local src = table.concat(vim.api.nvim_buf_get_lines(bufnr, 0, -1, false), "\n")

    local dest = argparse.dest_at(node, src)
    if not dest then return end

    local class = ts.ancestor("class_definition", node)
    if not class or not subclasses_tap(class, src) then return end

    local left = attribute(class, dest, src)
    if not left then return end
    local arow, acol = left:range()
    return arow, acol
end

return M
