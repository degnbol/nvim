---@diagnostic disable: undefined-global
local python = require "docstring.lang.python"

--- Parse Python `src`, returning the first `function_definition` node.
local function first_function(src)
    local root = vim.treesitter.get_string_parser(src, "python"):parse()[1]:root()
    local q = vim.treesitter.query.parse("python", "(function_definition) @f")
    for _, node in q:iter_captures(root, src) do
        return node
    end
end

describe("docstring python signature_params", function()
    local function params_of(src)
        return python.signature_params(first_function(src), src)
    end

    it("drops self and keeps annotated + unannotated params in order", function()
        local params = params_of("def m(self, a, b: int):\n    pass")
        assert.are.same({
            { name = "a", type = nil },
            { name = "b", type = "int" },
        }, params)
    end)

    it("drops cls", function()
        local params = params_of("def m(cls, x):\n    pass")
        assert.are.same({ { name = "x", type = nil } }, params)
    end)

    it("keeps *args and **kwargs with their stars", function()
        local params = params_of("def f(a, *args, **kwargs):\n    pass")
        assert.are.same({
            { name = "a", type = nil },
            { name = "*args", type = nil },
            { name = "**kwargs", type = nil },
        }, params)
    end)

    it("skips / and * separators", function()
        local params = params_of("def f(a, /, b, *, c):\n    pass")
        assert.are.same({
            { name = "a", type = nil },
            { name = "b", type = nil },
            { name = "c", type = nil },
        }, params)
    end)

    it("ignores defaults for the name and keeps the type", function()
        local params = params_of("def f(a=3, b: str = 'x'):\n    pass")
        assert.are.same({
            { name = "a", type = nil },
            { name = "b", type = "str" },
        }, params)
    end)

    it("returns an empty list for an empty signature", function()
        assert.are.same({}, params_of("def f():\n    pass"))
    end)
end)

describe("docstring python locate_docstring", function()
    it("returns dedented text with the first line inline", function()
        local src = table.concat({
            "def f(a):",
            '    """Summary.',
            "",
            "    Args:",
            "        a: A.",
            '    """',
            "    return a",
        }, "\n")
        local node, text = python.locate_docstring(first_function(src), src)
        assert.is_not_nil(node)
        assert.are.equal("Summary.\n\nArgs:\n    a: A.", text)
    end)

    it("returns nil for a function with no docstring", function()
        local src = "def f(a):\n    return a"
        local node = python.locate_docstring(first_function(src), src)
        assert.is_nil(node)
    end)

    it("rejects an f-string opening statement (not a docstring)", function()
        local src = 'def f(a):\n    f"""Value {a}."""\n    return a'
        local node = python.locate_docstring(first_function(src), src)
        assert.is_nil(node)
    end)
end)

describe("docstring python body_row", function()
    it("returns the first statement row of a normal indented body", function()
        local src = "def f(a):\n    return a"
        assert.are.equal(1, python.body_row(first_function(src), src))
    end)

    it("rejects a one-line body", function()
        local src = "def f(a): return a"
        assert.is_nil(python.body_row(first_function(src), src))
    end)

    it("rejects an inline body after a multi-line signature", function()
        local src = "def f(\n    a,\n): return a"
        assert.is_nil(python.body_row(first_function(src), src))
    end)
end)

describe("docstring python locate_target", function()
    it("finds the enclosing function from an inner node", function()
        local src = "def f(a):\n    return a"
        local root = vim.treesitter.get_string_parser(src, "python"):parse()[1]:root()
        local q = vim.treesitter.query.parse("python", "(identifier) @i")
        local inner
        for _, node in q:iter_captures(root, src) do
            inner = node -- last identifier: the `a` in `return a`
        end
        local target = assert(python.locate_target(inner))
        assert.are.equal("function_definition", target:type())
    end)
end)
