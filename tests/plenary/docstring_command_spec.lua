---@diagnostic disable: undefined-global
local command = require "docstring.command"

--- Run `:Docstring` `times` times on a buffer of `src` with the cursor on line
--- `cursor_row` (1-indexed), returning the resulting lines.
local function run(src, cursor_row, times)
    local buf = vim.api.nvim_create_buf(false, true)
    vim.api.nvim_buf_set_lines(buf, 0, -1, false, vim.split(src, "\n", { plain = true }))
    vim.bo[buf].filetype = "python"
    vim.api.nvim_set_current_buf(buf)
    vim.treesitter.get_parser(buf, "python"):parse()
    vim.api.nvim_win_set_cursor(0, { cursor_row, 4 })
    for _ = 1, times or 1 do
        command.reconcile_at_cursor()
        vim.treesitter.get_parser(buf, "python"):parse()
    end
    local lines = vim.api.nvim_buf_get_lines(buf, 0, -1, false)
    vim.api.nvim_buf_delete(buf, { force = true })
    return lines
end

describe("docstring command", function()
    it("adds a missing param to an existing docstring, re-indented", function()
        local src = table.concat({
            "def f(a, b):",
            '    """Summary.',
            "",
            "    Args:",
            "        a: First.",
            '    """',
            "    return a",
        }, "\n")
        assert.are.same({
            "def f(a, b):",
            '    """Summary.',
            "",
            "    Args:",
            "        a: First.",
            "        b:",
            '    """',
            "    return a",
        }, run(src, 1))
    end)

    it("reorders doc entries to signature order", function()
        local src = table.concat({
            "def f(b, a):",
            '    """S.',
            "",
            "    Args:",
            "        a: First.",
            "        b: Second.",
            '    """',
            "    return a",
        }, "\n")
        assert.are.same({
            "def f(b, a):",
            '    """S.',
            "",
            "    Args:",
            "        b: Second.",
            "        a: First.",
            '    """',
            "    return a",
        }, run(src, 1))
    end)

    it("inserts a fresh docstring when the function has none", function()
        local src = table.concat({
            "def g(x, y):",
            "    return x",
        }, "\n")
        assert.are.same({
            "def g(x, y):",
            '    """Args:',
            "        x:",
            "        y:",
            '    """',
            "    return x",
        }, run(src, 1))
    end)

    it("is idempotent: a second run leaves the buffer unchanged", function()
        local src = table.concat({
            "def f(a, b):",
            '    """Summary.',
            "",
            "    Args:",
            "        a: First.",
            '    """',
            "    return a",
        }, "\n")
        assert.are.same(run(src, 1, 1), run(src, 1, 2))
    end)
end)
