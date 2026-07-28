---@diagnostic disable: undefined-global
-- Tests for the $expansion repaint in after/queries/zsh/highlights.scm.
-- Inside interpreter injections (double-quoted strings, heredoc bodies) the
-- base zsh captures for $var are re-asserted at priority 101 so the outer
-- language wins over the injected language (default 100). See
-- queries/zsh/injections.scm for the injections themselves.

vim.treesitter.language.register("zsh", "sh.zsh")
-- minimal_init.lua only adds the config dir to rtp; its after/ subdir (where the
-- repaint query lives) is not searched unless appended explicitly.
vim.opt.runtimepath:append(vim.fn.getcwd() .. "/after")
-- Load the #inject-interp! / #inject-interp-cmd! directives: parsing injections
-- needs them, and plugin/* isn't reliably sourced before PlenaryBustedFile runs.
vim.cmd.source(vim.fn.getcwd() .. "/plugin/treesitter.lua")

--- Highlight captures landing on a node whose text equals `text`, parsing
--- `src` as zsh. Each entry is { group = "...", priority = <number|nil> }.
--- @param src string
--- @param text string
--- @return { group: string, priority: number|nil }[]
local function captures_on(src, text)
    local buf = vim.api.nvim_create_buf(true, false)
    vim.api.nvim_set_current_buf(buf)
    vim.api.nvim_buf_set_lines(buf, 0, -1, false,
        vim.split(src, "\n", { plain = true }))
    vim.bo[buf].filetype = "sh.zsh"
    local parser = vim.treesitter.get_parser(buf)
    parser:parse(true)
    local root = parser:trees()[1]:root()
    local query = vim.treesitter.query.get("zsh", "highlights")
    local captures = {}
    for id, node, meta in query:iter_captures(root, buf, 0, -1) do
        if vim.treesitter.get_node_text(node, buf) == text then
            captures[#captures + 1] = {
                group = query.captures[id],
                priority = tonumber(meta.priority),
            }
        end
    end
    pcall(vim.api.nvim_buf_delete, buf, { force = true })
    return captures
end

--- True if `captures` holds a capture of `group` at `priority`.
--- @param captures { group: string, priority: number|nil }[]
--- @param group string
--- @param priority number|nil
--- @return boolean
local function has_capture(captures, group, priority)
    for _, capture in ipairs(captures) do
        if capture.group == group and capture.priority == priority then
            return true
        end
    end
    return false
end

--- Assert the $var repaint fires over `src`: the `$` sigil is
--- @punctuation.special and the name is @variable, both at priority 101.
--- @param src string
--- @param name string the expansion's variable name node text
local function assert_repaint(src, name)
    assert.is_true(
        has_capture(captures_on(src, "$"), "punctuation.special", 101),
        "expected @punctuation.special@101 on $ sigil")
    assert.is_true(
        has_capture(captures_on(src, name), "variable", 101),
        "expected @variable@101 on " .. name)
end

describe("zsh $expansion repaint over injections", function()
    it("repaints $var in a double-quoted interpreter string", function()
        assert_repaint([[gnuplot -e "load '$ROOT/x'"]], "ROOT")
    end)

    it("repaints $var in a heredoc body", function()
        assert_repaint("gnuplot <<GP\nload '$ROOT/x'\nGP", "ROOT")
    end)

    it("repaints ${…} sigils and name in a double-quoted string", function()
        -- tree-sitter-zsh tokenises the `${` sigil as two 1-char tokens.
        local src = [[gnuplot -e "set title '${TITLE}'"]]
        assert.is_true(
            has_capture(captures_on(src, "$"), "punctuation.special", 101),
            "expected @punctuation.special@101 on $")
        assert.is_true(
            has_capture(captures_on(src, "}"), "punctuation.special", 101),
            "expected @punctuation.special@101 on }")
        assert.is_true(
            has_capture(captures_on(src, "TITLE"), "variable", 101),
            "expected @variable@101 on TITLE")
    end)

    it("does not repaint $var outside strings/heredocs", function()
        assert.is_false(
            has_capture(captures_on("echo $ROOT", "$"),
                "punctuation.special", 101),
            "bare $ should keep base priority, not 101")
        assert.is_false(
            has_capture(captures_on("echo $ROOT", "ROOT"), "variable", 101),
            "bare name should keep base priority, not 101")
    end)
end)
