---@diagnostic disable: undefined-global
local util = require "utils/init"

describe("script_command", function()
    it("executes via shebang when no interpreter is given", function()
        assert.are.equal([[cd '/tmp/proj' && './s.py']], util.script_command('/tmp/proj/s.py'))
    end)

    it("prefixes the interpreter when one is given", function()
        assert.are.equal([[cd '/tmp/proj' && python 's.py']],
            util.script_command('/tmp/proj/s.py', 'python'))
    end)

    it("keeps a multi-word interpreter unquoted", function()
        assert.are.equal([[cd '/tmp/proj' && uv run 's.py']],
            util.script_command('/tmp/proj/s.py', 'uv run'))
    end)

    -- fnameescape leaves these to the shell; shellescape must not.
    it("quotes shell metacharacters in the directory", function()
        assert.are.equal([[cd '/tmp/a;b (c)' && './s.py']], util.script_command('/tmp/a;b (c)/s.py'))
    end)

    it("quotes shell metacharacters in the basename", function()
        assert.are.equal([[cd '/tmp/proj' && './a;b (c).py']],
            util.script_command('/tmp/proj/a;b (c).py'))
    end)

    -- `:!` expands unescaped ! % # regardless of shell quoting.
    it("escapes the characters :! expands", function()
        assert.are.equal([[cd '/tmp/w\!\%\#$q' && './s.py']], util.script_command('/tmp/w!%#$q/s.py'))
    end)
end)
