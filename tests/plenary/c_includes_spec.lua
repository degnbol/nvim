---@diagnostic disable: undefined-global
local c_includes = require("utils/c_includes")

describe("c_includes.includes", function()
    it("splits angle and quoted headers, in order, each once", function()
        local lines = {
            "#include <stdio.h>",
            '#include "a.h"',
            "  #  include   <glib/gtypes.h>",
            "#include<stdio.h>",
            '#include "a.h" // again',
            "int x; // #include <not.h>",
            '#include "b/c.h"',
        }
        assert.are.same(
            { angle = { "stdio.h", "glib/gtypes.h" }, quoted = { "a.h", "b/c.h" } },
            c_includes.includes(lines)
        )
    end)

    it("reports nothing for no directives", function()
        assert.are.same({ angle = {}, quoted = {} }, c_includes.includes({ "int main(void) { return 0; }" }))
    end)
end)

describe("c_includes.parse_search_dirs", function()
    it("keeps the angle search dirs without framework dirs", function()
        local output = table.concat({
            "clang version 23.1.1",
            "ignoring nonexistent directory \"/x\"",
            '#include "..." search starts here:',
            " /quoted/only",
            "#include <...> search starts here:",
            " /opt/homebrew/Cellar/llvm/23.1.1_1/lib/clang/23/include",
            " /Library/Developer/CommandLineTools/SDKs/MacOSX26.sdk/usr/include",
            " /Library/Developer/CommandLineTools/SDKs/MacOSX26.sdk/System/Library/Frameworks (framework directory)",
            "End of search list.",
            "# 1 \"<stdin>\"",
        }, "\n")
        assert.are.same({
            "/opt/homebrew/Cellar/llvm/23.1.1_1/lib/clang/23/include",
            "/Library/Developer/CommandLineTools/SDKs/MacOSX26.sdk/usr/include",
        }, c_includes.parse_search_dirs(output))
    end)
end)

describe("c_includes.search_dirs", function()
    it("reports the compiler's search dirs", function()
        if vim.fn.executable("clang") == 0 then return pending("clang not installed") end
        local dirs = require("vim._async").run(function() return c_includes.search_dirs("clang", "c") end):wait()
        assert.truthy(c_includes.find_header("stdio.h", dirs))
    end)
end)

describe("c_includes.find_header", function()
    local tmp = vim.fn.tempname()
    vim.fn.mkdir(vim.fs.joinpath(tmp, "b", "sub"), "p")
    vim.fn.writefile({}, vim.fs.joinpath(tmp, "b", "sub", "x.h"))

    it("returns the path in the first dir that holds the header", function()
        assert.are.equal(
            vim.fs.joinpath(tmp, "b", "sub", "x.h"),
            c_includes.find_header("sub/x.h", { vim.fs.joinpath(tmp, "a"), vim.fs.joinpath(tmp, "b") })
        )
    end)

    it("returns nil when no dir holds it", function()
        assert.is_nil(c_includes.find_header("x.h", { tmp }))
    end)
end)
