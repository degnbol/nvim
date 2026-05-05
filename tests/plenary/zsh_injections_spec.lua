---@diagnostic disable: undefined-global
-- Tests for queries/zsh/injections.scm

vim.treesitter.language.register("zsh", "sh.zsh")
-- Mirror plugin/treesitter.lua: bash/sh injections resolve to zsh parser
-- (only zsh.so is installed; bash.so is not). Also load the `#trim!` query
-- directive from plugin/treesitter.lua. plugin/* isn't reliably sourced
-- before PlenaryBustedFile runs, so source it explicitly.
vim.treesitter.language.register("zsh", "bash")
vim.treesitter.language.register("zsh", "sh")
vim.cmd.source(vim.fn.getcwd() .. "/plugin/treesitter.lua")

--- Parse source as zsh and return a table of injected language names.
--- Each entry is { lang = "...", text = "..." } for the injected region.
--- Recurses into nested injections (e.g. zsh → vim → lua).
local function get_injections(src)
    local buf = vim.api.nvim_create_buf(true, false)
    vim.api.nvim_set_current_buf(buf)
    local lines = vim.split(src, "\n", { plain = true })
    vim.api.nvim_buf_set_lines(buf, 0, -1, false, lines)
    vim.bo[buf].filetype = "sh.zsh"
    local parser = vim.treesitter.get_parser(buf)
    parser:parse(true)
    local results = {}
    local function collect(p)
        for lang, child in pairs(p:children()) do
            for _, tree in ipairs(child:trees()) do
                local root = tree:root()
                results[#results + 1] = {
                    lang = lang,
                    text = vim.treesitter.get_node_text(root, buf),
                }
            end
            collect(child)
        end
    end
    collect(parser)
    pcall(vim.api.nvim_buf_delete, buf, { force = true })
    return results
end

--- Return all injection entries for a specific language.
local function injections_for(src, lang)
    local all = get_injections(src)
    local filtered = {}
    for _, inj in ipairs(all) do
        if inj.lang == lang then
            filtered[#filtered + 1] = inj
        end
    end
    return filtered
end

--- Assert at least one injection of the given language contains the expected text.
local function assert_injection(src, expected_lang, expected_text)
    local inj = injections_for(src, expected_lang)
    assert.is_true(#inj >= 1,
        "expected at least 1 " .. expected_lang .. " injection, got 0")
    for _, entry in ipairs(inj) do
        if entry.text == expected_text then return end
    end
    local got = {}
    for _, e in ipairs(inj) do got[#got + 1] = vim.inspect(e.text) end
    assert.is_true(false,
        "no " .. expected_lang .. " injection with text "
        .. vim.inspect(expected_text)
        .. ". Got: " .. table.concat(got, ", "))
end

local function assert_no_injection(src, lang)
    local inj = injections_for(src, lang)
    assert.are.equal(0, #inj,
        "expected no " .. lang .. " injections, got " .. #inj)
end

describe("zsh injections", function()
    describe("miller", function()
        it("injects miller into single-quoted put string", function()
            assert_injection(
                "mlr put '$x = 1'",
                "miller", "$x = 1"
            )
        end)

        it("injects miller into single-quoted filter string", function()
            assert_injection(
                "mlr filter '$score > 0.5'",
                "miller", "$score > 0.5"
            )
        end)

        it("does not inject for non-DSL verbs", function()
            assert_no_injection("mlr head -n 5", "miller")
        end)

        it("does not inject double-quoted strings", function()
            -- double-quoted strings have zsh variable expansion
            assert_no_injection('mlr put "$x = 1"', "miller")
        end)
    end)

    describe("python", function()
        it("injects python into single-quoted python3 -c", function()
            assert_injection(
                "python3 -c 'print(1)'",
                "python", "print(1)"
            )
        end)

        it("injects python into double-quoted python3 -c", function()
            assert_injection(
                'python3 -c "print(1)"',
                "python", "print(1)"
            )
        end)

        it("injects python into single-quoted python -c", function()
            assert_injection(
                "python -c 'import sys; print(sys.argv)'",
                "python", "import sys; print(sys.argv)"
            )
        end)

        it("works with flags before -c", function()
            assert_injection(
                "python3 -u -c 'print(1)'",
                "python", "print(1)"
            )
        end)

        it("does not inject without -c flag", function()
            assert_no_injection("python3 script.py", "python")
        end)

        it("does not inject for unrelated commands", function()
            assert_no_injection("echo 'print(1)'", "python")
        end)
    end)

    describe("julia", function()
        it("injects julia into single-quoted julia -e", function()
            assert_injection("julia -e 'println(1)'", "julia", "println(1)")
        end)

        it("injects julia into double-quoted julia -e", function()
            assert_injection('julia -e "println(1)"', "julia", "println(1)")
        end)

        it("works with flags before -e", function()
            assert_injection(
                "julia --threads=4 -e 'println(1)'",
                "julia", "println(1)"
            )
        end)

        it("does not inject without -e flag", function()
            assert_no_injection("julia script.jl", "julia")
        end)

        it("does not inject for other commands", function()
            assert_no_injection("echo -e 'println(1)'", "julia")
        end)
    end)

    describe("R", function()
        it("injects r into single-quoted Rscript -e", function()
            assert_injection("Rscript -e 'print(1)'", "r", "print(1)")
        end)

        it("injects r into double-quoted Rscript -e", function()
            assert_injection('Rscript -e "print(1)"', "r", "print(1)")
        end)

        it("injects r into R -e", function()
            assert_injection("R -e 'print(1)'", "r", "print(1)")
        end)

        it("works with flags before -e", function()
            assert_injection(
                "R --no-save -e 'print(1)'", "r", "print(1)"
            )
        end)

        it("does not inject without -e flag", function()
            assert_no_injection("Rscript script.R", "r")
        end)

        it("does not inject for other commands", function()
            assert_no_injection("echo -e 'print(1)'", "r")
        end)
    end)

    describe("javascript", function()
        it("injects javascript into single-quoted node -e", function()
            assert_injection(
                "node -e 'console.log(1)'",
                "javascript", "console.log(1)"
            )
        end)

        it("injects javascript into double-quoted node -e", function()
            assert_injection(
                'node -e "console.log(1)"',
                "javascript", "console.log(1)"
            )
        end)

        it("does not inject without -e flag", function()
            assert_no_injection("node app.js", "javascript")
        end)

        it("does not inject for other commands", function()
            assert_no_injection("echo -e 'console.log(1)'", "javascript")
        end)
    end)

    describe("shell -c", function()
        it("injects zsh into single-quoted zsh -c", function()
            assert_injection(
                "zsh -c 'echo hello'",
                "zsh", "echo hello"
            )
        end)

        it("injects zsh into double-quoted zsh -c", function()
            assert_injection(
                'zsh -c "echo hello"',
                "zsh", "echo hello"
            )
        end)

        it("injects zsh for bash -c", function()
            assert_injection("bash -c 'echo x'", "zsh", "echo x")
        end)

        it("injects zsh for sh -c", function()
            assert_injection("sh -c 'echo x'", "zsh", "echo x")
        end)

        it("injects zsh into each fragment of concatenated raw_strings", function()
            -- nvim's injection processor trims leading whitespace from
            -- injected content, so use content without leading spaces.
            local inj = injections_for(
                "zsh -c 'prefix='$ROOT';suffix'", "zsh")
            local texts = {}
            for _, entry in ipairs(inj) do texts[#texts + 1] = entry.text end
            assert.is_true(vim.tbl_contains(texts, "prefix="))
            assert.is_true(vim.tbl_contains(texts, ";suffix"))
        end)

        it("does not inject without -c flag", function()
            assert_no_injection("zsh script.sh", "zsh")
        end)

        it("does not inject for unrelated commands", function()
            assert_no_injection("echo -c 'hi'", "zsh")
        end)
    end)

    describe("jq", function()
        it("injects jq into single-quoted jq filter", function()
            assert_injection("jq '.a'", "jq", ".a")
        end)

        it("injects jq into double-quoted jq filter", function()
            assert_injection('jq ".a"', "jq", ".a")
        end)

        it("injects jq with flags before filter", function()
            assert_injection("jq -r '.items[]'", "jq", ".items[]")
        end)

        it("injects jq for gojq", function()
            assert_injection("gojq '.foo // \"x\"'", "jq", ".foo // \"x\"")
        end)

        it("injects the filter even when --arg is also passed", function()
            -- `jq --arg name "Alice" '.user = $name'`: the trailing filter
            -- must be highlighted. The `Alice` value is also injected (it
            -- parses as a harmless jq string literal); not worth the query
            -- complexity to suppress it.
            local inj = injections_for(
                "jq --arg name \"Alice\" '.user = $name'", "jq")
            local texts = {}
            for _, e in ipairs(inj) do texts[#texts + 1] = e.text end
            assert.is_true(vim.tbl_contains(texts, ".user = $name"))
        end)

        it("does not inject for unrelated commands", function()
            assert_no_injection("echo '.a'", "jq")
        end)
    end)

    describe("nvim lua", function()
        it("injects lua into single-quoted nvim -c 'lua ...'", function()
            assert_injection(
                "nvim -c 'lua print(1)'",
                "lua", "print(1)"
            )
        end)

        it("injects lua into double-quoted nvim -c", function()
            assert_injection(
                'nvim -c "lua print(1)"',
                "lua", "print(1)"
            )
        end)

        it("works with --headless before -c", function()
            assert_injection(
                "nvim --headless -c 'lua print(1)'",
                "lua", "print(1)"
            )
        end)

        it("works with vim command name", function()
            assert_injection(
                "vim -c 'lua print(1)'",
                "lua", "print(1)"
            )
        end)

        -- Multiline single-quoted lua verified manually via get_string_parser
        -- (buffer-based test is fragile due to #offset! across line boundaries)

        it("does not inject non-lua -c strings", function()
            assert_no_injection("nvim -c 'set number'", "lua")
        end)

        it("only injects the lua -c string, not others", function()
            local inj = injections_for(
                "nvim --headless -c 'lua print(1)' -c 'qa'", "lua")
            -- Should have the lua injection but not inject 'qa'
            local texts = {}
            for _, entry in ipairs(inj) do texts[#texts + 1] = entry.text end
            assert.is_true(vim.tbl_contains(texts, "print(1)"))
            assert.is_false(vim.tbl_contains(texts, "qa"))
        end)

        it("injects vim into non-lua -c args", function()
            -- The vim parser gets the whole arg; vim's own injections
            -- handle :lua / :python / etc. inside.
            assert_injection("nvim -c 'set number'", "vim", "set number")
        end)

        it("injects vim for double-quoted -c args", function()
            assert_injection('nvim -c "set number"', "vim", "set number")
        end)

        it("injects vim for + form", function()
            assert_injection("nvim +'set number'", "vim", "set number")
        end)

        it("injects vim for --cmd form", function()
            assert_injection("nvim --cmd 'set number'", "vim", "set number")
        end)

        it("injects lua for multi-line --cmd 'lua\\n...\\n'", function()
            assert_injection("nvim --cmd 'lua\nprint(1)\n'", "lua", "print(1)")
        end)

        it("injects lua directly for multi-line raw_string -c 'lua\\n...\\n'",
            function()
                local src = "nvim -c 'lua\nprint(1)\n'"
                assert_injection(src, "lua", "print(1)")
            end)

        it("injects lua directly for multi-line double-quoted -c \"lua\\n...\\n\"",
            function()
                local src = 'nvim -c "lua\nprint(1)\n"'
                assert_injection(src, "lua", "print(1)")
            end)

        it("injects lua directly for multi-line + form +'lua\\n...\\n'",
            function()
                local src = "nvim +'lua\nprint(1)\n'"
                assert_injection(src, "lua", "print(1)")
            end)
    end)

    describe("heredoc by file-redirect extension", function()
        local function heredoc(dest, tag, body)
            tag = tag or "EOF"
            local close = tag:gsub("['\"]", "")
            return table.concat({
                "cat > " .. dest .. " <<" .. tag,
                body,
                close,
            }, "\n")
        end

        it("injects lua for .lua redirect target", function()
            assert_injection(heredoc("/tmp/a.lua", "'EOF'", "local x = 1"),
                "lua", "local x = 1")
        end)

        it("injects python for .py redirect target", function()
            assert_injection(heredoc("/tmp/a.py", "'EOF'", "print(1)"),
                "python", "print(1)")
        end)

        it("injects julia for .jl redirect target", function()
            assert_injection(heredoc("/tmp/a.jl", "EOF", "println(1)"),
                "julia", "println(1)")
        end)

        it("injects r for .R redirect target", function()
            assert_injection(heredoc("/tmp/a.R", "'EOF'", "x <- 1"),
                "r", "x <- 1")
        end)

        it("injects javascript for .js redirect target", function()
            assert_injection(heredoc("/tmp/a.js", "'EOF'", "const x = 1;"),
                "javascript", "const x = 1;")
        end)

        it("injects zsh (via bash alias) for .sh redirect target", function()
            assert_injection(heredoc("/tmp/a.sh", "'EOF'", "echo hi"),
                "zsh", "echo hi")
        end)

        it("injects zsh for .zsh redirect target", function()
            assert_injection(heredoc("/tmp/a.zsh", "'EOF'",
                    "typeset -A map"),
                "zsh", "typeset -A map")
        end)

        it("handles quoted redirect destinations", function()
            assert_injection(
                'cat > "/tmp/a.lua" <<EOF\nlocal x = 1\nEOF',
                "lua", "local x = 1")
            assert_injection(
                "cat > '/tmp/a.lua' <<EOF\nlocal x = 1\nEOF",
                "lua", "local x = 1")
        end)

        it("handles heredoc before file-redirect", function()
            assert_injection(
                "cat <<EOF > /tmp/a.lua\nlocal x = 1\nEOF",
                "lua", "local x = 1")
        end)

        it("still supports heredoc-tag-as-language (base query)", function()
            assert_injection(
                "cat <<LUA\nlocal x = 1\nLUA",
                "lua", "local x = 1")
        end)

        it("does not inject for unknown extension", function()
            assert_no_injection(
                "cat > /tmp/a.xyz <<EOF\nblah\nEOF",
                "lua")
        end)

        it(".bash maps to bash (zsh alias), not sh", function()
            -- .bash uses the bash/zsh parser, not anything else
            assert_injection(heredoc("/tmp/a.bash", "'EOF'", "echo hi"),
                "zsh", "echo hi")
        end)

        it(".mjs and .cjs also inject javascript", function()
            assert_injection(heredoc("/tmp/a.mjs", "'EOF'", "export const x = 1;"),
                "javascript", "export const x = 1;")
            assert_injection(heredoc("/tmp/a.cjs", "'EOF'", "module.exports = 1;"),
                "javascript", "module.exports = 1;")
        end)
    end)
end)
