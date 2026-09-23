---@diagnostic disable: undefined-global
local stub_patches = require "autocmds/stub_patches"

--- Write text to a path, creating parent directories.
local function write(path, text)
    vim.fn.mkdir(vim.fs.dirname(path), "p")
    vim.fn.writefile(vim.split(text, "\n"), path)
end

describe("stub_patches.imports", function()
    local function imports(lines)
        local buf = vim.api.nvim_create_buf(false, true)
        vim.api.nvim_buf_set_lines(buf, 0, -1, false, lines)
        return stub_patches.imports(buf)
    end

    it("finds each package of a multi-name import", function()
        assert.are.same({ a = { line = 0, character = 7 }, b = { line = 0, character = 10 } },
            imports { "import a, b.c" })
    end)

    it("finds the top-level package of a from-import", function()
        assert.are.same({ a = { line = 0, character = 5 } }, imports { "from a.b import c" })
    end)

    it("finds an aliased import", function()
        assert.are.same({ a = { line = 1, character = 11 } },
            imports { "try:", "    import a.b as x" })
    end)

    it("skips relative imports", function()
        assert.are.same({}, imports { "from . import a", "from .b import c" })
    end)

    it("keeps the first of duplicate imports", function()
        assert.are.same({ a = { line = 0, character = 7 } },
            imports { "import a", "from a import b" })
    end)
end)

describe("stub_patches.env_of", function()
    it("takes the first component below site-packages", function()
        local env = stub_patches.env_of("/p/.venv/lib/python3.12/site-packages/numpy/__init__.pyi")
        assert.are.equal("/p/.venv/lib/python3.12/site-packages/numpy", env.stubs)
        assert.are.equal("/p/.venv/bin/python3", env.python)
    end)

    it("keeps a -stubs directory as it is", function()
        assert.are.equal("/p/lib/python3.12/site-packages/rdkit-stubs",
            stub_patches.env_of("/p/lib/python3.12/site-packages/rdkit-stubs/Chem/__init__.pyi").stubs)
    end)

    it("prefers bin/python and falls back to bin/python3", function()
        local root = vim.fn.tempname()
        local path = root .. "/lib/python3.13/site-packages/pkg/__init__.pyi"
        write(root .. "/bin/python3", "")
        assert.are.equal(root .. "/bin/python3", stub_patches.env_of(path).python)
        write(root .. "/bin/python", "")
        assert.are.equal(root .. "/bin/python", stub_patches.env_of(path).python)
        vim.fs.rm(root, { recursive = true })
    end)

    it("returns nil outside a site-packages layout", function()
        assert.is_nil(stub_patches.env_of("/home/me/proj/src/rdkit.py"))
    end)
end)

--- Replace a table field for the duration of a function.
local function with(tbl, key, value, fn)
    local original = tbl[key]
    tbl[key] = value
    local ok, err = pcall(fn)
    tbl[key] = original
    assert(ok, err)
end

--- A candidate environment on disk: a python and a tree with an `__init__.pyi`.
--- @return StubPatch.Env
local function fake_env()
    local root = vim.fn.tempname()
    local env = { stubs = root .. "/site/pkg", python = root .. "/bin/python" }
    write(env.stubs .. "/__init__.pyi", "")
    write(env.python, "")
    return env
end

describe("stub_patches spawn gates", function()
    --- The state `consider` leaves, and how many runs it started. Each run is
    --- then finished, so the next test finds none in progress.
    local function consider(env)
        local runs = {}
        local state
        with(vim, "system", function(_, _, on_exit) table.insert(runs, on_exit) end, function()
            with(vim, "notify", function() end, function()
                stub_patches.consider(env)
                state = stub_patches.state[env.stubs]
                for _, on_exit in ipairs(runs) do on_exit { code = 3 } end
                vim.wait(1000, function() return stub_patches.state[env.stubs] ~= "patching" end, 10)
            end)
        end)
        return state, #runs
    end

    it("skips a tree without an __init__.pyi", function()
        local env = fake_env()
        vim.fs.rm(env.stubs .. "/__init__.pyi")
        assert.are.same({ "ok", 0 }, { consider(env) })
    end)

    it("fails without an interpreter", function()
        local env = fake_env()
        vim.fs.rm(env.python)
        assert.are.same({ "failed", 0 }, { consider(env) })
    end)

    it("skips an unwritable tree", function()
        local env = fake_env()
        vim.uv.fs_chmod(env.stubs, tonumber("555", 8))
        local state, runs = consider(env)
        vim.uv.fs_chmod(env.stubs, tonumber("755", 8))
        assert.are.same({ "ok", 0 }, { state, runs })
    end)

    it("fails without uv on PATH", function()
        local env = fake_env()
        with(vim.fn, "executable", function() return 0 end, function()
            assert.are.same({ "failed", 0 }, { consider(env) })
        end)
    end)

    it("queues a run otherwise", function()
        assert.are.same({ "patching", 1 }, { consider(fake_env()) })
    end)
end)

describe("stub_patches.enqueue", function()
    it("runs one at a time and maps each exit status to a state", function()
        --- @type { cmd: string[], opts: table, on_exit: fun(out: table) }[]
        local runs = {}
        local envs = { fake_env(), fake_env(), fake_env() }
        local messages = {}
        with(vim, "system", function(cmd, opts, on_exit)
            table.insert(runs, { cmd = cmd, opts = opts, on_exit = on_exit })
        end, function()
            with(vim, "notify", function(msg, level) table.insert(messages, { msg, level }) end,
                function()
                    for _, env in ipairs(envs) do stub_patches.enqueue(env) end
                    assert.are.equal(1, #runs)
                    assert.are.equal(envs[1].stubs, runs[1].cmd[#runs[1].cmd])
                    assert.are.equal(vim.fn.stdpath("run"), runs[1].opts.cwd)
                    assert.are.equal("queued", stub_patches.state[envs[2].stubs])

                    for i, code in ipairs { 0, 3, 1 } do
                        runs[i].on_exit { code = code, stdout = "pkg.mod (boom)\n", stderr = "why" }
                        assert.is_true(vim.wait(1000, function()
                            return #runs == math.min(i + 1, 3)
                                and stub_patches.state[envs[i].stubs] ~= "patching"
                        end, 10))
                    end
                end)
        end)
        assert.are.same({ "patched", "ok", "failed" }, vim.tbl_map(function(env)
            return stub_patches.state[env.stubs]
        end, envs))
        assert.are.equal(2, #messages)
        assert.is_truthy(messages[1][1]:find("pkg.mod (boom)", 1, true))
        assert.are.equal(vim.log.levels.ERROR, messages[2][2])
        assert.is_truthy(messages[2][1]:find("why", 1, true))
    end)
end)

local MAIN = [[
import fakepkg
fakepkg.sqrt(2.0)
]]

--- Response of one server to a position request. Other clients may attach to
--- the fixture buffer too, so requests are scoped to the spec's client.
--- @param client vim.lsp.Client
--- @param buf number
--- @param method string
--- @param position lsp.Position
--- @return table|nil result
local function request(client, buf, method, position)
    local response = client:request_sync(method, {
        textDocument = vim.lsp.util.make_text_document_params(buf),
        position = position,
    }, 10000, buf)
    return response and response.result
end

--- Hover text at `fakepkg.sqrt` in the fixture, nil until analysed.
--- @param client vim.lsp.Client
--- @param buf number
--- @return string|nil
local function sqrt_hover(client, buf)
    local result = request(client, buf, "textDocument/hover", { line = 1, character = 9 })
    return result and result.contents.value
end

--- Env discovery, in-place patch and re-analysis against a running server.
--- @param buf number buffer of the fixture's `main.py`
--- @param client vim.lsp.Client
--- @param site string site-packages of the fixture env
local function assert_documents(buf, client, site)
    -- Waiting for the signature is also the wait for initial analysis.
    assert.is_true(vim.wait(20000, function()
        local text = sqrt_hover(client, buf)
        return text ~= nil and text:find("-> float", 1, true) ~= nil
    end, 200))
    assert.is_nil(sqrt_hover(client, buf):find("square root", 1, true))

    local position = assert(stub_patches.imports(buf).fakepkg)
    local result = assert(request(client, buf, "textDocument/declaration", position))
    local location = vim.islist(result) and result[1] or result
    local env = assert(stub_patches.env_of(vim.uri_to_fname(location.uri or location.targetUri)))
    assert.are.equal(vim.uv.fs_realpath(site .. "/fakepkg-stubs"), vim.uv.fs_realpath(env.stubs))

    stub_patches.consider(env)
    assert.is_true(vim.wait(60000, function()
        return stub_patches.state[env.stubs] ~= "patching"
    end, 200))
    assert.are.equal("patched", stub_patches.state[env.stubs])
    assert.is_true(vim.wait(20000, function()
        local text = sqrt_hover(client, buf)
        return text ~= nil and text:find("square root", 1, true) ~= nil
    end, 200))
end

describe("stub_patches", function()
    local server = vim.fn.exepath("basedpyright-langserver")

    it("documents a stub after patching it in its environment", function()
        if server == "" or vim.fn.executable("uv") == 0 then
            pending("basedpyright-langserver or uv not installed")
            return
        end
        -- The config's own handler would patch the ambient environment too.
        pcall(vim.api.nvim_del_augroup_by_name, "my.stub_patches")
        local root = vim.fn.tempname()
        local venv = vim.system({ "uv", "venv", "-q", root .. "/env" }, { text = true }):wait()
        assert(venv.code == 0, venv.stderr)
        local site = assert(vim.fs.find("site-packages",
            { path = root .. "/env/lib", type = "directory" })[1])
        write(site .. "/fakepkg/__init__.py", "from math import sqrt")
        write(site .. "/fakepkg-stubs/__init__.pyi", "def sqrt(x: float) -> float: ...")
        write(root .. "/proj/main.py", MAIN)
        local buf = vim.fn.bufadd(root .. "/proj/main.py")
        vim.fn.bufload(buf)
        local client = assert(vim.lsp.get_client_by_id(assert(vim.lsp.start({
            name = "basedpyright",
            cmd = { server, "--stdio" },
            root_dir = root .. "/proj",
            settings = { python = { pythonPath = root .. "/env/bin/python" } },
        }, { bufnr = buf, reuse_client = function() return false end }))))

        local ok, err = pcall(assert_documents, buf, client, site)
        client:stop(true)
        vim.fs.rm(root, { recursive = true, force = true })
        assert(ok, err)
    end)
end)
