---@diagnostic disable: undefined-global
local uv_script = require "autocmds/uv_script_env"

--- Whether a buffer of these lines is seen as carrying script metadata.
--- @param lines string[]
--- @return boolean
local function detects(lines)
    local buf = vim.api.nvim_create_buf(false, true)
    vim.api.nvim_buf_set_lines(buf, 0, -1, false, lines)
    return uv_script.has_script_block(buf)
end

describe("uv_script_env.has_script_block", function()
    it("finds a block below a shebang", function()
        assert.is_true(detects {
            "#!/usr/bin/env -S uv run --script",
            "# /// script",
            '# dependencies = ["rdkit"]',
            "# ///",
            "import rdkit",
        })
    end)

    it("accepts a comment line with nothing after the hash", function()
        assert.is_true(detects { "# /// script", "#", "# ///" })
    end)

    it("finds a block inside a docstring, as uv does", function()
        assert.is_true(detects { '"""', "# /// script", "# ///", '"""' })
    end)

    it("ignores a python file without metadata", function()
        assert.is_false(detects { "import os", "# a comment", "print(os.getcwd())" })
    end)

    it("ignores an unterminated block", function()
        assert.is_false(detects { "# /// script", '# dependencies = ["rdkit"]' })
    end)

    it("ignores a block interrupted by code", function()
        assert.is_false(detects { "# /// script", "import os", "# ///" })
    end)

    it("ignores a block of another type", function()
        assert.is_false(detects { "# /// pyproject", "# ///" })
    end)
end)

--- A config as `vim.lsp.start` would receive it.
--- @param root string|nil
--- @param python string|nil environment a script's client is pinned to
--- @return table config
local function config(root, python)
    return { name = "basedpyright", root_dir = root, uv_script_python = python }
end

--- A stand-in for a running client, holding the fields `same_env` reads.
--- @param root string|nil
--- @param python string|nil
--- @return table client
local function client(root, python)
    return {
        name = "basedpyright",
        config = config(root, python),
        workspace_folders = root and vim.lsp._get_workspace_folders(root) or nil,
        is_stopped = function() return false end,
    }
end

describe("uv_script_env.same_env", function()
    it("reuses a client of the same root and environment", function()
        assert.is_true(uv_script.same_env(client("/p", "/env/bin/python"),
            config("/p", "/env/bin/python")))
    end)

    it("keeps clients of one root but different environments apart", function()
        assert.is_false(uv_script.same_env(client("/p", "/a/bin/python"),
            config("/p", "/b/bin/python")))
    end)

    it("keeps a script client from serving the project's own root", function()
        assert.is_false(uv_script.same_env(client("/p", "/env/bin/python"), config("/p")))
        assert.is_false(uv_script.same_env(client("/p"), config("/p", "/env/bin/python")))
    end)

    it("rejects a client of another root", function()
        assert.is_false(uv_script.same_env(client("/other", "/env/bin/python"),
            config("/p", "/env/bin/python")))
    end)

    it("rejects a client of another server", function()
        local other = client("/p", "/env/bin/python")
        other.name = "ruff"
        assert.is_false(uv_script.same_env(other, config("/p", "/env/bin/python")))
    end)

    it("reuses a rootless client only for a rootless config", function()
        assert.is_true(uv_script.same_env(client(nil, nil), config(nil, nil)))
        assert.is_false(uv_script.same_env(client("/p", nil), config(nil, nil)))
    end)

    it("still reuses a project client whose interpreter was switched", function()
        -- `LspPyrightSetPythonPath` rewrites the live settings, which is why the
        -- pin, not `settings.python.pythonPath`, decides.
        local switched = client("/p")
        switched.settings = { python = { pythonPath = "/conda/bin/python" } }
        assert.is_true(uv_script.same_env(switched, config("/p")))
    end)
end)

describe("uv_script_env.resolve", function()
    --- Result of resolving a script, waiting out the two uv calls.
    --- @param path string
    --- @return string|nil python
    --- @return string|nil err
    local function resolve(path)
        local python, err, done
        uv_script.resolve(path, function(p, e)
            python, err, done = p, e, true
        end)
        assert.is_true(vim.wait(60000, function() return done end, 50))
        return python, err
    end

    it("reports a script that is not on disk", function()
        local python, err = resolve(vim.fn.tempname() .. ".py")
        assert.is_nil(python)
        assert.is_truthy(assert(err):find("not on disk", 1, true))
    end)

    it("answers with the script's own environment", function()
        if vim.fn.executable("uv") == 0 then
            pending("uv not installed")
            return
        end
        local path = vim.fn.tempname() .. ".py"
        vim.fn.writefile({ "# /// script", '# requires-python = ">=3.11"',
            "# dependencies = []", "# ///" }, path)
        local python, err = resolve(path)
        assert.is_nil(err)
        python = assert(python)
        -- The point of the sync: without it uv answers with a system interpreter.
        assert.is_truthy(python:find("/environments-v2/", 1, true))
        assert.is_truthy(vim.uv.fs_stat(python))
        vim.fs.rm(vim.fs.dirname(vim.fs.dirname(python)), { recursive = true })
        vim.fs.rm(path, {})
    end)
end)

--- The uv environment a client resolves imports in.
--- @param client vim.lsp.Client
--- @return string|nil
local function env_of(client)
    return client.config.uv_script_python
end

--- Clients of the server under test attached to a buffer.
--- @param bufnr integer
--- @return vim.lsp.Client[]
local function clients(bufnr)
    return vim.lsp.get_clients({ bufnr = bufnr, name = "basedpyright" })
end

--- Buffer of a dependency-free PEP 723 script, written under a directory.
--- @param dir string
--- @param name string
--- @return integer bufnr
local function script_buf(dir, name)
    local path = dir .. "/" .. name
    vim.fn.writefile({ "# /// script", '# requires-python = ">=3.11"',
        "# dependencies = []", "# ///", "import os" }, path)
    local buf = vim.fn.bufadd(path)
    vim.fn.bufload(buf)
    return buf
end

--- Run the `root_dir` hook for a buffer and wait until it settles on one client.
--- @param bufnr integer
--- @return vim.lsp.Client
local function route(bufnr)
    local fell_back = false
    uv_script.root_dir(bufnr, function() fell_back = true end)
    assert.is_true(vim.wait(60000, function()
        return fell_back or #clients(bufnr) == 1
    end, 100))
    assert.is_false(fell_back)
    return clients(bufnr)[1]
end

--- One client per script, each in its own environment, each stopped with its
--- last buffer. Two scripts of one directory share a root, so only the pin keeps
--- them apart.
--- @param dir string to write the scripts under
--- @param envs string[] collects each environment created, for the caller to remove
local function assert_client_per_script(dir, envs)
    local one, two = script_buf(dir, "one.py"), script_buf(dir, "two.py")
    local c1, c2 = route(one), route(two)
    for _, c in ipairs({ c1, c2 }) do envs[#envs + 1] = env_of(c) end
    -- uv keys an environment by script path, so identical metadata still yields
    -- two of them.
    assert.are_not.equal(env_of(c1), env_of(c2))
    assert.are.same({ one }, vim.tbl_keys(c1.attached_buffers))
    assert.are.same({ two }, vim.tbl_keys(c2.attached_buffers))

    -- A client of the wrong environment goes when the hook runs again.
    assert.is_true(vim.lsp.buf_attach_client(one, c2.id))
    assert.are.equal(c1.id, route(one).id)
    assert.is_false(c2:is_stopped())

    vim.api.nvim_buf_delete(one, { force = true })
    vim.api.nvim_buf_delete(two, { force = true })
    assert.is_true(vim.wait(10000, function()
        return c1:is_stopped() and c2:is_stopped()
    end, 100))
end

describe("uv_script_env clients", function()
    it("serves each script from its own environment", function()
        local server = vim.fn.exepath("basedpyright-langserver")
        if server == "" or vim.fn.executable("uv") == 0 then
            pending("basedpyright-langserver or uv not installed")
            return
        end
        -- nvim-lspconfig, which supplies `cmd`, is not on the test runtimepath.
        vim.lsp.config("basedpyright", { cmd = { server, "--stdio" } })
        local dir = vim.fn.tempname()
        vim.fn.mkdir(dir, "p")

        local envs = {}
        local ok, err = pcall(assert_client_per_script, dir, envs)
        for _, client in ipairs(vim.lsp.get_clients({ name = "basedpyright" })) do
            client:stop(true)
        end
        for _, env in ipairs(envs) do
            vim.fs.rm(vim.fs.dirname(vim.fs.dirname(env)), { recursive = true })
        end
        vim.fs.rm(dir, { recursive = true })
        assert(ok, err)
    end)
end)
