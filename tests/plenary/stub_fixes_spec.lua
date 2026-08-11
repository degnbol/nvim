---@diagnostic disable: undefined-global
local stub_fixes = require "autocmds/stub_fixes"

--- The configured fix for a package, as the autocmd would use it.
--- @param pkg string
--- @return StubFix
local function fix_for(pkg)
    for _, fix in ipairs(stub_fixes.fixes) do
        if fix.pkg == pkg then return fix end
    end
    error("no fix configured for " .. pkg)
end

describe("stub_fixes.import_position", function()
    local function position(lines)
        local buf = vim.api.nvim_create_buf(false, true)
        vim.api.nvim_buf_set_lines(buf, 0, -1, false, lines)
        return stub_fixes.import_position(buf, "rdkit")
    end

    it("finds the token in a plain import", function()
        assert.are.same({ line = 0, character = 7 }, position { "import rdkit" })
    end)

    it("finds the token in a from-import below other code", function()
        assert.are.same({ line = 1, character = 5 },
            position { "import os", "from rdkit import Chem" })
    end)

    it("finds an indented import", function()
        assert.are.same({ line = 1, character = 11 },
            position { "try:", "    import rdkit.Chem" })
    end)

    it("ignores a package that merely starts with the name", function()
        assert.is_nil(position { "import rdkit_extras" })
    end)

    it("ignores buffers without the import", function()
        assert.is_nil(position { "import numpy as np", "# rdkit is unused" })
    end)
end)

describe("stub_fixes.env_of", function()
    it("derives the stub dir and interpreter from a site-packages path", function()
        local env = stub_fixes.env_of(
            "/p/.venv/lib/python3.12/site-packages/rdkit/__init__.py", "rdkit")
        assert.are.equal("/p/.venv/lib/python3.12/site-packages/rdkit-stubs", env.stubs)
        assert.are.equal("/p/.venv/bin/python", env.python)
    end)

    it("returns nil outside a site-packages layout", function()
        assert.is_nil(stub_fixes.env_of("/home/me/proj/src/rdkit.py", "rdkit"))
    end)
end)

describe("stub_fixes.is_unpatched", function()
    local defects = fix_for("rdkit").defects

    it("detects a C++ signature block", function()
        assert.is_true(stub_fixes.is_unpatched('    """\n        C++ signature :\n    """\n',
            defects))
    end)

    it("detects an untyped property", function()
        assert.is_true(stub_fixes.is_unpatched(
            "    @property\n    def allowCXSMILES(*args, **kwargs):\n        ...\n", defects))
    end)

    it("accepts a patched stub with untyped plain methods", function()
        assert.is_false(stub_fixes.is_unpatched(
            "    @property\n    def allowCXSMILES(self) -> bool: ...\n" ..
            "    @staticmethod\n    def flush(*args, **kwargs): ...\n", defects))
    end)
end)

-- Mirrors the real defect: an untyped property in an installed `rdkit-stubs/`,
-- fixed in place and re-typed by the running server without a restart.
local STUB = [[
class SmilesParserParams:
    """
    Parameters controlling SMILES Parsing
    """
    def __init__(self) -> None:
        """
            C++ signature :
                void __init__(_object*)
        """
    @property
    def allowCXSMILES(*args, **kwargs):
        """
        controls whether or not the CXSMILES extensions are parsed
        """
    @allowCXSMILES.setter
    def allowCXSMILES(*args, **kwargs):
        ...
]]

local SOURCE = [[
class SmilesParserParams:
    @property
    def allowCXSMILES(self):
        return True
]]

local MAIN = [[
from rdkit.Chem.rdmolfiles import SmilesParserParams
params = SmilesParserParams()
params.allowCXSMILES
]]

--- Write text to a path, creating parent directories.
local function write(path, text)
    vim.fn.mkdir(vim.fs.dirname(path), "p")
    vim.fn.writefile(vim.split(text, "\n"), path)
end

--- Build an environment holding a fake rdkit and its unpatched stubs.
--- @param root string directory to build under
--- @return string site_packages
local function fake_env(root)
    local site = root .. "/env/lib/python3.13/site-packages"
    write(site .. "/rdkit/__init__.py", "")
    write(site .. "/rdkit/Chem/__init__.py", "")
    write(site .. "/rdkit/Chem/rdmolfiles.py", SOURCE)
    write(site .. "/rdkit-stubs/__init__.pyi", "")
    write(site .. "/rdkit-stubs/Chem/__init__.pyi", "")
    write(site .. "/rdkit-stubs/Chem/rdmolfiles.pyi", STUB)
    -- The stubs are patched by the env's own interpreter, which must import the
    -- package to read property types off a live instance.
    local python = root .. "/env/bin/python"
    write(python, ("#!/bin/sh\nPYTHONPATH=%s exec %s \"$@\"\n")
        :format(site, vim.fn.exepath("python3")))
    vim.uv.fs_chmod(python, tonumber("755", 8))
    return site
end

--- Response of one server to a position request. Other clients attach to the
--- fixture buffer too (the config enables basedpyright for python), and they
--- resolve rdkit elsewhere, so requests must be scoped to the spec's client.
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

--- Hover text at `params.allowCXSMILES` in the fixture, nil until analysed.
--- @param client vim.lsp.Client
--- @param buf number
--- @return string|nil
local function property_hover(client, buf)
    local result = request(client, buf, "textDocument/hover", { line = 2, character = 8 })
    return result and result.contents.value
end

--- Whether the hovered property is reported with a given type.
--- @param client vim.lsp.Client
--- @param buf number
--- @param type_name string
--- @return boolean
local function hover_shows(client, buf, type_name)
    local text = property_hover(client, buf)
    return text ~= nil and text:find(type_name, 1, true) ~= nil
end

--- Env discovery, in-place patch and re-analysis against a running server.
--- @param fix StubFix
--- @param buf number buffer of the fixture's `main.py`
--- @param client vim.lsp.Client
--- @param site string site-packages of the fixture env
local function assert_retypes(fix, buf, client, site)
    -- Waiting for the untyped form is also the wait for initial analysis.
    assert.is_true(vim.wait(20000, function() return hover_shows(client, buf, "*args") end, 200))

    local result = assert(request(client, buf, "textDocument/definition",
        assert(stub_fixes.import_position(buf, fix.pkg))))
    local location = vim.islist(result) and result[1] or result
    local env = assert(stub_fixes.env_of(vim.uri_to_fname(location.uri or location.targetUri),
        fix.pkg))
    assert.are.equal(site .. "/rdkit-stubs", env.stubs)
    assert.is_true(stub_fixes.is_unpatched(
        assert(require("utils/init").readtext(env.stubs .. "/" .. fix.probe)), fix.defects))

    stub_fixes.patch(fix, env, client)
    assert.is_true(vim.wait(20000, function() return hover_shows(client, buf, "bool") end, 200))
end

describe("stub_fixes", function()
    local server = vim.fn.exepath("basedpyright-langserver")

    it("re-types a property after patching the installed stubs", function()
        if server == "" then
            pending("basedpyright-langserver not installed")
            return
        end
        -- The config's own handler would prompt about the ambient python env,
        -- which no headless test can answer. The test drives the patch itself.
        vim.api.nvim_del_augroup_by_name("my.stub_fixes")
        local fix = fix_for("rdkit")
        assert.is_truthy(vim.uv.fs_stat(fix.script))
        local root = vim.fn.tempname()
        local site = fake_env(root)
        write(root .. "/proj/main.py", MAIN)
        local buf = vim.fn.bufadd(root .. "/proj/main.py")
        vim.fn.bufload(buf)
        local client = assert(vim.lsp.get_client_by_id(assert(vim.lsp.start({
            name = "stub_fixes_spec",
            cmd = { server, "--stdio" },
            root_dir = root .. "/proj",
            settings = { basedpyright = { analysis = { extraPaths = { site } } } },
        }, { bufnr = buf }))))

        local ok, err = pcall(assert_retypes, fix, buf, client, site)
        client:stop(true)
        vim.fs.rm(root, { recursive = true, force = true })
        assert(ok, err)
    end)
end)
