--- Fix defective type stubs where a python package installs them.
---
--- Some binary packages ship a `<package>-stubs/` tree inside their own wheel
--- with stub-generator defects. `lsp_ext/python_stubs/` holds fixed copies, but
--- basedpyright discards the language server's `stubPath` for any project owning
--- a `pyrightconfig.json`/`[tool.basedpyright]`, so those projects only ever
--- read the installed stubs. Rewrite those in place, then nudge the server to
--- re-analyse the rewritten files.
---
--- Cover another package by appending to `M.fixes`; nothing else here is
--- package-specific.
local util = require "utils/init"

local M = {}

--- @class StubFix
--- @field pkg string import name; gates the buffer scan and names `<package>-stubs`
--- @field probe string path inside the stub tree whose head carries the marker
--- @field marker string comment the script writes at the head of each file it rewrote
--- @field script string fixer, run as `<env python> <script> <package> --in-place`

--- @class StubFix.Env
--- @field stubs string the installed `<package>-stubs` directory
--- @field python string interpreter of the environment holding it

--- @type StubFix[]
M.fixes = {
    {
        pkg = "rdkit",
        -- Boost.Python defeats pybind11-stubgen: it leaves the raw C++ signature
        -- in docstrings and types every property `(*args, **kwargs)`. This module
        -- holds SmilesParserParams, which shows both.
        probe = "Chem/rdmolfiles.pyi",
        marker = "# fix_pybind_stubs:",
        script = vim.fn.stdpath("config") .. "/lsp_ext/python_stubs/fix_pybind_stubs.py",
    },
}

--- Last decision per stub directory, so a session prompts at most once for each.
--- The values (`ok`/`prompting`/`patching`/`patched`/`declined`/`failed`) are for
--- inspection; only presence is acted on.
--- @type table<string, string>
local state = {}

--- LSP position of a package's name in the first line importing it.
--- @param bufnr number
--- @param pkg string
--- @return lsp.Position|nil
function M.import_position(bufnr, pkg)
    for row, line in ipairs(vim.api.nvim_buf_get_lines(bufnr, 0, -1, false)) do
        -- `%f[^%w_]`, not `%f[%W]`: Lua's `%w` excludes `_`, so the latter would
        -- accept `<package>_extras`.
        local col = line:match("^%s*import ()" .. pkg .. "%f[^%w_]")
            or line:match("^%s*from ()" .. pkg .. "%f[^%w_]")
        if col then return { line = row - 1, character = col - 1 } end
    end
end

--- Stub directory and interpreter of the environment a module path belongs to.
--- @param path string file under `<env>/lib/python*/site-packages/`
--- @param pkg string
--- @return StubFix.Env|nil env nil if the path is not that layout
function M.env_of(path, pkg)
    local root, site = path:match("^(.*)(/lib/python[^/]*/site%-packages)/")
    if not root then return nil end
    return { stubs = root .. site .. "/" .. pkg .. "-stubs", python = root .. "/bin/python" }
end

--- Tell a server that every stub file under a directory changed on disk.
--- @param client vim.lsp.Client
--- @param stubs string stub directory
function M.notify_changed(client, stubs)
    local changes = {}
    for _, path in ipairs(vim.fs.find(
        function(name) return name:match("%.pyi$") end,
        { path = stubs, type = "file", limit = math.huge })
    ) do
        table.insert(changes, {
            uri = vim.uri_from_fname(path),
            type = vim.lsp.protocol.FileChangeType.Changed,
        })
    end
    client:notify("workspace/didChangeWatchedFiles", { changes = changes })
end

--- Rewrite an environment's stubs, then have the server re-analyse them.
--- @param fix StubFix
--- @param env StubFix.Env
--- @param client vim.lsp.Client
function M.patch(fix, env, client)
    state[env.stubs] = "patching"
    vim.system({ env.python, fix.script, fix.pkg, "--in-place" }, { text = true },
        function(out)
            vim.schedule(function()
                if out.code ~= 0 then
                    state[env.stubs] = "failed"
                    vim.notify(("Patching %s failed (exit %d):\n%s")
                        :format(env.stubs, out.code, util.strip(out.stderr or "")),
                        vim.log.levels.ERROR)
                    return
                end
                state[env.stubs] = "patched"
                if not client:is_stopped() then M.notify_changed(client, env.stubs) end
                local untyped = util.strip(out.stdout or "")
                vim.notify("Patched " .. env.stubs .. (untyped ~= "" and "\n" .. untyped or ""))
            end)
        end)
end

--- Offer to patch an environment's stubs unless the fixer's marker is already there.
--- @param fix StubFix
--- @param env StubFix.Env
--- @param client vim.lsp.Client
local function consider(fix, env, client)
    local probe = util.readtext(env.stubs .. "/" .. fix.probe)
    -- No stubs, or the fixer's own record of having rewritten them, is nothing to
    -- do. The marker rather than leftover defects, because the fixer legitimately
    -- leaves some behind (properties it cannot type); it still self-heals when a
    -- reinstall restores the bad stubs, which drops the marker with them.
    if not probe or vim.startswith(probe, fix.marker) then
        state[env.stubs] = "ok"
        return
    end
    if not vim.uv.fs_stat(env.python) then
        state[env.stubs] = "failed"
        vim.notify(("No interpreter at %s, cannot patch %s"):format(env.python, env.stubs),
            vim.log.levels.WARN)
        return
    end
    -- Set before the blocking prompt so a second buffer's callback bails.
    state[env.stubs] = "prompting"
    vim.schedule(function()
        if vim.fn.confirm("Patch " .. env.stubs .. "?", "&Yes\n&No") == 1 then
            M.patch(fix, env, client)
        else
            state[env.stubs] = "declined"
        end
    end)
end

--- Path of the first location in a definition response.
--- @param result lsp.Location|lsp.Location[]|lsp.LocationLink[]|nil
--- @return string|nil
local function definition_path(result)
    local loc = result and (vim.islist(result) and result[1] or result)
    local uri = loc and (loc.uri or loc.targetUri)
    return uri and vim.uri_to_fname(uri)
end

--- Check a package's stubs in whichever environment the server resolves it from.
--- basedpyright does not expose its resolved interpreter, but the definition it
--- answers with is inside the very site-packages whose stubs need fixing.
--- @param fix StubFix
--- @param bufnr number
--- @param position lsp.Position of the package name in an import
--- @param client vim.lsp.Client
local function check_package(fix, bufnr, position, client)
    client:request("textDocument/definition", {
        textDocument = vim.lsp.util.make_text_document_params(bufnr),
        position = position,
    }, function(err, result)
        if err then
            return vim.notify(("Locating %s failed: %s"):format(fix.pkg, err.message),
                vim.log.levels.WARN)
        end
        local path = definition_path(result)
        local env = path and M.env_of(path, fix.pkg)
        if env and not state[env.stubs] then consider(fix, env, client) end
    end, bufnr)
end

vim.api.nvim_create_autocmd("LspAttach", {
    group = vim.api.nvim_create_augroup("my.stub_fixes", { clear = true }),
    callback = function(args)
        local client = vim.lsp.get_client_by_id(args.data.client_id)
        if not client or client.name ~= "basedpyright" then return end
        for _, fix in ipairs(M.fixes) do
            local position = M.import_position(args.buf, fix.pkg)
            if position then check_package(fix, args.buf, position, client) end
        end
    end,
})

return M
