--- Patch the stubs of every package a buffer imports, in the environment the
--- language server resolves each from.
---
--- basedpyright discards the language server's `stubPath` for any project owning
--- a `pyrightconfig.json`/`[tool.basedpyright]`, so those projects only ever read
--- an environment's own stubs. `lsp_ext/stub_patches/patch_stubs.py` rewrites
--- those in place and decides everything about a tree: which package it
--- describes, whether it applies, whether it is current, which repair runs. This
--- module only finds candidate trees, checks a run can start, queues runs one at
--- a time, and nudges the servers to re-analyse what changed.
local M = {}

local SCRIPTS = vim.fn.stdpath("config") .. "/lsp_ext/stub_patches"
--- A backstop behind the script's own watchdog, `patch_stubs.WATCHDOG_S`.
local TIMEOUT_MS = 660 * 1000
--- `patch_stubs.NOTHING_WRITTEN`
local NOTHING_WRITTEN = 3

--- @class StubPatch.Env
--- @field stubs string `site-packages/<first component>` of a declaration
--- @field python string interpreter of the environment holding `stubs`

--- State per stub directory for the session. Set before any async step, so a
--- second buffer's callback bails.
--- @type table<string, "queued"|"patching"|"patched"|"ok"|"failed">
M.state = {}

--- @type StubPatch.Env[]
local queue = {}
local running = false
local warned_no_uv = false

local IMPORTS = [[
(import_statement name: (dotted_name . (identifier) @package))
(import_statement name: (aliased_import name: (dotted_name . (identifier) @package)))
(import_from_statement module_name: (dotted_name . (identifier) @package))
]]

--- Top-level packages a Python buffer's absolute imports name.
--- @param bufnr number
--- @return table<string, lsp.Position> packages top-level package -> position of its name
function M.imports(bufnr)
    local parser = vim.treesitter.get_parser(bufnr, "python", { error = false })
    if not parser then return {} end
    local query = vim.treesitter.query.parse("python", IMPORTS)
    local packages = {}
    for _, node in query:iter_captures(parser:parse()[1]:root(), bufnr) do
        local name = vim.treesitter.get_node_text(node, bufnr)
        if not packages[name] then
            local row, col = node:range()
            packages[name] = { line = row, character = col }
        end
    end
    return packages
end

--- The environment whose site-packages holds a declaration.
--- @param path string declaration path
--- @return StubPatch.Env|nil env nil outside `<env>/lib/python*/site-packages/`
function M.env_of(path)
    local root, site, first = path:match("^(.*)(/lib/python[^/]*/site%-packages)/([^/]+)")
    if not root then return nil end
    local python = root .. "/bin/python"
    if not vim.uv.fs_stat(python) then python = root .. "/bin/python3" end
    return { stubs = root .. site .. "/" .. first, python = python }
end

--- Tell every running basedpyright client that the `.pyi` files under a
--- directory changed on disk.
--- @param stubs string
function M.notify_changed(stubs)
    local changes = {}
    for _, path in ipairs(vim.fs.find(
        function(name) return vim.endswith(name, ".pyi") end,
        { path = stubs, type = "file", limit = math.huge })
    ) do
        table.insert(changes, {
            uri = vim.uri_from_fname(path),
            type = vim.lsp.protocol.FileChangeType.Changed,
        })
    end
    for _, client in ipairs(vim.lsp.get_clients { name = "basedpyright" }) do
        client:notify("workspace/didChangeWatchedFiles", { changes = changes })
    end
end

--- Record how a run ended and report it.
--- @param env StubPatch.Env
--- @param out vim.SystemCompleted
local function finish(env, out)
    if out.code == 0 then
        M.state[env.stubs] = "patched"
        M.notify_changed(env.stubs)
        local skipped = vim.trim(out.stdout or "")
        vim.notify("Patched " .. env.stubs .. (skipped ~= "" and "\nskipped:\n" .. skipped or ""))
    elseif out.code == NOTHING_WRITTEN then
        M.state[env.stubs] = "ok"
    else
        M.state[env.stubs] = "failed"
        vim.notify(("Patching %s failed (exit %d):\n%s")
            :format(env.stubs, out.code, vim.trim(out.stderr or "")), vim.log.levels.ERROR)
    end
end

--- Start the next queued run, if any.
local function run_next()
    local env = table.remove(queue, 1)
    running = env ~= nil
    if not env then return end
    M.state[env.stubs] = "patching"
    vim.system({
        "uv", "run", "--no-project", "--python", env.python,
        "--with-requirements", SCRIPTS .. "/requirements.txt",
        SCRIPTS .. "/patch_stubs.py", env.stubs,
    }, {
        text = true,
        -- An import writing relative files writes here, not into the project.
        cwd = vim.fn.stdpath("run"),
        timeout = TIMEOUT_MS,
    }, vim.schedule_wrap(function(out)
        finish(env, out)
        run_next()
    end))
end

--- Queue a run; start it now if none is running.
--- @param env StubPatch.Env
function M.enqueue(env)
    M.state[env.stubs] = "queued"
    table.insert(queue, env)
    if not running then run_next() end
end

--- Queue a run, unless the tree has no `__init__.pyi` or is unwritable (state
--- `ok`), or the interpreter or uv is missing (state `failed`, with a warning).
--- @param env StubPatch.Env
function M.consider(env)
    if not vim.uv.fs_stat(env.stubs .. "/__init__.pyi") or not vim.uv.fs_access(env.stubs, "W") then
        M.state[env.stubs] = "ok"
    elseif not vim.uv.fs_stat(env.python) then
        M.state[env.stubs] = "failed"
        vim.notify(("No interpreter at %s, cannot patch %s"):format(env.python, env.stubs),
            vim.log.levels.WARN)
    elseif vim.fn.executable("uv") == 0 then
        M.state[env.stubs] = "failed"
        if not warned_no_uv then
            warned_no_uv = true
            vim.notify("No uv on PATH, cannot patch stubs", vim.log.levels.WARN)
        end
    else
        M.enqueue(env)
    end
end

--- Path of the first location in a declaration response.
--- @param result lsp.Location|lsp.Location[]|lsp.LocationLink[]|nil
--- @return string|nil
local function location_path(result)
    local loc = result and (vim.islist(result) and result[1] or result)
    local uri = loc and (loc.uri or loc.targetUri)
    return uri and vim.uri_to_fname(uri)
end

--- Consider a package's stubs in whichever environment the server resolves it
--- from, unless they already have a state.
--- @param package string
--- @param bufnr number
--- @param position lsp.Position of the package name in an import
--- @param client vim.lsp.Client
local function locate(package, bufnr, position, client)
    -- basedpyright does not expose its resolved interpreter, but the declaration
    -- it answers with is inside the very site-packages holding the stubs. A
    -- declaration, not a definition: pyright answers a definition request with
    -- the source a stub shadows where it finds one.
    client:request("textDocument/declaration", {
        textDocument = vim.lsp.util.make_text_document_params(bufnr),
        position = position,
    }, function(err, result)
        if err then
            return vim.notify(("Locating %s failed: %s"):format(package, err.message),
                vim.log.levels.WARN)
        end
        local path = location_path(result)
        local env = path and M.env_of(path)
        if env and not M.state[env.stubs] then M.consider(env) end
    end, bufnr)
end

vim.api.nvim_create_autocmd("LspAttach", {
    group = vim.api.nvim_create_augroup("my.stub_patches", { clear = true }),
    callback = function(args)
        local client = vim.lsp.get_client_by_id(args.data.client_id)
        if not client or client.name ~= "basedpyright" then return end
        for package, position in pairs(M.imports(args.buf)) do
            locate(package, args.buf, position, client)
        end
    end,
})

return M
