--- Give a PEP 723 uv script its own basedpyright, pointed at the script's env.
---
--- A standalone `uv run --script` file resolves its dependencies against an
--- environment of its own under `~/.cache/uv/environments-v2/`, which
--- basedpyright knows nothing about, so every import in it comes up
--- unresolved. `settings.python.pythonPath` fixes that, but it is per-client and
--- one directory can hold several scripts with different environments — hence a
--- client per script, started from basedpyright's `root_dir` hook (see
--- `lsp/basedpyright.lua`) and stopped again when its last buffer detaches.
---
--- The clients keep the name `basedpyright`, so `stub_fixes.lua` covers them
--- too; the `uv_script_python` pin below is what tells them apart from the
--- project's own.
local util = require "utils/init"

local M = {}

local SERVER = "basedpyright"

--- @class UvScriptConfig : vim.lsp.ClientConfig
--- @field uv_script_python? string the environment a script's client is pinned
---   to. Set only here, so it also marks which clients are ours. `settings`
---   cannot serve either purpose: it is the live table, which
---   `LspPyrightSetPythonPath` rewrites in place.

--- The environment a client is pinned to, nil for the project's own client.
--- @param client vim.lsp.Client
--- @return string|nil
local function pin(client)
    return (client.config --[[@as UvScriptConfig]]).uv_script_python
end

--- Whether a buffer carries a PEP 723 `script` metadata block, which is what uv
--- reads to build the script's environment — not the shebang, which can declare
--- dependencies (`uv run --with rdkit`) that get no environment of their own.
--- Follows uv's reading of https://peps.python.org/pep-0723/#specification: by
--- text alone, so a block inside a docstring counts, as it does for uv.
--- @param bufnr integer
--- @return boolean
function M.has_script_block(bufnr)
    local open = false
    for _, line in ipairs(vim.api.nvim_buf_get_lines(bufnr, 0, -1, false)) do
        if open then
            if line == "# ///" then return true end
            open = line == "#" or line:match("^# ") ~= nil
        end
        if not open then open = line == "# /// script" end
    end
    return false
end

--- @param what string
--- @param out vim.SystemCompleted
--- @return string message
local function failure(what, out)
    return ("%s failed (exit %d):\n%s"):format(what, out.code, util.strip(out.stderr or ""))
end

--- Interpreter of the uv environment belonging to a PEP 723 script.
---
--- The sync is not optional: on an environment that has never been materialised
--- `uv python find --script` answers with a *system* interpreter and exit 0 — a
--- wrong answer rather than an error. uv reads the script from disk, so an
--- unsaved buffer cannot resolve.
--- @param path string script file
--- @param cb fun(python: string|nil, err: string|nil) on the main loop
function M.resolve(path, cb)
    cb = vim.schedule_wrap(cb)
    -- vim.system raises on ENOENT rather than reporting it.
    if vim.fn.executable("uv") == 0 then return cb(nil, "uv is not installed") end
    if not vim.uv.fs_stat(path) then
        return cb(nil, ("%s is not on disk; save it to resolve its uv environment")
            :format(path))
    end
    vim.system({ "uv", "sync", "--script", path }, { text = true }, function(sync)
        if sync.code ~= 0 then return cb(nil, failure("uv sync --script", sync)) end
        vim.system({ "uv", "python", "find", "--script", path }, { text = true },
            function(found)
                if found.code ~= 0 then
                    return cb(nil, failure("uv python find --script", found))
                end
                cb(util.strip(found.stdout or ""))
            end)
    end)
end

--- Whether a running client already serves a config, script environment included.
---
--- Mirrors neovim's own `reuse_client_default` — name and workspace folders —
--- and adds the pin, which it never compares: a plain `.py` buffer whose project
--- root happens to be a script's directory would otherwise be handed that
--- script's environment. Two unpinned sides compare equal, so the project's own
--- clients keep the default's behaviour exactly.
--- @param client vim.lsp.Client
--- @param config UvScriptConfig
--- @return boolean
function M.same_env(client, config)
    if client.name ~= config.name or client:is_stopped() then return false end
    if pin(client) ~= config.uv_script_python then return false end
    local wanted = vim.lsp._get_workspace_folders(config.workspace_folders or config.root_dir)
    if not wanted or not next(wanted) then
        -- Reuse only a client that was itself configured without folders.
        local configured = vim.lsp._get_workspace_folders(
            client.config.workspace_folders or client.config.root_dir)
        return not configured or not next(configured)
    end
    local uris = {}
    for _, folder in ipairs(client.workspace_folders or {}) do uris[folder.uri] = true end
    return vim.iter(wanted):all(function(folder) return uris[folder.uri] ~= nil end)
end

--- Detach every basedpyright a buffer should no longer be served by.
---
--- `root_dir` runs again on each `:e`, and `lsp_enable_callback` only detaches a
--- client that can no longer *start* for the buffer — which one of the same name
--- and filetype never is. So an environment that changed under the buffer, or a
--- script block gained or lost, would otherwise leave two clients on it.
--- @param bufnr integer
--- @param python string|nil environment the buffer's own client is pinned to
local function detach_others(bufnr, python)
    for _, client in ipairs(vim.lsp.get_clients({ bufnr = bufnr, name = SERVER })) do
        if pin(client) ~= python then vim.lsp.buf_detach_client(bufnr, client.id) end
    end
end

--- Start a basedpyright for one script, resolving its imports in `python`.
--- @param bufnr integer
--- @param python string interpreter of the script's uv environment
local function start(bufnr, python)
    if not vim.api.nvim_buf_is_loaded(bufnr) then return end
    -- The resolved config, not a `dofile` of `lsp/basedpyright.lua`: that file is
    -- merely one input to the merge and lacks nvim-lspconfig's cmd and on_attach.
    --- @type UvScriptConfig
    local config = vim.deepcopy(assert(vim.lsp.config[SERVER]))
    -- A workspace folder of its own keeps this client apart from the project's
    -- even where `same_env` is bypassed. Nothing is lost by it: basedpyright
    -- searches upward for `pyrightconfig.json`, so a parent's config still applies.
    config.root_dir = vim.fs.dirname(vim.api.nvim_buf_get_name(bufnr))
    config.uv_script_python = python
    config.settings = vim.tbl_deep_extend("force", config.settings or {},
        { python = { pythonPath = python } })
    if vim.lsp.start(config, { bufnr = bufnr, reuse_client = M.same_env }) then
        detach_others(bufnr, python)
    end
end

--- basedpyright `root_dir` hook, see |lsp-root_dir()|.
---
--- Every other Python buffer takes the ordinary path — `on_dir` with no argument
--- — which leaves the root to the config's `root_markers`.
--- @param bufnr integer
--- @param on_dir fun(root_dir: string|nil)
function M.root_dir(bufnr, on_dir)
    if not M.has_script_block(bufnr) then
        detach_others(bufnr, nil)
        return on_dir()
    end
    M.resolve(vim.api.nvim_buf_get_name(bufnr), function(python, err)
        if not python then
            vim.notify(assert(err), vim.log.levels.WARN)
            -- The project client resolves nothing here, but beats no LSP at all.
            return on_dir()
        end
        start(bufnr, python)
    end)
end

vim.api.nvim_create_autocmd("LspDetach", {
    group = vim.api.nvim_create_augroup("my.uv_script_env", { clear = true }),
    callback = function(args)
        local client = vim.lsp.get_client_by_id(args.data.client_id)
        -- Neovim stops no client of its own accord, and the project's must live
        -- on: stopping it would force a cold re-index of the whole workspace.
        if not client or not pin(client) then return end
        -- The event comes before the buffer leaves `attached_buffers`, so the
        -- count only reaches zero with it excluded. Checking on a `vim.schedule`
        -- instead would race with a re-attach.
        for bufnr in pairs(client.attached_buffers) do
            if bufnr ~= args.buf then return end
        end
        client:stop()
    end,
})

return M
