-- Alts are pyright, pylsp, jedi_language_server, etc.
-- https://old.reddit.com/r/neovim/comments/1bh0kba/psa_new_python_lsp_that_supports_inlay_hints_and/
-- https://github.com/neovim/nvim-lspconfig/blob/master/doc/configs.md#basedpyright
local uv_script = require "autocmds/uv_script_env"
local util = require "utils/init"

---Run `conda <args> --json` and decode its output, notifying on failure.
---@param args string[] conda subcommand and its arguments, without `--json`
---@return table|nil output decoded JSON, nil on failure
local function conda_json(args)
    local cmd = { "conda", unpack(args) }
    table.insert(cmd, "--json")
    local obj = vim.system(cmd, { text = true }):wait()
    if obj.code ~= 0 then
        util.notify_failure(cmd, obj)
        return nil
    end
    local ok, data = pcall(vim.json.decode, obj.stdout)
    if not ok then
        vim.notify("Failed to parse conda output: " .. data, vim.log.levels.ERROR)
        return nil
    end
    return data
end

-- :Conda command to switch conda environment for Python LSP
local function get_conda_envs()
    local data = conda_json({ "env", "list" })
    if not data or not data.envs then return {} end
    local envs = {}
    for _, path in ipairs(data.envs) do
        table.insert(envs, vim.fn.fnamemodify(path, ":t"))
    end
    return envs
end

vim.api.nvim_create_user_command("Conda", function(opts)
    local env = opts.args
    if env == "" then
        if vim.env.CONDA_PREFIX then
            vim.notify("Current: " .. vim.fn.fnamemodify(vim.env.CONDA_PREFIX, ":t"), vim.log.levels.INFO)
        else
            vim.notify("No conda env active", vim.log.levels.INFO)
        end
        return
    end
    -- Get conda info to find envs directory
    local info = conda_json({ "info" })
    if not info then return end
    -- Find the env path
    local env_path
    for _, path in ipairs(info.envs or {}) do
        if vim.fn.fnamemodify(path, ":t") == env then
            env_path = path
            break
        end
    end
    if not env_path then
        vim.notify("Conda env not found: " .. env, vim.log.levels.ERROR)
        return
    end
    local python_path = env_path .. "/bin/python"
    -- nvim-lspconfig creates this buffer-locally on attach, so run it from a
    -- python buffer.
    vim.cmd.LspPyrightSetPythonPath(python_path)
    vim.notify("Set Python: " .. python_path, vim.log.levels.INFO)
end, {
    nargs = "?",
    complete = function() return get_conda_envs() end,
    desc = "Set conda environment for Python LSP",
})

local pythonPath
-- Use conda python if available.
-- Actually no, this will switch to conda when uv would be found nicely by default.
-- TODO: maybe have a way of auto setting it to conda only if uv is not active?
-- There is also :LspPyrightSetPythonPath to set it manually.
if vim.env.CONDA_PREFIX ~= nil then
    pythonPath = vim.env.CONDA_PREFIX .. '/bin/python'
end

return {
    filetypes = { "python", "python.blender" },
    -- A PEP 723 script gets a client of its own, pointed at its uv environment.
    root_dir = uv_script.root_dir,
    -- A field rather than only an argument at our own `vim.lsp.start`:
    -- `start_config` forwards `config.reuse_client`, and without it the *project*
    -- client falls back to neovim's default, which never compares interpreters.
    reuse_client = uv_script.same_env,
    settings = {
        python = {
            pythonPath = pythonPath
        },
        basedpyright = {
            analysis = {
                -- defaults to complaining about unknown types, and we don't want to be reminded to specify types.
                -- Plus when using other's code that we can't change there will also be warnings about their lack of type declaration.
                -- https://detachhead.github.io/basedpyright/#/configuration
                typeCheckingMode = "standard",
                stubPath = vim.fn.stdpath("config") .. "/lsp_ext/python_stubs/",
                extraPaths = (function()
                    local paths = { "src" }
                    vim.list_extend(paths, vim.fn.glob(
                        vim.fn.stdpath("config") .. "/lsp_ext/extraPaths/*/", false, true))
                    -- PEP 561 stub packages (scipy-stubs, pandas-stubs, ...)
                    -- installed via lsp_ext/python_stubs_pypi/RUNME.sh.
                    local pypi_lib = vim.fn.expand("~/.local/share/python-stubs/lib")
                    for _, py in ipairs(vim.fs.find(
                        function(name) return name:match("^python%d+%.%d+$") end,
                        { path = pypi_lib, type = "directory", limit = math.huge })
                    ) do
                        table.insert(paths, py .. "/site-packages")
                    end
                    return paths
                end)(),
                -- Prevent goto-definition from landing in build/ directories.
                -- https://docs.basedpyright.com/v1.20.0/configuration/language-server-settings/
                exclude = { "**/build" },
                diagnosticSeverityOverrides = {
                    reportUnusedCallResult = "none",
                    reportMissingModuleSource = "none",
                },
            }
        }
    }
}
