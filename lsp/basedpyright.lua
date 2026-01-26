-- Alts are pyright, pylsp, jedi_language_server, etc.
-- https://old.reddit.com/r/neovim/comments/1bh0kba/psa_new_python_lsp_that_supports_inlay_hints_and/
-- https://github.com/neovim/nvim-lspconfig/blob/master/doc/configs.md#basedpyright

-- :Conda command to switch conda environment for Python LSP
local function get_conda_envs()
    local handle = io.popen("conda env list --json 2>/dev/null")
    if not handle then return {} end
    local output = handle:read("*a")
    handle:close()
    local ok, data = pcall(vim.json.decode, output)
    if not ok or not data.envs then return {} end
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
    local handle = io.popen("conda info --json 2>/dev/null")
    if not handle then
        vim.notify("Failed to run conda", vim.log.levels.ERROR)
        return
    end
    local output = handle:read("*a")
    handle:close()
    local ok, info = pcall(vim.json.decode, output)
    if not ok then
        vim.notify("Failed to parse conda info", vim.log.levels.ERROR)
        return
    end
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
    vim.cmd("PyrightSetPythonPath " .. python_path)
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
-- There is also :PyrightSetPythonPath to set it manually.
if vim.env.CONDA_PREFIX ~= nil then
    pythonPath = vim.env.CONDA_PREFIX .. '/bin/python'
end

return {
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
                stubPath = "~/.config/nvim/lsp/python_stubs/",
            }
        }
    }
}
