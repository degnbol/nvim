-- Alts are pyright, pylsp, jedi_language_server, etc.
-- https://old.reddit.com/r/neovim/comments/1bh0kba/psa_new_python_lsp_that_supports_inlay_hints_and/
-- https://github.com/neovim/nvim-lspconfig/blob/master/doc/configs.md#basedpyright

local pythonPath
-- Use conda python if available.
if vim.env.CONDA_PREFIX ~= nil then
    pythonPath = vim.env.CONDA_PREFIX .. '/bin/python'
end
-- There is also :PyrightSetPythonPath to set it on the fly

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
