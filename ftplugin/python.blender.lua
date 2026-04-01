-- Add blender stubs (fake-bpy-module) to basedpyright extraPaths.
-- Triggered by compound filetype "python.blender" (set in ftdetect/blender.lua).
-- Reinstall: uv pip install --target lsp_ext/blender-stubs fake-bpy-module-latest
local stubs = vim.fn.stdpath("config") .. "/lsp_ext/blender-stubs" -- no API equivalent

vim.api.nvim_create_autocmd("LspAttach", {
    buffer = 0,
    once = true,
    callback = function(args)
        local client = vim.lsp.get_client_by_id(args.data.client_id)
        if not client or client.name ~= "basedpyright" then return end
        local settings = client.settings
        local paths = settings.basedpyright.analysis.extraPaths or {}
        if not vim.tbl_contains(paths, stubs) then
            table.insert(paths, stubs)
            settings.basedpyright.analysis.extraPaths = paths
            client:notify("workspace/didChangeConfiguration", { settings = settings })
        end
    end,
})
