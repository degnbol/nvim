local map = require "utils/keymap"

-- Ours replaces nvim-lspconfig's on_attach (its `Lsp*` export/pin commands) in
-- the lsp/ merge, so chain it.
local this_file = vim.uv.fs_realpath(debug.getinfo(1, "S").source:sub(2))
local base_on_attach
for _, f in ipairs(vim.api.nvim_get_runtime_file("lsp/tinymist.lua", true)) do
    if vim.uv.fs_realpath(f) ~= this_file then
        base_on_attach = dofile(f).on_attach or base_on_attach
    end
end

-- https://github.com/Myriad-Dreamin/tinymist/blob/main/editors/neovim/Configuration.md
return {
    on_attach = function (client, bufnr)
        if base_on_attach then base_on_attach(client, bufnr) end
        -- Pinning a main file is load-bearing, not a nicety: tinymist's ref/label
        -- features (hover on @cite/@fig/@tbl/@heading, goto-def on an uncompiled
        -- doc) read `ctx.success_doc()`, which is only populated for a pinned
        -- main -- with no pin they silently return nil.
        local function pinMain(fname)
            client:exec_cmd { title = "Pin main", command = "tinymist.pinMain", arguments = { fname } }
        end
        map.n('<LocalLeader>p', function ()
            return pinMain(vim.api.nvim_buf_get_name(0))
        end, "Pin buffer as main", { buffer=true})
        -- main.typ upward for a multi-file project, else the buffer itself so a
        -- standalone doc's refs still resolve.
        local mainfile = vim.fs.find("main.typ", {type="file", upward=true})[1]
        pinMain(mainfile or vim.api.nvim_buf_get_name(bufnr))
    end,
}
