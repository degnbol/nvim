local map = require "utils/keymap"

-- https://github.com/Myriad-Dreamin/tinymist/blob/main/editors/neovim/Configuration.md
-- NOTE: this on_attach is dead on its own -- nvim-lspconfig ships its own
-- lsp/tinymist.lua whose on_attach wins the runtimepath merge over ours (see the
-- `tbl_deep_extend` note in the neovim skill). It is re-applied, chained after
-- nvim-lspconfig's, in lua/plugins/lsp.lua's mason-lspconfig `after` block.
return {
    on_attach = function (client, bufnr)
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
