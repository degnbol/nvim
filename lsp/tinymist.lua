local map = require "utils/keymap"

-- https://github.com/Myriad-Dreamin/tinymist/blob/main/editors/neovim/Configuration.md
-- for typst
return {
    on_attach = function (client, bufnr)
        -- Don't override default on_attach
        require"lspconfig".util.default_config.on_attach(client, bufnr)

        local function pinMain(fname)
            for _, client in ipairs(vim.lsp.get_clients{name="tinymist"}) do
                client:exec_cmd { command = 'tinymist.pinMain', arguments = { fname } }
            end
        end
        map.n('<LocalLeader>p', function ()
            return pinMain(vim.api.nvim_buf_get_name(0))
        end, "Pin buffer as main", { buffer=true})
        -- search upwards for a main.typ
        local mainfile = vim.fs.find("main.typ", {type="file", upward=true})[1]
        if mainfile ~= nil then pinMain(mainfile) end
    end,
}
