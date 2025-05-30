#!/usr/bin/env lua

-- replace the default lsp diagnostic letters with prettier symbols
vim.fn.sign_define("LspDiagnosticsSignError", { text = "", numhl = "LspDiagnosticsDefaultError" })
vim.fn.sign_define("LspDiagnosticsSignWarning", { text = "", numhl = "LspDiagnosticsDefaultWarning" })
vim.fn.sign_define("LspDiagnosticsSignInformation", { text = "", numhl = "LspDiagnosticsDefaultInformation" })
vim.fn.sign_define("LspDiagnosticsSignHint", { text = "", numhl = "LspDiagnosticsDefaultHint" })

local function nmap(lhs, rhs, desc)
    vim.keymap.set('n', lhs, rhs, { desc = desc })
end

-- LSP and completion status, overall conf etc.
nmap("<leader>li", "<Cmd>LspInfo<CR>", "Info")
-- can't use backspace since it is hardcoded by mini.clue for up one level
nmap("<leader>l0", "<Cmd>LspStop<CR>", "Stop")
nmap("<leader>l1", "<Cmd>LspStart<CR>", "Start")
nmap("<leader>l!", "<Cmd>LspRestart<CR>", "Restart")
nmap("<leader>lL", "<Cmd>LspLog<CR>", "Log")
-- match other completion related entries under x
nmap("<leader>xS", "<Cmd>CmpStatus<CR>", "Cmp status")

-- See `:help vim.diagnostic.*` for documentation on any of the below functions
nmap('<leader>dd', vim.diagnostic.open_float, "Line diagnostic")
-- can't use backspace since it is hardcoded by mini.clue for up one level
nmap('<leader>d0', function() vim.diagnostic.enable(false) end, "Disable diagnostics")
nmap('<leader>d1', vim.diagnostic.enable, "Enable diagnostics")
nmap('[d', function() vim.diagnostic.jump { count = -1, float = true } end, "Diagnostic")
nmap(']d', function() vim.diagnostic.jump { count = 1, float = true } end, "Diagnostic")
nmap('<leader>dl', vim.diagnostic.setloclist, "Loclist diagnostics")

local grp = vim.api.nvim_create_augroup("my.lsp", {})

vim.api.nvim_create_autocmd('LspAttach', {
    group = grp,
    callback = function(args)
        local function map(desc, keys, func, mode)
            mode = mode or 'n'
            vim.keymap.set(mode, keys, func, { buffer = args.buf, desc = desc })
        end
        local function mapfzf(desc, keys, funcname, mode)
            map(desc, keys, function()
                require "fzf-lua"["lsp_" .. funcname]()
            end, mode)
        end

        local client = assert(vim.lsp.get_client_by_id(args.data.client_id))

        -- See `:help vim.lsp.*` for documentation on any of the below functions
        -- TODO: have treesitter fallback for things like go to references for when LSP is not attached.

        if client:supports_method('textDocument/references') then
            -- `vim.lsp.buf.references`
            mapfzf("Goto references", "gr", "references")
        end

        if client:supports_method('textDocument/definition') then
            -- `vim.lsp.buf.definition`
            mapfzf("Goto definition", "gd", "definitions")
        end

        if client:supports_method('textDocument/declaration') then
            -- `vim.lsp.buf.declaration`
            mapfzf("Goto declaration", "gD", "declaration")
        end

        -- `vim.lsp.buf.type_definition`
        mapfzf("Goto type definition", "g<C-d>", "typedefs")

        -- gi is for goto last insert and switch to insert mode, and gI is to insert at first column <count> times.
        -- We rarely use go to implementation though
        mapfzf("Goto implementation", "g<C-i>", "implementations")
        mapfzf("Goto defs+refs+impl+...", "g<C-S-d>", "finder")

        -- `vim.lsp.buf.code_action`
        mapfzf("Code action", '<leader>la', "code_actions")

        mapfzf("Symbols", '<leader>ls', "document_symbols")
        mapfzf("Symbols", '<leader>ws', "workspace_symbols")
        -- a for all diagnostics, no filtering
        mapfzf("Diagnostics", '<leader>da', "document_diagnostics")
        mapfzf("Diagnostics", '<leader>wd', "workspace_diagnostics")

        map("Hover", 'K', vim.lsp.buf.hover)
        map("Signature", 'gh', vim.lsp.buf.signature_help)
        map("Rename", '<leader>rn', vim.lsp.buf.rename)
        -- currently autoformatting on save (see below). Maybe TODO have a way to toggle this.
        vim.keymap.set({ 'n', 'v' }, '<leader>lf', vim.lsp.buf.format, { buffer = args.buf, desc = "Format" })

        -- Builtin auto-completion disabled for now since blink.cmp claims to be faster.
        -- Enable auto-completion. Note: Use CTRL-Y to select an item. |complete_CTRL-Y|
        -- if client:supports_method('textDocument/completion') then
        --     -- Optional: trigger autocompletion on EVERY keypress. May be slow!
        --     -- local chars = {}; for i = 32, 126 do table.insert(chars, string.char(i)) end
        --     -- client.server_capabilities.completionProvider.triggerCharacters = chars
        --     vim.lsp.completion.enable(true, client.id, args.buf, {autotrigger = true})
        -- end

        -- Auto-format ("lint") on save.
        -- Usually not needed if server supports "textDocument/willSaveWaitUntil".
        if not client:supports_method('textDocument/willSaveWaitUntil')
            and client:supports_method('textDocument/formatting') then
            vim.api.nvim_create_autocmd('BufWritePre', {
                group = grp,
                buffer = args.buf,
                callback = function()
                    vim.lsp.buf.format({ bufnr = args.buf, id = client.id, timeout_ms = 1000 })
                end,
            })
        end
    end,
})
