local map = require "utils/keymap"

-- TODO: we have fzf and qf versions at gr* and <leader>l.
-- Should we default to qf and have a quick qf specific keymap that moves to fzf?
-- fzf is for quickly finding needle in haystack, and qf is to navigate between things.
-- qf is better for references.
-- fzf is better for finding word in repo, or finding file with certain name.
-- Maybe just use most appropriate tool for the case then.
-- And then have convenient keymap to switch between them.
-- Check out fzf-lua alts.

local grp = vim.api.nvim_create_augroup("my.lsp", { clear = true })

vim.api.nvim_create_autocmd('LspAttach', {
    group = grp,
    callback = function(args)
        local opts = { buffer = args.buf }

        local function map_fzf(lhs, funcname, desc)
            map.n(lhs, function()
                require "fzf-lua"["lsp_" .. funcname]()
            end, desc, opts)
        end

        local client = assert(vim.lsp.get_client_by_id(args.data.client_id))

        -- See `:help vim.lsp.*` for documentation on any of the below functions
        -- TODO: have treesitter fallback for things like go to references for when LSP is not attached.

        if client:supports_method('textDocument/references') then
            -- `vim.lsp.buf.references`
            map_fzf("<leader>lr", "references", "References")
        end

        if client:supports_method('textDocument/definition') then
            -- `vim.lsp.buf.definition`
            map_fzf("gd", "definitions", "Definition")
        end

        if client:supports_method('textDocument/declaration') then
            -- `vim.lsp.buf.declaration`
            map_fzf("gD", "declaration", "Declaration")
        end

        -- `vim.lsp.buf.type_definition`
        map_fzf("g<C-d>", "typedefs", "Type definition")

        -- gi is for goto last insert and switch to insert mode, and gI is to insert at first column <count> times.
        -- We rarely use go to implementation though
        map_fzf("g<C-i>", "implementations", "Implementation")
        map_fzf("g<C-S-d>", "finder", "Defs+refs+impl+...")

        -- `vim.lsp.buf.code_action`
        map_fzf('<leader>la', "code_actions", "Code action")

        map_fzf('<leader>ls', "document_symbols", "Symbols")
        map_fzf('<leader>ws', "workspace_symbols", "Symbols")
        -- a for all diagnostics, no filtering
        map_fzf('<leader>da', "document_diagnostics", "Diagnostics")
        map_fzf('<leader>wd', "workspace_diagnostics", "Diagnostics")

        map.n('gh', vim.lsp.buf.signature_help, "Signature", opts)
        -- TODO: maybe have a way to enable auto format on write or make an explicit keybind that does format and write both.
        map.nx('grf', vim.lsp.buf.format, "Format", opts)

        -- Builtin auto-completion disabled for now since blink.cmp claims to be faster.
        -- Enable auto-completion. Note: Use CTRL-Y to select an item. |complete_CTRL-Y|
        -- if client:supports_method('textDocument/completion') then
        --     -- Optional: trigger autocompletion on EVERY keypress. May be slow!
        --     -- local chars = {}; for i = 32, 126 do table.insert(chars, string.char(i)) end
        --     -- client.server_capabilities.completionProvider.triggerCharacters = chars
        --     vim.lsp.completion.enable(true, client.id, args.buf, {autotrigger = true})
        -- end

        -- Auto-format on save.
        -- Usually not needed if server supports "textDocument/willSaveWaitUntil".
        -- if not client:supports_method('textDocument/willSaveWaitUntil')
        --     and client:supports_method('textDocument/formatting') then
        --     vim.api.nvim_create_autocmd('BufWritePre', {
        --         group = grp,
        --         buffer = args.buf,
        --         callback = function()
        --             -- Short timeout so it doesn't hang.
        --             -- If it needs more time then invoke formatting manually.
        --             -- pcall so it doesn't complain if it didn't finish in time.
        --             pcall(vim.lsp.buf.format, { bufnr = args.buf, id = client.id, timeout_ms = 400 })
        --         end,
        --     })
        -- end
    end,
})
