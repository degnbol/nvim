#!/usr/bin/env lua
return {
    -- scala (not yet in Mason)
    {
        "scalameta/nvim-metals",
        ft = "scala",
        dependencies = { "nvim-lua/plenary.nvim" }
    },
    {
        "neovim/nvim-lspconfig",
        dependencies = "folke/neodev.nvim", -- signature on nvim lua calls which helps in messing with plugins etc
        config = function()
            -- default config copied from https://github.com/neovim/nvim-lspconfig
            -- inspiration from https://vonheikemen.github.io/devlog/tools/setup-nvim-lspconfig-plus-nvim-cmp/

            -- call before rest of lspconfig
            require "neodev".setup {}

            -- hide diagnostics for hints and information.
            vim.lsp.handlers["textDocument/publishDiagnostics"] = vim.lsp.with(
                vim.lsp.diagnostic.on_publish_diagnostics, {
                    signs = {
                        severity_limit = 'Warning',
                    },
                    underline = false,
                    update_in_insert = false,
                    virtual_text = {
                        spacing = 4,
                        severity_limit = 'Error',
                    },
                }
            )

            -- replace the default lsp diagnostic letters with prettier symbols
            vim.fn.sign_define("LspDiagnosticsSignError", { text = "", numhl = "LspDiagnosticsDefaultError" })
            vim.fn.sign_define("LspDiagnosticsSignWarning", { text = "", numhl = "LspDiagnosticsDefaultWarning" })
            vim.fn.sign_define("LspDiagnosticsSignInformation", { text = "", numhl = "LspDiagnosticsDefaultInformation" })
            vim.fn.sign_define("LspDiagnosticsSignHint", { text = "", numhl = "LspDiagnosticsDefaultHint" })

            -- now for adding the language servers
            local lsp = require "lspconfig"


            -- See `:help vim.diagnostic.*` for documentation on any of the below functions
            vim.keymap.set('n', '<space>dd', vim.diagnostic.open_float, {desc = "Line diagnostic"})
            vim.keymap.set('n', '<space>D', "<cmd>Telescope diagnostics<CR>", { desc = "Diagnostics" })
            vim.keymap.set('n', '[d', vim.diagnostic.goto_prev, { desc = "Diagnostic" })
            vim.keymap.set('n', ']d', vim.diagnostic.goto_next, { desc = "Diagnostic" })
            vim.keymap.set('n', '<space>dl', vim.diagnostic.setloclist, { desc = "Loclist diagnostics" })

            -- Use an on_attach function to only map the following keys
            -- after the language server attaches to the current buffer
            local on_attach = function(client, bufnr)
                -- Enable completion triggered by <c-x><c-o>
                vim.api.nvim_buf_set_option(bufnr, 'omnifunc', 'v:lua.vim.lsp.omnifunc')
                -- See `:help vim.lsp.*` for documentation on any of the below functions
                vim.keymap.set('n', 'gD', vim.lsp.buf.declaration, { buffer = bufnr, desc = "Goto declaration" })
                vim.keymap.set('n', 'gd', vim.lsp.buf.definition, { buffer = bufnr, desc = "Goto definition" })
                vim.keymap.set('n', 'K', vim.lsp.buf.hover, { buffer = bufnr, desc = "Hover" })
                -- gi is for goto last insert and switch to insert mode so we use gI
                vim.keymap.set('n', 'gI', vim.lsp.buf.implementation, { buffer = bufnr, desc = "Goto implementation" })
                vim.keymap.set({'i', 'n'}, '<C-S-k>', vim.lsp.buf.signature_help, { buffer = bufnr, desc = "Signature" })
                vim.keymap.set('n', '<space>wa', vim.lsp.buf.add_workspace_folder, { buffer = bufnr, desc = "Add" })
                vim.keymap.set('n', '<space>wr', vim.lsp.buf.remove_workspace_folder, { buffer = bufnr, desc = "Remove" })
                vim.keymap.set('n', '<space>wl', function() print(vim.inspect(vim.lsp.buf.list_workspace_folders())) end, { buffer = bufnr, desc = "List" })
                vim.keymap.set('n', '<space>dt', vim.lsp.buf.type_definition, { buffer = bufnr, desc = "Type def" })
                vim.keymap.set('n', '<space>rn', vim.lsp.buf.rename, { buffer = bufnr, desc = "Rename" })
                vim.keymap.set('n', '<space>ca', vim.lsp.buf.code_action, { buffer = bufnr, desc = "Code action" })
                -- vim.keymap.set('n', 'gr', vim.lsp.buf.references, { buffer=bufnr, desc="Goto references" })
                vim.keymap.set('n', 'gr', "<cmd>Telescope lsp_references<CR>", { buffer = bufnr, desc = "Goto references" })
                vim.keymap.set('n', '<space>F', vim.lsp.buf.format, { buffer = bufnr, desc = "Format" })
            end

            local capabilities = require('cmp_nvim_lsp').default_capabilities(vim.lsp.protocol.make_client_capabilities())

            -- added for https://github.com/kevinhwang91/nvim-ufo see ufo-conf etc
            capabilities.textDocument.foldingRange = {
                dynamicRegistration = false,
                lineFoldingOnly = true
            }

            -- add to lsp default config
            lsp.util.default_config = vim.tbl_deep_extend('force', lsp.util.default_config, {
                capabilities = capabilities,
                on_attach = on_attach
            })

            -- naming: https://github.com/williamboman/mason-lspconfig.nvim/blob/main/doc/server-mapping.md
            -- config help: https://github.com/neovim/nvim-lspconfig/blob/master/doc/server_configurations.md

            lsp.bashls.setup { filetypes = { "sh", "bash", "zsh" } }
            -- lsp.pyright.setup { }
            -- lsp.pylsp.setup { }
            lsp.jedi_language_server.setup {}
            lsp.julials.setup {
                on_new_config = function(new_config, _)
                    local julia = vim.fn.expand("~/.julia/environments/nvim-lspconfig/bin/julia")
                    -- check if we have made the dedicated julia env which should have a custom system image
                    if require 'lspconfig'.util.path.is_file(julia) then
                        -- regular julia seem to work
                        -- new_config.cmd[1] = julia
                    end
                end
            }
            lsp.r_language_server.setup {}
            lsp.vimls.setup {}
            lsp.lua_ls.setup {
                settings = {
                    Lua = {
                        diagnostics = {
                            -- https://www.reddit.com/r/neovim/comments/11mdxex/lspconfig_cannot_access_configuration_for_lua_ls/
                            -- Get the language server to recognize the `vim` global
                            globals = { 'vim' },
                        },
                        workspace = {
                            -- https://www.reddit.com/r/neovim/comments/11mdxex/lspconfig_cannot_access_configuration_for_lua_ls/
                            -- Make the server aware of Neovim runtime files
                            library = vim.api.nvim_get_runtime_file("", true),
                            -- stopped a prompt at startup asking to configure project as luv.
                            checkThirdParty = false
                        },
                        completion = { callSnippet = "Replace" },
                    }
                }
            }
            lsp.marksman.setup { filetypes = { "markdown" } }
            -- seems ltex has more text description for functions but texlab has more functions so I use both in combination
            lsp.ltex.setup {
                filetypes = { "tex" }, -- active for markdown as well by default which crashes
                on_attach = function(client, bufnr)
                    on_attach(client, bufnr)
                    -- hacky. VimtexErrors puts errors found by Vimtex in quickfix (should be
                    -- running, use <leader>Lb) then cclose closes quickfix, and then Telescope
                    -- opens the quickfix in a nicer view.
                    vim.keymap.set('n', '<space>E', "<Cmd>VimtexErrors<CR>|:cclose|<Cmd>Telescope quickfix<CR>", opts)
                end }
            lsp.texlab.setup {}

            lsp.csharp_ls.setup {
                -- AutomaticWorkspaceInit = true,
            }

            lsp.rust_analyzer.setup {}

            -- MasonInstall typescript-language-server
            -- Used on javascript as well.
            lsp.tsserver.setup {}

            -- scala. Linked as minimal setup from on nvim-metals git:
            -- https://github.com/scalameta/nvim-metals/discussions/39
            local metals_config = require "metals".bare_config()
            metals_config.capabilities = capabilities
            -- Autocmd that will actually be in charging of starting the whole thing.
            -- Create a .sc file, then reopen it. There should be a MetalsInstall warning.
            -- Run MetalsInstall, wait patiently until it gives a new message and restart.
            local nvim_metals_group = vim.api.nvim_create_augroup("nvim-metals", { clear = true })
            vim.api.nvim_create_autocmd("FileType", {
                pattern = { "scala", "sbt" },
                callback = function()
                    require("metals").initialize_or_attach(metals_config)
                end,
                group = nvim_metals_group,
            })
        end
    },
    -- add :LspInstall <language> and :Mason for conveniently installing LSP language specific servers
    {
        "williamboman/mason.nvim",
        build = ":MasonUpdate",
        lazy = true,                                                   -- load as mason-lspconfig dep
        config = true,
    },
    {
        "williamboman/mason-lspconfig.nvim",
        -- NOTE: can't be lazy for some lsp servers to work, e.g. lua
        -- cmd = {"Mason", "MasonUpdate", "MasonInstall", "MasonLog", "MasonUninstall", "MasonUninstallAll"},
        dependencies = { "neovim/nvim-lspconfig", "williamboman/mason.nvim" },
        -- naming: https://github.com/williamboman/mason-lspconfig.nvim/blob/main/doc/server-mapping.md
        -- config help: https://github.com/neovim/nvim-lspconfig/blob/master/doc/server_configurations.md
        -- also see lsp.lua
        opts = {
            ensure_installed = {
                "bashls",
                -- https://github.com/neovim/nvim-lspconfig/blob/master/doc/server_configurations.md#omnisharp
                "csharp_ls",
                "jedi_language_server",
                "pyright",
                "pylsp",
                "julials",
                "ltex",
                "texlab",
                "r_language_server",
                "lua_ls",
                "vimls",
                "rust_analyzer",
                "tsserver",     -- javascript. MasonInstall typescript-language-server
                "sqlls",
                -- "latexindent",
                -- "matlab-language-server",
                -- "kotlin-language-server",
            }
        }
    },
}
