local ensure_installed
if vim.fn.executable("npm") == 1 and vim.fn.executable("cargo") == 1 then
    ensure_installed = {
        "bashls",
        -- https://github.com/neovim/nvim-lspconfig/blob/master/doc/configs.md#omnisharp
        -- "csharp_ls", -- instead see install/csharp.sh. There's also omnisharp on other LSPs on Mason.
        -- python:
        -- "jedi_language_server",
        -- "pyright",
        -- "pylsp",
        "basedpyright",
        -- "julials", -- Manual
        -- "ltex", -- grammar check for latex, markdown, etc
        "texlab",
        -- "r_language_server",
        "lua_ls",
        -- "vimls",
        "rust_analyzer",
        -- "ts_ls", -- javascript. MasonInstall typescript-language-server
        "sqlls",
        -- "latexindent",
        -- "matlab_ls",
        -- "kotlin-language-server",
        "yamlls",
        -- "awk_ls",
        -- "cypher_ls", -- neo4j
        -- spell check/autocorrectors:
        -- "misspell",
        -- "typos",
        -- "typos_lsp",
    }
end

return {
    -- scala (not yet in Mason)
    {
        "scalameta/nvim-metals",
        ft = "scala",
        dependencies = { "nvim-lua/plenary.nvim" }
    },
    {
        "neovim/nvim-lspconfig",
        config = function()
            -- TODO: understand thefoldingRange capabilities. What are they and where should they be added?
            -- Look at https://github.com/kevinhwang91/nvim-ufo and
            -- https://cmp.saghen.dev/installation.html#merging-lsp-capabilities
            local capabilities
            local using_blink, blink = pcall(require, "blink.cmp")
            if using_blink then
                capabilities = blink.get_lsp_capabilities({
                    textDocument = { completion = { completionItem = { snippetSupport = true } } },
                })
            else
                capabilities = require('cmp_nvim_lsp').default_capabilities(vim.lsp.protocol.make_client_capabilities())
            end

            -- added for https://github.com/kevinhwang91/nvim-ufo see ufo-conf etc
            capabilities.textDocument.foldingRange = {
                dynamicRegistration = false,
                lineFoldingOnly = true
            }

            -- now for adding the language servers
            local lsp = require "lspconfig"

            -- add to lsp default config
            lsp.util.default_config = vim.tbl_deep_extend('force', lsp.util.default_config, {
                capabilities = capabilities,
                on_attach = on_attach
            })

            -- naming: https://github.com/williamboman/mason-lspconfig.nvim/blob/main/doc/server-mapping.md
            -- config help: https://github.com/neovim/nvim-lspconfig/blob/master/doc/configs.md

            -- Doesn't work as well as tinymist, e.g. for multiple files and goto def.
            -- lsp.typst_lsp.setup {
            --     settings = {
            --         exportPdf = "never", -- Choose onType, onSave or never. Toggle in ftplugin/typst.lua
            --         -- serverPath = "", -- Normally, there is no need to uncomment it.
            --     },
            --     -- LSP can only understand multiple files if main.typ is opened first and then other files after.
            --     -- Tried pinning the main:
            --     -- https://github.com/nvarner/typst-lsp/issues/366
            -- }

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

            -- grammar. Not yet supporting latex but supports typst.
            -- lsp.harper_ls.setup {
            --     settings = {
            --         ["harper-ls"] = {
            --             userDictPath = vim.opt.runtimepath:get()[1] .. "/spell/custom.utf8.add",
            --         }
            --     },
            --     filetypes = {
            --         "tex", "typst"
            --     }
            -- }
        end
    },
    -- add :LspInstall <language> and :Mason for conveniently installing LSP language specific servers
    {
        "williamboman/mason.nvim",
        build = ":MasonUpdate",
        lazy = true, -- load as mason-lspconfig dep
        config = true,
    },
    {
        "williamboman/mason-lspconfig.nvim",
        -- NOTE: can't be lazy for some lsp servers to work, e.g. lua.
        -- cmd = {"Mason", "MasonUpdate", "MasonInstall", "MasonLog", "MasonUninstall", "MasonUninstallAll"},
        dependencies = { "neovim/nvim-lspconfig", "williamboman/mason.nvim" },
        -- naming: https://github.com/williamboman/mason-lspconfig.nvim/blob/main/doc/server-mapping.md
        -- config help: https://github.com/neovim/nvim-lspconfig/blob/master/doc/configs.md
        -- also see lsp.lua
        opts = {
            ensure_installed = ensure_installed
        }
    },
    -- flutter tools contains LSP for dart.
    -- The FlutterRun doesn't seem to work, so currently only using the flutter packaged dart lsp,
    -- and running flutter run --flavor dev --debug from the terminal, while
    -- android studio has opened an emulator.
    {
        "akinsho/flutter-tools.nvim",
        ft = "dart",
        dependencies = {
            "neovim/nvim-lspconfig",
            "nvim-lua/plenary.nvim",
            'stevearc/dressing.nvim', -- optional for vim.ui.select
        },
        config = function()
            -- https://github.com/akinsho/flutter-tools.nvim#full-configuration
            require "flutter-tools".setup {
                -- there are other fields under lsp than capabilities and
                -- on_attach, but this is a fine way to write it as long as we
                -- don't change the other settings.
                lsp = require "lspconfig".util.default_config
            }
        end
    },

    {
        'saecki/crates.nvim',
        event = { "BufRead Cargo.toml" },
        config = function()
            local crates = require "crates"
            crates.setup {
                completion = {
                    cmp = { enabled = false },
                    crates = {
                        enabled = true,
                        max_results = 8, -- The maximum number of search results to display
                        min_chars = 3,   -- The minimum number of charaters to type before completions begin appearing
                    },
                },
                lsp = {
                    enabled = true,
                    on_attach = function(client, bufnr)
                        -- the same on_attach function as for your other lsp's
                    end,
                    actions = true,
                    completion = true,
                    hover = true,
                },
            }

            local function localmap(key, func, desc)
                vim.keymap.set("n", "<LocalLeader>" .. key, crates[func], { desc = desc })
            end

            localmap("t", "toggle", "Toggle crates")
            localmap("r", "reload", "Reload crates")

            localmap("v", "show_versions_popup", "Show versions")
            localmap("f", "show_features_popup", "Show features")
            localmap("d", "show_dependencies_popup", "Show dependencies")

            localmap("u", "update_crate", "Update")
            localmap("u", "update_crates", "Update")
            localmap("a", "update_all_crates", "Update all")
            localmap("U", "upgrade_crate", "Upgrade")
            localmap("U", "upgrade_crates", "Upgrade")
            localmap("A", "upgrade_all_crates", "Upgrade all")

            localmap("x", "expand_plain_crate_to_inline_table", "Plain crate -> inline table")
            localmap("X", "extract_crate_into_table", "Crate -> table")

            localmap("H", "open_homepage", "Homepage")
            localmap("R", "open_repository", "Repo")
            localmap("D", "open_documentation", "Documentation")
            localmap("C", "open_crates_io", "crates.io")
            localmap("L", "open_lib_rs", "lib.rs")
        end,
    },
}
