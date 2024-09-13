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

            -- replace the default lsp diagnostic letters with prettier symbols
            vim.fn.sign_define("LspDiagnosticsSignError", { text = "", numhl = "LspDiagnosticsDefaultError" })
            vim.fn.sign_define("LspDiagnosticsSignWarning", { text = "", numhl = "LspDiagnosticsDefaultWarning" })
            vim.fn.sign_define("LspDiagnosticsSignInformation", { text = "", numhl = "LspDiagnosticsDefaultInformation" })
            vim.fn.sign_define("LspDiagnosticsSignHint", { text = "", numhl = "LspDiagnosticsDefaultHint" })

            -- now for adding the language servers
            local lsp = require "lspconfig"

            -- See `:help vim.diagnostic.*` for documentation on any of the below functions
            vim.keymap.set('n', '<space>dd', vim.diagnostic.open_float, { desc = "Line diagnostic" })
            vim.keymap.set('n', '<space>dk', vim.diagnostic.disable, { desc = "Disable diagnostics" })
            vim.keymap.set('n', '<space>ds', vim.diagnostic.enable, { desc = "Enable diagnostics" })
            vim.keymap.set('n', '<space>fd', "<cmd>Telescope diagnostics<CR>", { desc = "Diagnostics" })
            vim.keymap.set('n', '[d', vim.diagnostic.goto_prev, { desc = "Diagnostic" })
            vim.keymap.set('n', ']d', vim.diagnostic.goto_next, { desc = "Diagnostic" })
            vim.keymap.set('n', '<space>dl', vim.diagnostic.setloclist, { desc = "Loclist diagnostics" })
            -- LSP shorthand keymaps
            vim.keymap.set('n', '<leader>li', "<Cmd>LspInfo<CR>", { buffer=bufnr, desc="Info" })
            vim.keymap.set('n', '<leader>lr', "<Cmd>LspRestart<CR>", { buffer=bufnr, desc="Restart" })
            vim.keymap.set('n', '<leader>lk', "<Cmd>LspStop<CR>", { buffer=bufnr, desc="Stop" })
            vim.keymap.set('n', '<leader>ls', "<Cmd>LspStart<CR>", { buffer=bufnr, desc="Start" })
            vim.keymap.set('n', '<leader>ll', "<Cmd>LspLog<CR>", { buffer=bufnr, desc="Log" })
            vim.keymap.set('n', '<leader>lc', "<Cmd>CmpStatus<CR>", { buffer=bufnr, desc="Cmp status" })

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
                vim.keymap.set('n', 'gh', vim.lsp.buf.signature_help, { buffer = bufnr, desc = "Signature" })
                vim.keymap.set('n', '<leader>wa', vim.lsp.buf.add_workspace_folder, { buffer = bufnr, desc = "Add" })
                vim.keymap.set('n', '<leader>wr', vim.lsp.buf.remove_workspace_folder, { buffer = bufnr, desc = "Remove" })
                vim.keymap.set('n', '<leader>wl', function() print(vim.inspect(vim.lsp.buf.list_workspace_folders())) end, { buffer = bufnr, desc = "List" })
                vim.keymap.set('n', '<leader>dt', vim.lsp.buf.type_definition, { buffer = bufnr, desc = "Type def" })
                vim.keymap.set('n', '<leader>rn', vim.lsp.buf.rename, { buffer = bufnr, desc = "Rename" })
                vim.keymap.set('n', '<leader>ca', vim.lsp.buf.code_action, { buffer = bufnr, desc = "Code action" })
                -- vim.keymap.set('n', 'gr', vim.lsp.buf.references, { buffer=bufnr, desc="Goto references" })
                vim.keymap.set('n', 'gr', "<cmd>Telescope lsp_references<CR>", { buffer = bufnr, desc = "Goto references" })
                vim.keymap.set({'n', 'v'}, '<leader>lf', vim.lsp.buf.format, { buffer = bufnr, desc = "Format" })
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
            -- lsp.jedi_language_server.setup {}
            -- https://old.reddit.com/r/neovim/comments/1bh0kba/psa_new_python_lsp_that_supports_inlay_hints_and/
            -- https://github.com/neovim/nvim-lspconfig/blob/master/doc/server_configurations.md#basedpyright
            lsp.basedpyright.setup {
                settings = {
                    python = {
                        -- Has :PyrightSetPythonPath to set it on the fly
                        -- pythonPath = "~/bin/mambaforge/bin/python",
                    },
                    basedpyright = {
                        analysis = {
                            -- defaults to complaining about unknown types, and we don't want to be reminded to specify types.
                            -- Plus when using other's code that we can't change there will also be warnings about their lack of type declaration.
                            -- https://detachhead.github.io/basedpyright/#/configuration
                            typeCheckingMode = "standard",
                        }
                    }
                }
            }

            lsp.julials.setup {
                -- cmd = {
                --     "julia", "--startup-file=no", "--history-file=no",
                --     "-J", "/Users/cdmadsen/.julia/environments/nvim-lspconfig/languageserver.dylib",
                --     "-e",
                --     '# Load LanguageServer.jl: attempt to load from ~/.julia/environments/nvim-lspconfig\
                --     # with the regular load path as a fallback\
                --     ls_install_path = joinpath(\
                --     get(DEPOT_PATH, 1, joinpath(homedir(), ".julia")),\
                --     "environments", "nvim-lspconfig"\
                --     )\
                --     pushfirst!(LOAD_PATH, ls_install_path)\
                --     using LanguageServer\
                --     popfirst!(LOAD_PATH)\
                --     depot_path = get(ENV, "JULIA_DEPOT_PATH", "")\
                --     project_path = let\
                --     dirname(something(\
                --     ## 1. Finds an explicitly set project (JULIA_PROJECT)\
                --     Base.load_path_expand((\
                --     p = get(ENV, "JULIA_PROJECT", nothing);\
                --     p === nothing ? nothing : isempty(p) ? nothing : p\
                --     )),\
                --     ## 2. Look for a Project.toml file in the current working directory,\
                --     ##    or parent directories, with $HOME as an upper boundary\
                --     Base.current_project(),\
                --     ## 3. First entry in the load path\
                --     get(Base.load_path(), 1, nothing),\
                --     ## 4. Fallback to default global environment,\
                --     ##    this is more or less unreachable\
                --     Base.load_path_expand("@v#.#"),\
                --     ))\
                --     end\
                --     @info "Running language server" VERSION pwd() project_path depot_path\
                --     server = LanguageServer.LanguageServerInstance(stdin, stdout, project_path, depot_path)\
                --     server.runlinter = true\
                --     run(server)',
                -- },
                -- on_new_config = function(new_config, _)
                --     local julia = vim.fn.expand("~/.julia/environments/nvim-lspconfig/bin/julia")
                --     local dylib = vim.fn.expand("~/.julia/environments/nvim-lspconfig/languageserver.dylib")
                --     -- check if we have made the dedicated julia env which should have a custom system image
                --     -- if require 'lspconfig'.util.path.is_file(julia) then
                --         new_config.cmd[1] = julia
                --         -- new_config.cmd = {julia, "-J", dylib}
                --     -- end
                -- end
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
            lsp.matlab_ls.setup {}
            lsp.marksman.setup { filetypes = { "markdown" } }
            -- https://github.com/neovim/nvim-lspconfig/blob/master/doc/server_configurations.md#texlab
            lsp.texlab.setup {
                settings = { texlab = { experimental = {
                    -- custom citation function completion.
                    -- In some projects I use \citea{...} as \citeauthor{...}~\cite{...}
                    citationCommands = {"citea"}
                }}}
            }

            lsp.csharp_ls.setup {
                -- AutomaticWorkspaceInit = true,
            }

            lsp.rust_analyzer.setup {}

            -- doesn't seem to do anything
            lsp.wgsl_analyzer.setup {}

            -- MasonInstall typescript-language-server
            -- Used on javascript as well.
            lsp.tsserver.setup {}

            lsp.yamlls.setup {
                settings = { yaml = { schemas = {
                    [vim.opt.runtimepath:get()[1] .. "/lua/completion/asciidoc-theme.json"] = "*.adoc.yml",
                }}}
            }

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

            -- https://github.com/Myriad-Dreamin/tinymist/blob/main/editors/neovim/Configuration.md
            lsp.tinymist.setup {
                on_attach = function (client, bufnr)
                    -- Don't override default on_attach
                    lsp.util.default_config.on_attach(client, bufnr)

                    local function pinMain(fname)
                        return vim.lsp.buf.execute_command({ command = 'tinymist.pinMain', arguments = { fname } })
                    end
                    vim.keymap.set('n', '<leader><leader>p', function ()
                        return pinMain(vim.api.nvim_buf_get_name(0))
                    end, { desc="Pin buffer as main" })
                    -- search upwards for a main.typ
                    local mainfile = vim.fs.find("main.typ", {type="file", upward=true})[1]
                    if mainfile ~= nil then
                        vim.lsp.buf.execute_command({ command = 'tinymist.pinMain', arguments = { mainfile } })
                    end
                end
            }

            lsp.kotlin_language_server.setup {
                -- assume project root is at git root.
                -- Some imports don't work if unset.
                -- If you have a project without git root at project root then 
                -- search upwards for something else with this function.
                root_dir = function() return vim.env.ROOT end
            }

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
                -- "csharp_ls", -- instead see install/csharp.sh. There's also omnisharp on other LSPs on Mason.
                -- python:
                -- "jedi_language_server",
                -- "pyright",
                -- "pylsp",
                "basedpyright",

                "julials",
                -- "ltex", -- grammar check for latex, markdown, etc
                "texlab",
                -- "r_language_server",
                "lua_ls",
                -- "vimls",
                "rust_analyzer",
                "tsserver",     -- javascript. MasonInstall typescript-language-server
                "sqlls",
                -- "latexindent",
                "matlab_ls",
                -- "kotlin-language-server",
                "yamlls",
                "typst_lsp",
                -- "awk_ls",
                -- "cypher_ls", -- neo4j
                -- spell check/autocorrectors:
                -- "misspell",
                -- "typos",
                -- "typos_lsp",
            }
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
        config = function ()
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
        config = function ()
            local crates = require"crates"
            crates.setup {
                completion = {
                    cmp = {enabled = false},
                    crates = {
                        enabled = true,
                        max_results = 8, -- The maximum number of search results to display
                        min_chars = 3, -- The minimum number of charaters to type before completions begin appearing
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

            vim.keymap.set("n", "<leader><leader>t", crates.toggle, {desc="Toggle crates"})
            vim.keymap.set("n", "<leader><leader>r", crates.reload, {desc="Reload crates"})

            vim.keymap.set("n", "<leader><leader>v", crates.show_versions_popup, {desc="Show versions"})
            vim.keymap.set("n", "<leader><leader>f", crates.show_features_popup, {desc="Show features"})
            vim.keymap.set("n", "<leader><leader>d", crates.show_dependencies_popup, {desc="Show dependencies"})

            vim.keymap.set("n", "<leader><leader>u", crates.update_crate, {desc="Update"})
            vim.keymap.set("v", "<leader><leader>u", crates.update_crates, {desc="Update"})
            vim.keymap.set("n", "<leader><leader>a", crates.update_all_crates, {desc="Update all"})
            vim.keymap.set("n", "<leader><leader>U", crates.upgrade_crate, {desc="Upgrade"})
            vim.keymap.set("v", "<leader><leader>U", crates.upgrade_crates, {desc="Upgrade"})
            vim.keymap.set("n", "<leader><leader>A", crates.upgrade_all_crates, {desc="Upgrade all"})

            vim.keymap.set("n", "<leader><leader>x", crates.expand_plain_crate_to_inline_table, {desc="Plain crate -> inline table"})
            vim.keymap.set("n", "<leader><leader>X", crates.extract_crate_into_table, {desc="Crate -> table"})

            vim.keymap.set("n", "<leader><leader>H", crates.open_homepage, {desc="Homepage"})
            vim.keymap.set("n", "<leader><leader>R", crates.open_repository, {desc="Repo"})
            vim.keymap.set("n", "<leader><leader>D", crates.open_documentation, {desc="Documentation"})
            vim.keymap.set("n", "<leader><leader>C", crates.open_crates_io, {desc="crates.io"})
            vim.keymap.set("n", "<leader><leader>L", crates.open_lib_rs, {desc="lib.rs"})
        end,
    },

}
