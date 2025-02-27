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
        config = function()
            -- default config copied from https://github.com/neovim/nvim-lspconfig
            -- inspiration from https://vonheikemen.github.io/devlog/tools/setup-nvim-lspconfig-plus-nvim-cmp/

            -- replace the default lsp diagnostic letters with prettier symbols
            vim.fn.sign_define("LspDiagnosticsSignError", { text = "", numhl = "LspDiagnosticsDefaultError" })
            vim.fn.sign_define("LspDiagnosticsSignWarning", { text = "", numhl = "LspDiagnosticsDefaultWarning" })
            vim.fn.sign_define("LspDiagnosticsSignInformation", { text = "", numhl = "LspDiagnosticsDefaultInformation" })
            vim.fn.sign_define("LspDiagnosticsSignHint", { text = "", numhl = "LspDiagnosticsDefaultHint" })

            -- LSP and completion status, overall conf etc.
            vim.keymap.set('n', '<leader>li', "<Cmd>LspInfo<CR>", { desc="Info" })
            -- can't use backspace since it is hardcoded by mini.clue for up one level
            vim.keymap.set('n', '<leader>l<Del>', "<Cmd>LspStop<CR>", { desc="Stop" })
            vim.keymap.set('n', '<leader>l1', "<Cmd>LspStart<CR>", { desc="Start" })
            vim.keymap.set('n', '<leader>l!', "<Cmd>LspRestart<CR>", { desc="Restart" })
            vim.keymap.set('n', '<leader>lL', "<Cmd>LspLog<CR>", { desc="Log" })
            -- match other completion related entries under x
            vim.keymap.set('n', '<leader>xS', "<Cmd>CmpStatus<CR>", { desc="Cmp status" })

            -- See `:help vim.diagnostic.*` for documentation on any of the below functions
            vim.keymap.set('n', '<leader>dd', vim.diagnostic.open_float, { desc = "Line diagnostic" })
            -- can't use backspace since it is hardcoded by mini.clue for up one level
            vim.keymap.set('n', '<leader>d<Del>', function() vim.diagnostic.enable(false) end, { desc = "Disable diagnostics" })
            vim.keymap.set('n', '<leader>d1', vim.diagnostic.enable, { desc = "Enable diagnostics" })
            vim.keymap.set('n', '[d', vim.diagnostic.goto_prev, { desc = "Diagnostic" })
            vim.keymap.set('n', ']d', vim.diagnostic.goto_next, { desc = "Diagnostic" })
            vim.keymap.set('n', '<leader>dl', vim.diagnostic.setloclist, { desc = "Loclist diagnostics" })

            -- Use an on_attach function to only map the following keys
            -- after the language server attaches to the current buffer
            local on_attach = function(client, bufnr)
                local function map(desc, keys, func, mode)
                    mode = mode or 'n'
                    vim.keymap.set(mode, keys, func, { buffer = bufnr, desc = desc })
                end
                local function mapfzf(desc, keys, funcname, mode)
                    map(desc, keys, function ()
                        require"fzf-lua"["lsp_" .. funcname]()
                    end, mode)
                end
                -- See `:help vim.lsp.*` for documentation on any of the below functions
                -- TODO: have treesitter fallback for things like go to references for when LSP is not attached.
                -- `vim.lsp.buf.references`
                mapfzf("Goto references", "gr", "references")
                -- `vim.lsp.buf.definition`
                mapfzf("Goto definition", "gd", "definitions")
                -- `vim.lsp.buf.declaration`
                mapfzf("Goto declaration", "gD", "declarations")
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
                vim.keymap.set({'n', 'v'}, '<leader>lf', vim.lsp.buf.format, { buffer = bufnr, desc = "Format" })
            end

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

            -- https://github.com/neovim/nvim-lspconfig/blob/master/doc/configs.md#julials
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
            -- Don't autoshow completion in cmdline
            -- https://github.com/neovim/nvim-lspconfig/blob/master/doc/server_configurations.md#texlab
            lsp.texlab.setup {
                -- https://github.com/latex-lsp/texlab/wiki/Configuration
                settings = { texlab = {
                    experimental = {
                        -- custom citation function completion.
                        -- In some projects I use \citea{...} as \citeauthor{...}~\cite{...} and \citePDBlit etc.
                        citationCommands = {"citea", "citePDB", "citePDBlit"},
                        -- \subref from the subpcation package.
                        -- \footref from the footmisc package to refer to footnote a second (or more) time.
                        -- In some projects I use \see[...]{...} as (see~\cref{...}...)
                        labelReferenceCommands = {"footref", "subref", "see", "seename", "seefull"},
                    },
                }}
            }

            lsp.csharp_ls.setup {
                -- AutomaticWorkspaceInit = true,
            }

            lsp.rust_analyzer.setup {}

            -- doesn't seem to do anything
            lsp.wgsl_analyzer.setup {}

            -- MasonInstall typescript-language-server
            -- Used on javascript as well.
            lsp.ts_ls.setup {}

            lsp.yamlls.setup {
                settings = { yaml = { schemas = {
                    [vim.opt.runtimepath:get()[1] .. "/lua/completion/asciidoc/asciidoc-theme.json"] = "*.adoc.yml",
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
            -- for typst
            lsp.tinymist.setup {
                on_attach = function (client, bufnr)
                    -- Don't override default on_attach
                    lsp.util.default_config.on_attach(client, bufnr)

                    local function pinMain(fname)
                        return vim.lsp.buf.execute_command({ command = 'tinymist.pinMain', arguments = { fname } })
                    end
                    vim.keymap.set('n', '<LocalLeader>p', function ()
                        return pinMain(vim.api.nvim_buf_get_name(0))
                    end, { buffer=true, desc="Pin buffer as main" })
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
        lazy = true, -- load as mason-lspconfig dep
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
                "ts_ls",     -- javascript. MasonInstall typescript-language-server
                "sqlls",
                -- "latexindent",
                "matlab_ls",
                -- "kotlin-language-server",
                "yamlls",
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

            local function localmap(key, func, desc)
                vim.keymap.set("n", "<LocalLeader>" .. key, crates[func], {desc=desc})
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
