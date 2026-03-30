local util = require "utils/init"
local hi = require "utils/highlights"

return {
    -- add keybindings to toggle comments with motions etc.
    {
        "nvim-comment",
        after = function()
            require 'nvim_comment'.setup {
                -- default is ic which should instead be the treesitter @comment.inner.
                -- aC is almost the same and provided by a vim plugin that extends the kana plugin.
                comment_chunk_text_object = "iC",
            }
        end
    },
    -- julia support, colors and unicode substitution. CANNOT use ft=julia
    {
        "julia-vim",
        -- Let's see if we actually need this. It's not well maintained and is
        -- really agressively enforcing things like indentexpr that doesn't
        -- work at all. We have colours from treesitter and can write our own otherwise.
        -- We get unicode from lsp, e.g. \notin<C-space>
        enabled = false,
        after = function()
            -- this was necessary, random lhs rhs messages was appearing
            vim.g.latex_to_unicode_tab = "off"
            -- auto insert latex if moving on with e.g. space or other writing that is not part of a unicode char
            -- Super useful for e.g. completing \notin if LSP can't be used for completing it.
            vim.g.latex_to_unicode_auto = true
        end
    },
    -- python aware changes to [], [[, ]], ][, ]m, ]M, [m, [M for moving cursor to starts and ends of python functions. This should be covered by tree sitter in the future when they add support for visual mode
    { "vim-pythonsense", ft = 'python' },
    -- try it out with :Cheat <query> where the query should be search terms like you would search in StackOverflow for answers
    {
        "nvim-cheat.sh",
        after = function()
            -- default is float i.e. a floating window
            -- vim.g.cheat_default_window_layout = 'tab'
        end,
    },
    -- {"mrjones2014/dash.nvim", build='make install', dependencies='nvim-telescope/telescope.nvim'}, -- :DashWord with <leader>K. conf in telescope-conf.lua
    -- "tpope/vim-sleuth", -- sleuth that let's you autodetect if file is using 2 or 4 spaces. Mistakenly set noexpandtab
    -- {"preservim/vim-markdown", dependencies={"godlygeek/tabular"}}, -- conceal markdown expressions like _emphasis_ and folding. Overkill, see {after/,}syntax/markdown.vim
    -- :MarkdownPreview live in browser
    { "markdown-preview.nvim", ft = 'markdown' },
    -- for asciidoctor filetype use ".adoc"
    -- adoc argues to be a better replacement for markdown in general
    -- https://docs.asciidoctor.org/asciidoc/latest/asciidoc-vs-markdown/
    -- rst seems instead to be all about technical documentation.
    {
        -- a fork with rename asciidoctor -> asciidoc
        "vim-asciidoc",
        -- doesn't seem there is any asciidoc plugin that helps with goto definition
        ft = { 'asciidoc', 'chordpro' },
        after = function()
            -- conceal _ and * in _italic_ and *bold*
            vim.g.asciidoctor_syntax_conceal = 1
            -- If attribute `:pdf-theme: book` is written then with the
            -- following setting it will look for "book-theme.yml" in the same
            -- directory. Default is some themes folder buried in the
            -- installation for asciidoctor.
            vim.g.asciidoctor_pdf_themes_path = '.'
            -- somehow setting colorscheme in way that resets some highlight groups,
            -- so we set them again here:
            local grp = vim.api.nvim_create_augroup("adoc", { clear = true })
            vim.api.nvim_create_autocmd("ColorScheme", {
                buffer = 0,
                group = grp,
                callback = function()
                    hi.set("asciidoctorBold", { bold = true })
                    hi.set("asciidoctorItalic", { italic = true })
                    hi.set("asciidoctorBoldItalic", { bold = true, italic = true })
                    hi.set("asciidoctorBoldComment", { bold = true, fg = hi.get("Comment")['fg'] })
                    hi.set("asciidoctorItalicComment", { italic = true, fg = hi.get("Comment")['fg'] })
                    hi.set("asciidoctorBoldItalicComment", { bold = true, italic = true, fg = hi.get("Comment")['fg'] })
                    hi.link("asciidoctorTitleDelimiter", "Comment")
                    hi.link("asciidocPassthrough", "Constant") -- same as asciidoctorCode
                    hi.rev("asciidocHighlight")
                    hi.set("asciidocUnderline", { underline = true })
                    hi.set("asciidocBoldUnderline", { underline = true, bold = true })
                    hi.set("asciidocItalicUnderline", { underline = true, italic = true })
                    hi.set("asciidocStrikethrough", { strikethrough = true })
                    hi.set("asciidocBoldStrikethrough", { strikethrough = true, bold = true })
                    hi.set("asciidocItalicStrikethrough", { strikethrough = true, italic = true })
                    for i = 1, 6 do
                        hi.link("asciidoctorH" .. i .. "Delimiter", "Comment")
                    end
                    -- there are more hi groups that might be unset
                    -- ...
                    -- hide comment delim (hl group set in syntax file).
                    hi.setfg("commentDelimiter", hi.get("Normal")["bg"])
                    hi.link("filenameCommentNoSpell", "Comment")
                    hi.link("UrlCommentNoSpell", "Comment")
                    hi.link("linebreak", "Comment")
                    hi.set("Geo", { underline = true })
                    hi.set("Chord", { underline = true })
                end
            })
            -- use conversion to PDF as default for :make
            vim.cmd [[compiler asciidoctor2pdf]]
            -- double <CR> to auto-close after successful compilation.
            -- If this is not desired then use :make.
            vim.keymap.set('n', '<leader>cc', "<Cmd>Asciidoctor2PDF<CR><CR>", { buffer = true, desc = "Compile to PDF" })
            vim.keymap.set('n', '<leader>oo', "<Cmd>AsciidoctorOpenPDF<CR><CR>",
                { buffer = true, desc = "Open compiled PDF" })
            -- start autocompiling on save
            vim.keymap.set('n', '<leader>cC', function()
                vim.api.nvim_create_autocmd("BufWritePost", {
                    buffer = 0,
                    group = vim.api.nvim_create_augroup("asciidocCompile", { clear = true }),
                    command = "silent Asciidoctor2PDF"
                })
            end, { buffer = true, desc = "Compile on save" })
        end,
    },
    -- https://quarto.org/
    {
        "quarto-nvim",
        ft = "quarto",
    },
    -- flashing for code blocks
    {
        "headlines.nvim",
        enabled = false,
        ft = { 'markdown', 'norg', 'orgmode', 'rst', 'asciidoc', 'asciidoctor' },
        after = function()
            require("headlines").setup {}
        end,
    },
    -- "elzr/vim-json", -- json
    { "omnisharp-vim", ft = "cs" },
    {
        "neorg",
        -- Seems to throw errors at random times and I don't see pandoc
        -- conversion support so I think asciidoctor might cover the use
        -- case better.
        ft = "norg",
        after = function()
            require("neorg").setup {
                load = {
                    ["core.defaults"] = {},
                    ["core.dirman"] = {
                        config = {
                            workspaces = {
                                bio = "~/Documents/Bio/",
                                writing = "~/Documents/Writing/",
                            }
                        }
                    },
                    ["core.concealer"] = {},
                    ["core.completion"] = {
                        config = {
                            engine = "nvim-cmp",
                            -- a notebook icon to indicate neorg as completion source
                            name = "",
                        }
                    },
                    ["core.export"] = {},
                }
            }
        end,
    },
    {
        "typst.vim",
        ft = "typst",
        before = function()
            if util.is_mac() then
                vim.g.typst_pdf_viewer = "skim"
            else
                vim.g.typst_pdf_viewer = "zathura"
            end
        end,
    },
    -- "MrPicklePinosaur/typst-conceal.vim",

    -- syntax highlights for kitty conf
    {
        "vim-kitty",
        ft = "kitty",
        after = function ()
            local grp = vim.api.nvim_create_augroup("kitty", { clear = true })
            vim.api.nvim_create_autocmd("ColorScheme", {
                buffer = 0,
                group = grp,
                callback = function()
                    -- Italic is for the small set of rarely used keywords in most
                    -- languages, not for the large amount of "keywords" in some
                    -- language like the kitty.conf.
                    -- There's also kittyKW, why do both exist? They're both from this plugin.
                    hi.set("kittyKeyword", { italic = false, fg=hi.fg("Keyword") })
                    -- Maybe fine to leave italic?
                    -- hi.link("kittyMap", "kittyKeyword")
                    -- Was linked to Keyword
                    hi.link("kittyInclude", "Include")
                end
            })
        end,
    },

    -- basic kotlin support
    { "kotlin-vim" },

    {
        -- pretty decorations on help.
        -- Currently a bit buggy, might use later when you can toggle it without getting error.
        "helpview.nvim",
        enabled = false,
        lazy = false, -- Recommended
    },
    { "vim-applescript" },
    -- Consistent completion error in sql without this.
    {
        "dbext.vim",
        ft = "sql",
    },
    {
        "srt.nvim",
        after = function()
            require "srtnvim".setup {}
        end,
    },
    {
        "math-conceal.nvim",
        -- Project was renamed to reflect expanding to also cover typst.
        -- Code in lua/ etc for project still has old name but this might get fixed quickly.
        event = "DeferredUIEnter",
        enabled = false,
        --- @type LaTeXConcealOptions
        after = function()
            require("latex-conceal").setup {
                enabled = true,
                conceal = {
                    "greek",
                    "script",
                    "math",
                    "font",
                    "delim",
                    "phy",
                },
                ft = { "tex", "latex", "markdown", "typst" },
            }
        end,
    },
    -- edit jupyter notebook. NOTE: requires `pip install jupytext`
    {
        "jupytext.vim",
        after = function()
            -- conversion settings
            vim.g.jupytext_fmt = "py:percent"
            -- highlight blocks.
            -- Doesn't work since treesitter trumps. Need to add a custom query.
            vim.cmd 'syn match Block /^# %%/'
            vim.api.nvim_set_hl(0, "Block", { link = "LineNr" })
        end,
    },
    -- {
    --     "salkin-mada/openscad.nvim",
    --     config = true,
    -- },
}
