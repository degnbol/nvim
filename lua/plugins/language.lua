#!/usr/bin/env lua
local hl = require "utils/highlights"

return {
    -- add keybindings to toggle comments with motions etc.
    {"terrortylor/nvim-comment", config=function() require'nvim_comment'.setup {
        -- default is ic which should instead be the treesitter @comment.inner.
        -- aC is almost the same and provided by a vim plugin that extends the kana plugin.
        comment_chunk_text_object="iC",
    } end},
    -- has the useful gcO and gcA extra mappings, but the basic mappings aren't working as I like from terrortylor
    {"numToStr/Comment.nvim", opts={ mappings={ basic=false, } }},
    -- julia support, colors and unicode substitution. CANNOT use ft=julia
    {"JuliaEditorSupport/julia-vim",
        config=function()
            -- this was necessary, random lhs rhs messages was appearing 
            vim.g.latex_to_unicode_tab = "off"
            -- auto insert latex if moving on with e.g. space or other writing that is not part of a unicode char
            -- Super useful for e.g. completing \notin if LSP can't be used for completing it.
            vim.g.latex_to_unicode_auto = true
    end},
    -- python aware changes to [], [[, ]], ][, ]m, ]M, [m, [M for moving cursor to starts and ends of python functions. This should be covered by tree sitter in the future when they add support for visual mode
    {"jeetsukumaran/vim-pythonsense", ft='python'},
    -- try it out with :Cheat <query> where the query should be search terms like you would search in StackOverflow for answers
    {"RishabhRD/nvim-cheat.sh", config=function()
        -- default is float i.e. a floating window
        -- vim.g.cheat_default_window_layout = 'tab'
    end, dependencies={"RishabhRD/popfix"}},
    -- {"mrjones2014/dash.nvim", build='make install', dependencies='nvim-telescope/telescope.nvim'}, -- :DashWord with <leader>K. conf in telescope-conf.lua
    -- "tpope/vim-sleuth", -- sleuth that let's you autodetect if file is using 2 or 4 spaces. Mistakenly set noexpandtab
    -- {"preservim/vim-markdown", dependencies={"godlygeek/tabular"}}, -- conceal markdown expressions like _emphasis_ and folding. Overkill, see {after/,}syntax/markdown.vim
    -- :MarkdownPreview live in browser
    {"iamcco/markdown-preview.nvim", build=':call mkdp#util#install()', ft='markdown'},
    -- for asciidoctor filetype use ".adoc"
    -- adoc argues to be a better replacement for markdown in general
    -- https://docs.asciidoctor.org/asciidoc/latest/asciidoc-vs-markdown/
    -- rst seems instead to be all about technical documentation.
    {
        -- a fork with rename asciidoctor -> asciidoc
        "mkschreder/vim-asciidoc",
        branch = "bugfix/asciidoctor",
        -- doesn't seem there is any asciidoc plugin that helps with goto definition
        ft={'asciidoc', 'chordpro'},
        config=function()
            -- conceal _ and * in _italic_ and *bold*
            vim.g.asciidoctor_syntax_conceal = 1
            -- If attribute `:pdf-theme: book` is written then with the 
            -- following setting it will look for "book-theme.yml" in the same 
            -- directory. Default is some themes folder buried in the 
            -- installation for asciidoctor.
            vim.g.asciidoctor_pdf_themes_path = '.'
            -- somehow setting colorscheme in way that resets some highlight groups, 
            -- so we set them again here:
            local grp = vim.api.nvim_create_augroup("adoc", {clear=true})
            vim.api.nvim_create_autocmd("Colorscheme", {
                buffer = 0,
                group = grp,
                callback = function ()
                    hl.set("asciidoctorBold", {bold=true})
                    hl.set("asciidoctorItalic", {italic=true})
                    hl.set("asciidoctorBoldItalic", {bold=true, italic=true})
                    hl.set("asciidoctorBoldComment", {bold=true, fg=hl.get("Comment")['fg']})
                    hl.set("asciidoctorItalicComment", {italic=true, fg=hl.get("Comment")['fg']})
                    hl.set("asciidoctorBoldItalicComment", {bold=true, italic=true, fg=hl.get("Comment")['fg']})
                    hl.link("asciidoctorTitleDelimiter", "Comment")
                    hl.link("asciidocPassthrough", "Constant") -- same as asciidoctorCode
                    hl.rev("asciidocHighlight")
                    hl.set("asciidocUnderline", {underline=true})
                    hl.set("asciidocBoldUnderline", {underline=true, bold=true})
                    hl.set("asciidocItalicUnderline", {underline=true, italic=true})
                    hl.set("asciidocStrikethrough", {strikethrough=true})
                    hl.set("asciidocBoldStrikethrough", {strikethrough=true, bold=true})
                    hl.set("asciidocItalicStrikethrough", {strikethrough=true, italic=true})
                    for i = 1, 6 do
                        hl.link("asciidoctorH"..i.."Delimiter", "Comment")
                    end
                    -- there are more hi groups that might be unset
                    -- ...
                    -- hide comment delim (hl group set in syntax file).
                    hl.fg("commentDelimiter", hl.get("Normal")["bg"])
                    hl.link("filenameCommentNoSpell", "Comment")
                    hl.link("UrlCommentNoSpell", "Comment")
                    hl.link("linebreak", "Comment")
                    hl.set("Geo", {underline=true})
                    hl.set("Chord", {underline=true})
                end
            })
            -- use conversion to PDF as default for :make
            vim.cmd [[compiler asciidoctor2pdf]]
            -- double <CR> to auto-close after successful compilation.
            -- If this is not desired then use :make.
            vim.keymap.set('n', '<leader>cc', "<Cmd>Asciidoctor2PDF<CR><CR>", { buffer=true, desc="Compile to PDF" })
            vim.keymap.set('n', '<leader>oo', "<Cmd>AsciidoctorOpenPDF<CR><CR>", { buffer=true, desc="Open compiled PDF" })
            -- start autocompiling on save
            vim.keymap.set('n', '<leader>cC', function ()
                local grp = vim.api.nvim_create_augroup("asciidocCompile", {clear=true})
                vim.api.nvim_create_autocmd("BufWritePost", {
                    buffer = 0,
                    group = grp,
                    command = "silent Asciidoctor2PDF"
                })
            end, { buffer=true, desc="Compile on save" })
        end,
    },
    -- https://quarto.org/
    {"quarto-dev/quarto-nvim", ft="quarto"},
    -- flashing for code blocks
    {"lukas-reineke/headlines.nvim", enabled = false, ft = {'markdown', 'norg', 'orgmode', 'rst', 'asciidoc', 'asciidoctor'},
    dependencies = { "nvim-treesitter/nvim-treesitter" },
    opts = { }, },
    -- "elzr/vim-json", -- json
    {"OmniSharp/omnisharp-vim", ft="cs"},

    -- autoclose pairs.
    -- "m4xshen/autoclose.nvim" is too simple.
    -- "windwp/nvim-autopairs" doesn't delete properly even with the check_ts 
    -- (treesitter) setting on. I also tried pears.nvim. I haven't tried mini.pairs
    { "tmsvg/pear-tree", enabled = false, config=function ()
        -- vim.g.pear_tree_smart_openers = 1
        -- vim.g.pear_tree_smart_closers = 1 # bad in lua, type "{}" -> "{}}"
        vim.g.pear_tree_smart_backspace = 1
        -- uncomment to not hide the closing bracket on newline at the cost of 
        -- dot-repeat only performing the last part of the edit.
        -- vim.g.pear_tree_repeatable_expand = 0
        vim.g.pear_tree_ft_disabled = { "TelescopePrompt", "NvimTree", "qf", "tex", }
    end},
    {
        "nvim-neorg/neorg",
        -- Seems to throw errors at random times and I don't see pandoc 
        -- conversion support so I think asciidoctor might cover the use 
        -- case better.
        ft = "norg",
        opts = {
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
                ["core.completion"] = {config={
                    engine="nvim-cmp",
                    -- a notebook icon to indicate neorg as completion source
                    name="",
                }},
                ["core.export"] = {},
            }
        },
    },
    {
        "kaarmu/typst.vim",
        ft = "typst",
        init = function ()
            vim.g.typst_pdf_viewer = "sioyek"
        end,
    },
    -- "MrPicklePinosaur/typst-conceal.vim",

    "fladson/vim-kitty", -- syntax highlights for kitty conf

    -- basic kotlin support
    "udalov/kotlin-vim",

    {
        -- pretty decorations on help.
        -- Currently a bit buggy, might use later when you can toggle it without getting error.
        "OXY2DEV/helpview.nvim",
        enabled = false,
        lazy = false, -- Recommended
        dependencies = {
            "nvim-treesitter/nvim-treesitter"
        }
    },
    "mityu/vim-applescript",
}

