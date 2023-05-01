#!/usr/bin/env lua
local g = vim.g
return {
    -- add keybindings to toggle comments with motions etc.
    {"terrortylor/nvim-comment", config=function() require'nvim_comment'.setup() end},
    -- has the useful gcO and gcA extra mappings, but the basic mappings aren't working as I like from terrortylor
    {"numToStr/Comment.nvim", opts={ mappings={ basic=false, } }},
    -- "windwp/nvim-autopairs", -- auto add second parenthesis etc.
    -- julia support, colors and unicode substitution. CANNOT use ft=julia
    {"JuliaEditorSupport/julia-vim", config=function()
        -- this was necessary, random lhs rhs messages was appearing 
        g.latex_to_unicode_tab = "off"
        -- auto insert latex if moving on with e.g. space or other writing that is not part of a unicode char
        -- g.latex_to_unicode_auto = true
    end},
    {"jeetsukumaran/vim-pythonsense", ft='python'}, -- python aware changes to [], [[, ]], ][, ]m, ]M, [m, [M for moving cursor to starts and ends of python functions. This should be covered by tree sitter in the future when they add support for visual mode
    -- try it out with :Cheat <query> where the query should be search terms like you would search in StackOverflow for answers
    {"RishabhRD/nvim-cheat.sh", config=function()
        -- default is float i.e. a floating window
        -- vim.g.cheat_default_window_layout = 'tab'
    end, dependencies={"RishabhRD/popfix"}},
    -- {"mrjones2014/dash.nvim", build='make install', dependencies='nvim-telescope/telescope.nvim'}, -- :DashWord with <leader>K. conf in telescope-conf.lua
    -- "tpope/vim-sleuth", -- sleuth that let's you autodetect if file is using 2 or 4 spaces. Mistakenly set noexpandtab
    "tpope/vim-abolish", -- crs: snake_case, crm: MixedCase, crc: camelCase, cru: UPPER_CASE, cr-: dash-case, cr.: dot.case, cr<SPACE>: space case, crt: Title Case
    -- {"preservim/vim-markdown", dependencies={"godlygeek/tabular"}}, -- conceal markdown expressions like _emphasis_ and folding. Overkill, see {after/,}syntax/markdown.vim
    -- :MarkdownPreview live in browser
    {"iamcco/markdown-preview.nvim", build=':call mkdp#util#install()', ft='markdown'},
    {"habamax/vim-asciidoctor", config=function()
        -- conceal _ and * in _italic_ and *bold*
        vim.g.asciidoctor_syntax_conceal = 1
    end, ft='asciidoctor'},
    -- https://quarto.org/
    {"quarto-dev/quarto-vim", dependencies={"vim-pandoc/vim-pandoc-syntax"}, ft="quarto"},
    {"habamax/vim-rst"},
    "jbyuki/nabla.nvim", -- show pretty math in term
    -- "elzr/vim-json", -- json
    {"OmniSharp/omnisharp-vim", ft="cs"},
}

