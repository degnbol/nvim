#!/usr/bin/env lua
return {
    {"terrortylor/nvim-comment", config=function() require'nvim_comment'.setup() end}, -- add keybindings to toggle comments with motions etc. Consider alt: https://github.com/numToStr/Comment.nvim
    -- "windwp/nvim-autopairs", -- auto add second parenthesis etc.
    -- "lukas-reineke/indent-blankline.nvim", -- show "|" on indented lines
    {"JuliaEditorSupport/julia-vim", config=function() require'julia' end}, -- julia support, colors and unicode substitution. CANNOT use ft=julia
    {"jeetsukumaran/vim-pythonsense", ft='python'}, -- python aware changes to [], [[, ]], ][, ]m, ]M, [m, [M for moving cursor to starts and ends of python functions. This should be covered by tree sitter in the future when they add support for visual mode
    {"samirettali/shebang.nvim", config=function() require'shebang-nvim' end}, -- insert shebang on new file edit
    -- try it out with :Cheat <query> where the query should be search terms like you would search in StackOverflow for answers
    {"RishabhRD/nvim-cheat.sh", config=function() require'cheat' end, dependencies={"RishabhRD/popfix"}},
    -- {"mrjones2014/dash.nvim", build='make install', dependencies='nvim-telescope/telescope.nvim'}, -- :DashWord with <leader>K. conf in telescope-conf.lua
    -- "tpope/vim-sleuth", -- sleuth that let's you autodetect if file is using 2 or 4 spaces. Mistakenly set noexpandtab
    "tpope/vim-abolish", -- crs: snake_case, crm: MixedCase, crc: camelCase, cru: UPPER_CASE, cr-: dash-case, cr.: dot.case, cr<SPACE>: space case, crt: Title Case
    -- {"preservim/vim-markdown", dependencies={"godlygeek/tabular"}}, -- conceal markdown expressions like _emphasis_ and folding. Overkill, see {after/,}syntax/markdown.vim
    {"iamcco/markdown-preview.nvim", build=':call mkdp#util#install()', ft='markdown'}, -- :MarkdownPreview live in browser
    {"habamax/vim-asciidoctor", config=function() require'asciidoc' end, ft='asciidoctor'},
    {"quarto-dev/quarto-vim", dependencies={"vim-pandoc/vim-pandoc-syntax"}, ft="quarto"}, -- https://quarto.org/
    {"habamax/vim-rst"},
    "jbyuki/nabla.nvim", -- show pretty math in term
    {"goerz/jupytext.vim", config=function() require'jupytext-conf' end}, -- edit jupyter notebook. requires `pip install jupytext`
    -- "elzr/vim-json", -- json
    {"lervag/vimtex", config=function() require'vimtex' end}, -- :VimtexCompile. Adds so much more good stuff, e.g. dse, cse to delete or change surrounding env
    {"nvim-telescope/telescope-bibtex.nvim", dependencies={'nvim-telescope/telescope.nvim'}, config=function() require"telescope-bibtex-conf" end},
    -- bibliography references, mostly relevant for citations in .tex documents.
    -- { "jghauser/papis.nvim", 
    --     dependencies = { "kkharji/sqlite.lua", "nvim-lua/plenary.nvim", "MunifTanjim/nui.nvim", "nvim-treesitter/nvim-treesitter", "nvim-telescope/telescope.nvim", "hrsh7th/nvim-cmp" },
    --     rocks="lyaml", config = function() require("papis").setup() end,
    -- },
    {"OmniSharp/omnisharp-vim", ft="cs"},
    -- {"neoclide/coc.nvim", branch="release", build="npm install"},
}

