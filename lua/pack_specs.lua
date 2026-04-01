-- vim.pack.add registry — all remote plugins.
-- Dev plugins (modules/) are on rtp manually, not managed here.
local gh = function(x) return "https://github.com/" .. x end

vim.pack.add({
    -- lz.n (lazy-loader, bootstrapped via packadd before lz.n.load)
    gh("BirdeeHub/lz.n"),

    -- Core behaviour (init.lua)
    gh("tpope/vim-repeat"),
    gh("tpope/vim-abolish"),
    gh("andymass/vim-matchup"),
    gh("tpope/vim-unimpaired"),
    gh("tpope/vim-characterize"),
    gh("lewis6991/fileline.nvim"),
    gh("kana/vim-textobj-user"),
    gh("kana/vim-arpeggio"),
    gh("glts/vim-textobj-comment"),
    gh("chrisgrieser/nvim-various-textobjs"),
    gh("monaqa/dial.nvim"),
    gh("linrongbin16/gentags.nvim"),

    -- Icons
    gh("nvim-tree/nvim-web-devicons"),

    -- Alignment
    gh("junegunn/vim-easy-align"),
    gh("Vonr/align.nvim"),

    -- Colorschemes
    gh("EdenEast/nightfox.nvim"),
    gh("maxmx03/fluoromachine.nvim"),

    -- Clipboard
    gh("svermeulen/vim-subversive"),
    gh("svermeulen/vim-yoink"),
    gh("gbprod/cutlass.nvim"),
    gh("ibhagwan/smartyank.nvim"),
    gh("inkarkat/vim-ingo-library"),
    gh("inkarkat/vim-UnconditionalPaste"),
    gh("AndrewRadev/whitespaste.vim"),

    -- DAP
    gh("mfussenegger/nvim-dap"),
    gh("jay-babu/mason-nvim-dap.nvim"),
    gh("theHamsta/nvim-dap-virtual-text"),
    gh("rcarriga/nvim-dap-ui"),
    gh("nvim-neotest/nvim-nio"),

    -- Dashboard
    gh("nvimdev/dashboard-nvim"),

    -- Fidget
    { src = gh("j-hui/fidget.nvim"), version = "legacy" },

    -- Filetypes
    gh("terrortylor/nvim-comment"),
    gh("JuliaEditorSupport/julia-vim"),
    gh("jeetsukumaran/vim-pythonsense"),
    gh("RishabhRD/popfix"),
    gh("RishabhRD/nvim-cheat.sh"),
    gh("iamcco/markdown-preview.nvim"),
    { src = gh("mkschreder/vim-asciidoc"), version = "bugfix/asciidoctor" },
    gh("jmbuhr/otter.nvim"),
    gh("quarto-dev/quarto-nvim"),
    gh("lukas-reineke/headlines.nvim"),
    gh("OmniSharp/omnisharp-vim"),
    gh("nvim-neorg/neorg"),
    gh("kaarmu/typst.vim"),
    gh("fladson/vim-kitty"),
    gh("udalov/kotlin-vim"),
    gh("OXY2DEV/helpview.nvim"),
    gh("mityu/vim-applescript"),
    gh("vim-scripts/dbext.vim"),
    { src = gh("gert7/srt.nvim"), version = "main" },
    gh("pxwg/math-conceal.nvim"),
    gh("goerz/jupytext.vim"),

    -- Formatting
    gh("stevearc/conform.nvim"),

    -- Fuzzy finders
    gh("folke/snacks.nvim"),
    gh("ibhagwan/fzf-lua"),
    gh("nvim-telescope/telescope-fzf-native.nvim"),
    gh("nvim-telescope/telescope.nvim"),
    gh("sudormrfbin/cheatsheet.nvim"),
    gh("lalitmee/browse.nvim"),
    gh("nvim-telescope/telescope-bibtex.nvim"),
    gh("krissen/snacks-bibtex.nvim"),

    -- Git
    gh("tpope/vim-fugitive"),
    gh("NeogitOrg/neogit"),
    gh("rhysd/conflict-marker.vim"),
    gh("lewis6991/gitsigns.nvim"),
    gh("sindrets/diffview.nvim"),

    -- LSP
    gh("neovim/nvim-lspconfig"),
    gh("williamboman/mason.nvim"),
    gh("williamboman/mason-lspconfig.nvim"),
    gh("scalameta/nvim-metals"),
    gh("akinsho/flutter-tools.nvim"),
    gh("saecki/crates.nvim"),
    gh("stevearc/dressing.nvim"),

    -- Mini
    gh("echasnovski/mini.nvim"),

    -- Multicursor
    gh("mg979/vim-visual-multi"),
    gh("smoka7/hydra.nvim"),
    gh("smoka7/multicursors.nvim"),
    gh("jake-stewart/multicursor.nvim"),

    -- Notion
    gh("Al0den/notion.nvim"),

    -- Quickfix
    gh("stevearc/quicker.nvim"),
    gh("junegunn/fzf"),
    gh("kevinhwang91/nvim-bqf"),

    -- Search
    gh("google/vim-searchindex"),
    gh("kevinhwang91/nvim-hlslens"),
    gh("petertriho/nvim-scrollbar"),
    gh("haya14busa/vim-asterisk"),

    -- Snippets
    { src = gh("L3MON4D3/LuaSnip"), version = vim.version.range("2") },

    -- Splitjoin
    gh("Wansmer/treesj"),
    gh("AckslD/nvim-trevJ.lua"),
    gh("AndrewRadev/splitjoin.vim"),
    gh("sgur/vim-textobj-parameter"),
    gh("flwyd/vim-conjoin"),

    -- TeX
    gh("KeitaNakamura/tex-conceal.vim"),
    gh("lervag/vimtex"),

    -- Tree
    gh("nvim-tree/nvim-tree.lua"),

    -- Treesitter
    { src = gh("nvim-treesitter/nvim-treesitter"), version = "main" },
    { src = gh("nvim-treesitter/nvim-treesitter-textobjects"), version = "main" },
    gh("romgrk/nvim-treesitter-context"),

    -- UI
    gh("tyru/capture.vim"),
    gh("stevearc/oil.nvim"),
    gh("A7Lavinraj/fyler.nvim"),
    gh("jake-stewart/auto-cmdheight.nvim"),
    gh("b0o/incline.nvim"),

    -- Completion (nvim-cmp ecosystem)
    gh("hrsh7th/nvim-cmp"),
    gh("saghen/blink.compat"),
    gh("folke/lazydev.nvim"),
    gh("Kaiser-Yang/blink-cmp-dictionary"),
    { src = gh("saghen/blink.cmp"), version = vim.version.range("1") },
    gh("saadparwaiz1/cmp_luasnip"),
    gh("iguanacucumber/mag-nvim-lsp"),
    gh("iguanacucumber/mag-nvim-lua"),
    gh("iguanacucumber/mag-buffer"),
    gh("iguanacucumber/mag-cmdline"),
    gh("onsails/lspkind.nvim"),
    "https://codeberg.org/FelipeLema/cmp-async-path",
    gh("hrsh7th/cmp-nvim-lsp-signature-help"),
    gh("hrsh7th/cmp-omni"),
    gh("tamago324/cmp-zsh"),
    gh("hrsh7th/cmp-calc"),
    gh("ray-x/cmp-treesitter"),
    gh("chrisgrieser/cmp-nerdfont"),
    gh("KadoBOT/cmp-plugins"),
    gh("uga-rosa/cmp-dictionary"),
    gh("honza/vim-snippets"),
    gh("rafamadriz/friendly-snippets"),
    gh("degnbol/cmp-buffer"),
    gh("hrsh7th/cmp-nvim-lua"),
    gh("davidmh/cmp-nerdfonts"),

    -- Agents
    gh("coder/claudecode.nvim"),

    -- Shared dependencies
    gh("nvim-lua/plenary.nvim"),
})
