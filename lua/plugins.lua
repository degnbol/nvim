-- packer as package manager as opposed to packages.lua
return require("packer").startup(
    function()
        use "wbthomason/packer.nvim"

        -- color
        use "norcalli/nvim-colorizer.lua" -- when a hex or other color is defined, highlight the text with its color 
        use "siduck76/nvim-base16.lua"
        use "maxwells-daemons/base16-gigavolt-scheme"
        
        use "Pocco81/TrueZen.nvim" -- reduce visuals with TZ... commands to e.g. remove left and bottom element on the screen.

        -- language
        use "nvim-treesitter/nvim-treesitter" -- language coloring and ensuring of installation
        use "neovim/nvim-lspconfig" -- lsp
        use "kabouzeid/nvim-lspinstall" -- adds :LspInstall <language> for conveniently installing language support
        use "hrsh7th/nvim-compe"  -- adds autocompletion
        use "onsails/lspkind-nvim" -- VS code like pictograms for completion
        use "terrortylor/nvim-comment" -- Toggle commenting out code
        use "windwp/nvim-autopairs" -- auto add second parenthesis etc.
        use "lukas-reineke/indent-blankline.nvim" -- show | on indented lines
        use 'tpope/vim-fugitive' -- git
        use "lewis6991/gitsigns.nvim"  -- git decoration
        
        use "folke/which-key.nvim"  -- pop-up to help with keybindings that have been started

        -- UI
        use {'akinsho/nvim-bufferline.lua', requires = 'kyazdani42/nvim-web-devicons'} -- add a line at the top with all the files open in the buffer
        use "glepnir/galaxyline.nvim"
        -- Fuzzy finder
        use {
            'nvim-telescope/telescope.nvim',
            requires = {{'nvim-lua/popup.nvim'}, {'nvim-lua/plenary.nvim'}}
        }
        use "glepnir/dashboard-nvim" -- open to a dashboard for vi without a file selection, requires telescope or an alternative installed.

        -- file managing, picker etc
        use "ryanoasis/vim-devicons" -- adds icons to files
        use "kyazdani42/nvim-tree.lua" -- tree view window for file exploring with bug in .config files

        use "tweekmonster/startuptime.vim"  -- use :StartupTime to measure what things are affecting startup time
    end,
    {
        display = {
            border = {"┌", "─", "┐", "│", "┘", "─", "└", "│"}
        }
    }
)
