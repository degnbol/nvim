-- packer as package manager as opposed to packages.lua
-- Run PackerSync after making changes to this file to recompile the file in plugin
return require("packer").startup(
    function()
        use "wbthomason/packer.nvim"

        -- core behaviour
        use "tpope/vim-repeat" -- change . to repeat last native command to last "full" command, which feels more natural.
        use "tpope/vim-surround" -- press cs'" to change surrounding ' with ", ds' to delete surrounding ', ysiw) to surround word with ) and yss[ to surround line with [ ... ] (incl. spaces)
        -- use "tpope/vim-sensible" -- Y should yank to end of line which is consistent with other uppercase use, rather than yank whole line like yy which is for ancient vi compatibility.
        use "svermeulen/vim-subversive" -- add substitution functions to e.g. replace a word with clipboard content by writing siw
        use {"gbprod/cutlass.nvim", config=function() require'cutlass-conf' end} -- c(hange), d(elete) no longer copies, remapped in keymapping file so x will cut. Since we have added backspace and delete button support in normal mode there is no need for default x behavior
        use "svermeulen/vim-yoink" -- yank history that you can cycle with c-n and c-p
        -- use "mg979/vim-visual-multi" -- multi cursor TODO https://github.com/mg979/vim-visual-multi/wiki/Quick-start
        -- use {"moll/vim-bbye", config=function() require'bbye' end}
        -- use {"famiu/bufdelete.nvim", config=function() require'famiu_bufdelete' end}
        use "farmergreg/vim-lastplace" -- open file in last edited location
        use {"inkarkat/vim-UnconditionalPaste", requires='inkarkat/vim-ingo-library'} -- lots of ways to paste using g{c,C,l,b}{,i}{p,P} and may others 
        use "google/vim-searchindex" -- let's search result box show number of matches when there's >99 matches
        use "haya14busa/vim-asterisk" -- improvements to z* and visual *. See git for uses https://github.com/haya14busa/vim-asterisk
        use "kana/vim-textobj-user" -- easily define custom textobjects such as i( and a( to select in/an \left( \right) block in latex
        -- TODO: add from https://github.com/kana/vim-textobj-user and https://github.com/kana/vim-textobj-user/wiki

        -- color
        use "norcalli/nvim-colorizer.lua" -- when a hex or other color is defined, highlight the text with its color
        use "siduck76/nvim-base16.lua"
        use "maxwells-daemons/base16-gigavolt-scheme"
        use {"unblevable/quick-scope", config=function() require'quick-scope' end} -- highlight letters for jumping with f/F/t/T

        use "ryanoasis/vim-devicons" -- adds icons to files

        use "sakshamgupta05/vim-todo-highlight" -- highlight todos
        -- use {"Pocco81/TrueZen.nvim", config=function() require'truezen-nvim' end} -- reduce visuals with :TZ... commands to e.g. remove left and bottom element on the screen.
        -- use {"folke/twilight.nvim", config=function() require'twilight'.setup{dimming={alpha=0.5}, context=30} end} -- dim code that isn't currently being edited with :Twilight.
        -- use {"p00f/nvim-ts-rainbow", requires='nvim-treesitter/nvim-treesitter', config=function() require'treesitter-rainbow' end} -- tree sitter based rainbow color parenthesis to easily see the matching
        use {"p00f/nvim-ts-rainbow", requires='nvim-treesitter/nvim-treesitter', config=function() require'treesitter-rainbow' end} -- tree sitter based rainbow color parenthesis to easily see the matching

        -- language
        use {'nvim-treesitter/nvim-treesitter', run=':TSUpdate', config=function() require'treesitter' end} -- language coloring and ensuring of installation
        use {"nvim-treesitter/nvim-treesitter-refactor", requires='nvim-treesitter/nvim-treesitter', config=function() require'treesitter-refactor' end} -- refactor
        use {"nvim-treesitter/nvim-treesitter-textobjects", requires='nvim-treesitter/nvim-treesitter', config=function() require'treesitter-textobjects' end} -- selecting, moving functions etc.
        use {"RRethy/nvim-treesitter-textsubjects", requires='nvim-treesitter/nvim-treesitter', config=function() require'treesitter-textsubjects' end} -- in vis mode use . , ; i; to select based on treesitter 
        -- use "romgrk/nvim-treesitter-context" -- show the "context" at the top line, i.e. function name when in a function
        use {"andymass/vim-matchup", config=function() require'matchup' end} -- % jumps between matching coding blocks, not just single chars.
        use {"neovim/nvim-lspconfig"} -- lsp
        -- use {"neovim/nvim-lspconfig", -- lsp
        --     requires = {
        --         {'ms-jpq/coq_nvim', branch='coq', config=function() require'coq-nvim' end}, -- completion
        --         {'ms-jpq/coq.artifacts', branch='artifacts'}
        --     }
        -- }
        -- use {"kabouzeid/nvim-lspinstall", requires="neovim/nvim-lspconfig", config=function() require "lspinstall".setup() end} -- adds :LspInstall <language> for conveniently installing language support
        -- use {"hrsh7th/nvim-cmp", requires={"hrsh7th/cmp-path", "hrsh7th/cmp-nvim-lsp"}, config=function() require'nvim-cmp' end}  -- autocompletion
        use {"neoclide/coc.nvim", branch="release"} -- https://github.com/neoclide/coc.nvim/wiki/Language-servers e.g. :CocInstall coc-texlab
        -- -- use "ray-x/lsp_signature.nvim" -- hover signatures for function arguments. 
        -- use {"onsails/lspkind-nvim", config=function() require'lspkind'.init() end} -- VS code like pictograms for completion
        use {"terrortylor/nvim-comment", config=function() require'nvim_comment'.setup() end} -- Toggle commenting out code
        -- use "windwp/nvim-autopairs" -- auto add second parenthesis etc.
        -- use "lukas-reineke/indent-blankline.nvim" -- show "|" on indented lines
        use "tpope/vim-fugitive" -- git
        use {"lewis6991/gitsigns.nvim", requires='nvim-lua/plenary.nvim', config=function() require'gitsigns-nvim' end} -- git decoration to the left
        use {"JuliaEditorSupport/julia-vim", config=function() require'julia' end} -- julia support, colors and unicode substitution. CANNOT use ft=julia
        -- use "urbainvaes/vim-ripple" -- REPL with some indent and tab problems
        -- use {"hkupty/iron.nvim", config=function () require'iron-nvim' end} -- REPL that doesn't support bpython or radian
        -- use {"pappasam/nvim-repl", config=function() require'pappasam_repl' end} -- REPL that has to be started and can only send whole lines
        -- use {"kassio/neoterm", config=function() require'kassio_neoterm' end}
        -- use {"HiPhish/repl.nvim", config=function () require'HiPhish_repl' end}
        -- use {"jpalardy/vim-slime", config=function() require'slime' end, requires='rxi/json.lua'} -- send code to REPL that can even be in another window.
        -- use "metakirby5/codi.vim" -- scratchpad coding, see output of all lines to the right https://github.com/metakirby5/codi.vim
        use "rxi/json.lua" -- for kitty.lua
        use {"jeetsukumaran/vim-pythonsense", ft='python'} -- python aware changes to [], [[, ]], ][, ]m, ]M, [m, [M for moving cursor to starts and ends of python functions. This should be covered by tree sitter in the future when they add support for visual mode
        use {"samirettali/shebang.nvim", config=function() require'shebang-nvim' end} -- insert shebang on new file edit
        -- try it out with :Cheat <query> where the query should be search terms like you would search in StackOverflow for answers
        use {"RishabhRD/nvim-cheat.sh", config=function() require'cheat' end, requires="RishabhRD/popfix"}
        use {"lervag/vimtex", config=function() require'vimtex' end} -- :VimtexCompile
        -- use "tpope/vim-sleuth" -- sleuth that let's you autodetect if file is using 2 or 4 spaces. Mistakenly set noexpandtab
        -- use {"preservim/vim-markdown", ft='markdown'} -- conceal markdown expressions like _emphasis_ and folding. Overkill, see after/syntax/markdown.vim
        use {"iamcco/markdown-preview.nvim", run=':call mkdp#util#install()', ft='markdown'} -- :MarkdownPreview live in browser
        use {"habamax/vim-asciidoctor", config=function() require'asciidoc' end, ft='asciidoctor'}
        use {"habamax/vim-rst"}
        use "jbyuki/nabla.nvim" -- show pretty math in term

        -- use "elzr/vim-json" -- json

        -- UI
        use {"akinsho/nvim-bufferline.lua", tag = "*", requires="kyazdani42/nvim-web-devicons", config=function() require'top-bufferline' end} -- add a line at the top with all the files open in the buffer
        -- use {"glepnir/galaxyline.nvim", config=function() require'statusline' end}
        use {"nvim-telescope/telescope-fzf-native.nvim", run='make', opt=true} -- recommended compiled fuzzy finder for telescope
        use {"nvim-telescope/telescope.nvim", requires={"nvim-lua/plenary.nvim", "nvim-telescope/telescope-fzf-native.nvim"}, config=function() require'telescope-conf' end} -- Fuzzy finder
        use {"glepnir/dashboard-nvim", config=function() require'dashboard' end} -- open to a dashboard for vi without a file selection, requires telescope or an alternative installed.
        use {"kyazdani42/nvim-tree.lua", requires='kyazdani42/nvim-web-devicons', config=function() require'tree' end} -- tree file explorer to the left. A more featured alternative: https://github.com/ms-jpq/chadtree
        use {"ojroques/nvim-bufdel", config=function() require'ojroques_bufdel' end} -- :BufDel that deletes a buffer better than built-in :bdelete and :bwipeout, by preserving layout and closing terminal buffers better.
        use {"folke/which-key.nvim", config=function() require'whichkey' end} -- pop-up to help with keybindings that have been started
        use {'sudormrfbin/cheatsheet.nvim', requires='nvim-telescope/telescope.nvim'} -- <leader>? to give cheatsheet popup. 
        use {"lalitmee/browse.nvim", requires="nvim-telescope/telescope.nvim"} -- search stackoverflow quicker
        use {"kevinhwang91/nvim-hlslens", config=function() require'hlslens-conf' end} -- show search match numbers
        use {"petertriho/nvim-scrollbar", requires="kevinhwang91/nvim-hlslens", config=function() require'scrollbar-conf' end} -- requires hlslens to show search results in scrollbar
        use "tyru/capture.vim" -- :Capture hi to call :hi where you can search etc.
    end
)
