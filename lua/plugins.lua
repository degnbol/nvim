-- packer as package manager as opposed to packages.lua
-- Run PackerSync after making changes to this file to recompile the file in plugin
return require("packer").startup(function()
    use "wbthomason/packer.nvim"
    
    -- core behaviour
    use "tpope/vim-repeat" -- change . to repeat last native command to last "full" command, which feels more natural.
    -- problems: it doesn't have simple ds working at all.
    -- use {"kylechui/nvim-surround", config=function() require"surround-conf" end} -- press cs'" to change surrounding ' with ", ds' to delete surrounding ', ysiw) to surround word with ) and yss[ to surround line with [ ... ] (incl. spaces)
    use { 'echasnovski/mini.nvim', branch = 'stable', config=function() require'mini-conf' end } -- mini.surround
    -- use "tpope/vim-sensible" -- Y should yank to end of line which is consistent with other uppercase use, rather than yank whole line like yy which is for ancient vi compatibility.
    use "svermeulen/vim-subversive" -- add substitution functions to e.g. replace a word with clipboard content by writing siw
    use {"gbprod/cutlass.nvim", config=function() require'cutlass-conf' end} -- c(hange), d(elete) no longer copies, remapped in keymapping file so x will cut. Since we have added backspace and delete button support in normal mode there is no need for default x behavior
    use "svermeulen/vim-yoink" -- yank history that you can cycle with c-n and c-p
    -- use "mg979/vim-visual-multi" -- multi cursor TODO https://github.com/mg979/vim-visual-multi/wiki/Quick-start
    -- use {"moll/vim-bbye", config=function() require'bbye' end}
    use "farmergreg/vim-lastplace" -- open file in last edited location
    use {"inkarkat/vim-UnconditionalPaste", requires='inkarkat/vim-ingo-library', config=function() require'unconditionalPaste' end} -- lots of ways to paste using g{c,C,l,b}{,i}{p,P} and may others 
    use {"AndrewRadev/whitespaste.vim", config=function() require'whitespaste' end} -- paste without empty newlines
    use "google/vim-searchindex" -- let's search result box show number of matches when there's >99 matches
    use "haya14busa/vim-asterisk" -- improvements to z* and visual *. See git for uses https://github.com/haya14busa/vim-asterisk
    -- use "kana/vim-textobj-user" -- easily define custom textobjects such as i( and a( to select in/an \left( \right) block in latex
    -- TODO add from https://github.com/kana/vim-textobj-user and https://github.com/kana/vim-textobj-user/wiki
    use {"glts/vim-textobj-comment", requires="kana/vim-textobj-user"} -- not working?
    use {"AckslD/nvim-trevJ.lua", config=function() require'trevj-conf' end}
    use {"monaqa/dial.nvim", config=function() require'dial-conf' end} -- increment and decrement numbers, dates, color hex, even bool
    
    -- color
    use {"norcalli/nvim-colorizer.lua", config=function() require'colorizer'.setup() end} -- when a hex or other color is defined, highlight the text with its color
    use {"norcalli/nvim-base16.lua", requires="norcalli/nvim.lua"}
    use {"unblevable/quick-scope", config=function() require'quick-scope' end} -- highlight letters for jumping with f/F/t/T
    
    use "kyazdani42/nvim-web-devicons" -- adds icons to files
    
    use "sakshamgupta05/vim-todo-highlight" -- highlight todos
    -- use {"folke/twilight.nvim", config=function() require'twilight'.setup{dimming={alpha=0.5}, context=30} end} -- dim code that isn't currently being edited with :Twilight.
    -- use {"p00f/nvim-ts-rainbow", requires='nvim-treesitter/nvim-treesitter', config=function() require'treesitter-rainbow' end} -- tree sitter based rainbow color parenthesis to easily see the matching
    
    -- UI
    use {"akinsho/nvim-bufferline.lua", tag="*", requires="kyazdani42/nvim-web-devicons", config=function() require'top-bufferline' end} -- add a line at the top with all the files open in the buffer
    -- use {"glepnir/galaxyline.nvim", config=function() require'statusline' end}
    use {"nvim-telescope/telescope-fzf-native.nvim", run='make'} -- recommended compiled fuzzy finder for telescope. Cannot be opt=true when needed by tzachar/cmp-fuzzy-path
    use {"nvim-telescope/telescope.nvim", requires={"nvim-lua/plenary.nvim", "nvim-telescope/telescope-fzf-native.nvim"}, config=function() require'telescope-conf' end} -- Fuzzy finder
    use {"glepnir/dashboard-nvim", config=function() require'dashboard-conf' end} -- open to a dashboard for vi without a file selection, requires telescope or an alternative installed.
    use {"kyazdani42/nvim-tree.lua", requires='kyazdani42/nvim-web-devicons', config=function() require'tree' end} -- tree file explorer to the left. A more featured alternative: https://github.com/ms-jpq/chadtree
    use "ojroques/nvim-bufdel" -- :BufDel that deletes a buffer better than built-in :bdelete and :bwipeout, by preserving layout and closing terminal buffers better.
    use {"folke/which-key.nvim", config=function() require'whichkey' end} -- pop-up to help with keybindings that have been started
    use {'sudormrfbin/cheatsheet.nvim', requires='nvim-telescope/telescope.nvim'} -- <leader>? to give cheatsheet popup. 
    use {"lalitmee/browse.nvim", requires="nvim-telescope/telescope.nvim"} -- search stackoverflow quicker
    -- use {"kevinhwang91/nvim-hlslens", config=function() require'hlslens-conf' end} -- show search match numbers
    -- use {"petertriho/nvim-scrollbar", requires="kevinhwang91/nvim-hlslens", config=function() require'scrollbar-conf' end} -- requires hlslens to show search results in scrollbar
    use "tyru/capture.vim" -- :Capture hi to call :hi where you can search etc.
    use {"~/nvim/kittyREPL.nvim", config=function() require'kittyREPL-conf' end}
    use {'kevinhwang91/nvim-bqf', config=function() require'bqf-conf' end, requires='junegunn/fzf'} -- better quickfix window. zf to open fzf inside quickfix.
    
    -- treesitter
    use {'nvim-treesitter/nvim-treesitter', run=':TSUpdate', config=function() require'treesitter' end} -- language coloring and ensuring of installation
    use {"nvim-treesitter/nvim-treesitter-refactor", requires='nvim-treesitter/nvim-treesitter', config=function() require'treesitter-refactor' end} -- refactor
    use {"nvim-treesitter/nvim-treesitter-textobjects", requires='nvim-treesitter/nvim-treesitter', config=function() require'treesitter-textobjects' end} -- selecting, moving functions etc.
    use {"RRethy/nvim-treesitter-textsubjects", requires='nvim-treesitter/nvim-treesitter', config=function() require'treesitter-textsubjects' end} -- in vis mode use . , ; i; to select based on treesitter 
    -- use "romgrk/nvim-treesitter-context" -- show the "context" at the top line, i.e. function name when in a function
    use {"andymass/vim-matchup", requires='nvim-treesitter/nvim-treesitter', config=function() require'matchup' end} -- % jumps between matching coding blocks, not just single chars.
    use {"p00f/nvim-ts-rainbow", requires='nvim-treesitter/nvim-treesitter', config=function() require'treesitter-rainbow' end} -- tree sitter based rainbow color parenthesis to easily see the matching
    
    -- LSP. For a given file, either complete with cmp (builtin recommended, but e.g. jedi language servers is slow), coc (old, not using builtin LSP), or coq (hacks builtin LSP)
    use {"neovim/nvim-lspconfig", config=function() require'lsp' end}
    -- add :LspInstall <language> and :Mason for conveniently installing LSP language specific servers
    use {"williamboman/mason-lspconfig.nvim", requires={"neovim/nvim-lspconfig", "williamboman/mason.nvim"}, config=function() require "mason-conf" end}
    -- completion menu using builtin LSP
    use {"hrsh7th/nvim-cmp", requires = {
        'hrsh7th/cmp-nvim-lsp', 'hrsh7th/cmp-buffer',
        'hrsh7th/cmp-path',
        'hrsh7th/cmp-nvim-lsp-signature-help',
        'hrsh7th/cmp-nvim-lua', -- neovim Lua API
        'tamago324/cmp-zsh', -- neovim zsh completion
        'onsails/lspkind.nvim', -- pretty pictograms
        'hrsh7th/cmp-calc', -- quick math in completion
    }, config=function() require'cmp-conf' end}
    use {'saadparwaiz1/cmp_luasnip', requires={'L3MON4D3/LuaSnip', "hrsh7th/nvim-cmp"}, config=function() require"luasnip-conf" end } -- alts: hrsh7th/vim-vsnip, SirVer/ultisnips, ...
    -- custom dicts and spell check that doesn't require spell and spelllang (f3fora/cmp-spell)
    use {'uga-rosa/cmp-dictionary', config=function() require"cmp_dict" end}
    use {"j-hui/fidget.nvim", config=function() require"fidget-conf" end} -- corner print what LSP is running
    use "rafamadriz/friendly-snippets"
    
    -- language
    use {"terrortylor/nvim-comment", config=function() require'nvim_comment'.setup() end} -- add keybindings to toggle comments with motions etc.
    -- use "windwp/nvim-autopairs" -- auto add second parenthesis etc.
    -- use "lukas-reineke/indent-blankline.nvim" -- show "|" on indented lines
    use "tpope/vim-fugitive" -- git
    -- use {"lewis6991/gitsigns.nvim", requires='nvim-lua/plenary.nvim', config=function() require'gitsigns-nvim' end} -- git decoration to the left
    use {"JuliaEditorSupport/julia-vim", config=function() require'julia' end} -- julia support, colors and unicode substitution. CANNOT use ft=julia
    -- use "urbainvaes/vim-ripple" -- REPL with some indent and tab problems
    -- use {"hkupty/iron.nvim", config=function () require'repl/iron-nvim' end} -- REPL that doesn't support bpython or radian
    -- use {"pappasam/nvim-repl", config=function() require'repl/pappasam_repl' end} -- REPL that has to be started and can only send whole lines
    -- use {"kassio/neoterm", config=function() require'repl/kassio_neoterm' end}
    -- use {"HiPhish/repl.nvim", config=function () require'repl/HiPhish_repl' end}
    -- use {"jpalardy/vim-slime", config=function() require'repl/slime' end, requires='rxi/json.lua'} -- send code to REPL that can even be in another window.
    -- use "metakirby5/codi.vim" -- scratchpad coding, see output of all lines to the right https://github.com/metakirby5/codi.vim
    use {"jeetsukumaran/vim-pythonsense", ft='python'} -- python aware changes to [], [[, ]], ][, ]m, ]M, [m, [M for moving cursor to starts and ends of python functions. This should be covered by tree sitter in the future when they add support for visual mode
    use {"samirettali/shebang.nvim", config=function() require'shebang-nvim' end} -- insert shebang on new file edit
    -- try it out with :Cheat <query> where the query should be search terms like you would search in StackOverflow for answers
    use {"RishabhRD/nvim-cheat.sh", config=function() require'cheat' end, requires="RishabhRD/popfix"}
    -- use {"mrjones2014/dash.nvim", run='make install', requires='nvim-telescope/telescope.nvim'} -- :DashWord with <leader>K. conf in telescope-conf.lua
    use {"lervag/vimtex", config=function() require'vimtex' end} -- :VimtexCompile. Adds so much more good stuff, e.g. dse, cse to delete or change surrounding env
    -- use "tpope/vim-sleuth" -- sleuth that let's you autodetect if file is using 2 or 4 spaces. Mistakenly set noexpandtab
    use "tpope/vim-abolish" -- crs: snake_case, crm: MixedCase, crc: camelCase, cru: UPPER_CASE, cr-: dash-case, cr.: dot.case, cr<SPACE>: space case, crt: Title Case
    -- use {"preservim/vim-markdown", requires="godlygeek/tabular"} -- conceal markdown expressions like _emphasis_ and folding. Overkill, see {after/,}syntax/markdown.vim
    use {"iamcco/markdown-preview.nvim", run=':call mkdp#util#install()', ft='markdown'} -- :MarkdownPreview live in browser
    use {"habamax/vim-asciidoctor", config=function() require'asciidoc' end, ft='asciidoctor'}
    use {"quarto-dev/quarto-vim", requires="vim-pandoc/vim-pandoc-syntax", ft="quarto"} -- https://quarto.org/
    use {"habamax/vim-rst"}
    use "jbyuki/nabla.nvim" -- show pretty math in term
    use {"goerz/jupytext.vim", config=function() require'jupytext-conf' end} -- edit jupyter notebook. requires `pip install jupytext`
    -- use "elzr/vim-json" -- json
end)

