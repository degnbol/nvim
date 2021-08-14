-- Paq as package manager as opposed to Packer in plugins.lua
require "paq-nvim" {
	"savq/paq-nvim";

    -- typing behaviour
	"folke/which-key.nvim"; -- pop-up to help with keybindings that have been started
    "tpope/vim-repeat"; -- change . to repeat last native command to last "full" command, which feels more natural.
    "tpope/vim-surround"; -- press cs'" to change surrounding ' with ", ds' to delete surrounding ', ysiw) to surround word with ) and yss[ to surround line with [ ... ] (incl. spaces)
    "svermeulen/vim-subversive"; -- add substitution functions to e.g. replace a word with clipboard content by writing siw
    "svermeulen/vim-cutlass"; -- c(hange), d(elete) no longer copies, remapped in keymapping file so x will cut. Since we have added backspace and delete button support in normal mode there is no need for default x behavior
    "svermeulen/vim-yoink"; -- yank history that you can cycle with c-n and c-p
    "mg979/vim-visual-multi"; -- multi cursor TODO https://github.com/mg979/vim-visual-multi/wiki/Quick-start

	-- color
	"norcalli/nvim-colorizer.lua"; -- when a hex or other color is defined, highlight the text with its color
	"siduck76/nvim-base16.lua";
	"maxwells-daemons/base16-gigavolt-scheme";

	"Pocco81/TrueZen.nvim"; -- reduce visuals with :TZ... commands to e.g. remove left and bottom element on the screen.

	"ryanoasis/vim-devicons"; -- adds icons to files

    "sakshamgupta05/vim-todo-highlight"; -- highlight todos
    "p00f/nvim-ts-rainbow"; -- tree sitter based rainbow color parenthesis to easily see the matching
    "folke/twilight.nvim"; -- dim code that isn't currently being edited with :Twilight.

	-- language
	"nvim-treesitter/nvim-treesitter"; -- language coloring and ensuring of installation
    "nvim-treesitter/nvim-treesitter-refactor"; -- refactor
    "nvim-treesitter/nvim-treesitter-textobjects"; -- selecting, moving functions etc.
    -- "romgrk/nvim-treesitter-context"; -- show the "context" at the top line, i.e. function name when in a function
	"neovim/nvim-lspconfig"; -- lsp
	"kabouzeid/nvim-lspinstall"; -- adds :LspInstall <language> for conveniently installing language support
	"hrsh7th/nvim-compe";  -- adds autocompletion. It is an alternative to nvim-lua/completion-nvim which online discussions say is slower.
    "ray-x/lsp_signature.nvim"; -- hover signatures for function arguments. TODO make it work
	"onsails/lspkind-nvim"; -- VS code like pictograms for completion
	"terrortylor/nvim-comment"; -- Toggle commenting out code
	-- "windwp/nvim-autopairs"; -- auto add second parenthesis etc.
	-- "lukas-reineke/indent-blankline.nvim"; -- show "|" on indented lines
	"tpope/vim-fugitive"; -- git
	"lewis6991/gitsigns.nvim"; -- git decoration to the left
    "JuliaEditorSupport/julia-vim"; -- julia support, colors and unicode substitution.
    -- "urbainvaes/vim-ripple"; -- REPL with some indent and tab problems
    "hkupty/iron.nvim"; -- REPL that doesn't support bpython
    "pappasam/nvim-repl"; -- REPL that has to be started
    "jeetsukumaran/vim-pythonsense"; -- python aware changes to [], [[, ]], ][, ]m, ]M, [m, [M for moving cursor to starts and ends of python functions. This should be covered by tree sitter in the future when they add support for visual mode
    "samirettali/shebang.nvim"; -- insert shebang on new file edit
    -- 'metakirby5/codi.vim'; -- scratchpad coding, see output of all lines to the right https://github.com/metakirby5/codi.vim

    "RishabhRD/popfix"; -- required by nvim-cheat
    "RishabhRD/nvim-cheat.sh"; -- try it out with :Cheat <query> where the query should be search terms like you would search in StackOverflow for answers

    -- "elzr/vim-json"; -- json

	-- UI
	"kyazdani42/nvim-web-devicons"; -- required by nvim-bufferline
	"akinsho/nvim-bufferline.lua"; -- add a line at the top with all the files open in the buffer
	"glepnir/galaxyline.nvim";
	"nvim-lua/popup.nvim";
	"nvim-lua/plenary.nvim";
	"nvim-telescope/telescope.nvim"; -- Fuzzy finder
	"glepnir/dashboard-nvim"; -- open to a dashboard for vi without a file selection, requires telescope or an alternative installed.
    "kyazdani42/nvim-tree.lua"; -- tree file explorer to the left
    "ojroques/nvim-bufdel"; -- :BufDel that deletes a buffer better than built-in :bdelete and :bwipeout, by preserving layout and closing terminal buffers better.

    -- "khzaw/vim-conceal"; -- add some language specific conceal features
    -- "ehamberg/vim-cute-python"; 

	-- "tweekmonster/startuptime.vim";  -- use :StartupTime to measure what things are affecting startup time
}
