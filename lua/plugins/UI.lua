#!/usr/bin/env lua
return {
    {"akinsho/nvim-bufferline.lua", dependencies={"DaikyXendo/nvim-material-icon"}, config=function() require'top-bufferline' end}, -- add a line at the top with all the files open in the buffer
    -- {"glepnir/galaxyline.nvim", config=function() require'statusline' end},
    {"nvim-telescope/telescope-fzf-native.nvim", build='make'}, -- recommended compiled fuzzy finder for telescope. Cannot be opt=true when needed by tzachar/cmp-fuzzy-path
    -- {"nvim-telescope/telescope.nvim", dependencies={"nvim-lua/plenary.nvim", "nvim-telescope/telescope-fzf-native.nvim"}, config=function() require'telescope-conf' end}, -- Fuzzy finder
    {"glepnir/dashboard-nvim", config=function() require'dashboard-conf' end}, -- open to a dashboard for vi without a file selection, requires telescope or an alternative installed.
    {"DaikyXendo/nvim-tree.lua", dependencies={'DaikyXendo/nvim-material-icon'}, config=function() require'tree' end}, -- tree file explorer to the left. A more featured alternative: https://github.com/ms-jpq/chadtree
    "ojroques/nvim-bufdel", -- :BufDel that deletes a buffer better than built-in :bdelete and :bwipeout, by preserving layout and closing terminal buffers better.
    {'sudormrfbin/cheatsheet.nvim', dependencies={'nvim-telescope/telescope.nvim'}}, -- <leader>? to give cheatsheet popup. 
    {"lalitmee/browse.nvim", dependencies={"nvim-telescope/telescope.nvim"}}, -- search stackoverflow quicker
    -- {"kevinhwang91/nvim-hlslens", config=function() require'hlslens-conf' end}, -- show search match numbers
    -- {"petertriho/nvim-scrollbar", dependencies={"kevinhwang91/nvim-hlslens"}, config=function() require'scrollbar-conf' end}, -- requires hlslens to show search results in scrollbar
    "tyru/capture.vim", -- :Capture hi to call :hi where you can search etc.
    {dir="~/nvim/kittyREPL.nvim", config=function() require'kittyREPL-conf' end},
    {'kevinhwang91/nvim-bqf', config=function() require'bqf-conf' end, dependencies={'junegunn/fzf'}}, -- better quickfix window. zf to open fzf inside quickfix.
    {'kevinhwang91/nvim-ufo', dependencies={'kevinhwang91/promise-async'}, config=function() require'ufo-conf' end}, -- ultra fold
}
