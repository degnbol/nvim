local map = require "utils/keymap"

return {
    {
        -- :Git and similar commands
        "tpope/vim-fugitive",
        cmd = { "Git", "G", "Gedit", "Gsplit", "Gread", "Gwrite", "Ggrep", "GMove", "GDelete", "GRemove" },
        init = function()
            -- using same naming that I have in the terminal.
            -- abbrev since a cmd will need to start with uppercase and cannot be modified after expansion.
            -- a "c" keymap is also possible but would expand at any mention and make typing anything with g a bit odd.
            -- There's no abbrev lua function yet.
            vim.cmd [[
                cnoreabbrev gts Git status -uno
                cnoreabbrev gtl Git pull
                cnoreabbrev gtp Git push
                cnoreabbrev gta Git add %
                cnoreabbrev gtc Git commit -m
            ]]
        end,
    },
    {
        "NeogitOrg/neogit",
        cmd = "Neogit",
        dependencies = {
            "nvim-lua/plenary.nvim",         -- required
            "nvim-telescope/telescope.nvim", -- optional
            "sindrets/diffview.nvim",        -- optional
            "ibhagwan/fzf-lua",              -- optional
        },
        init = function()
            map.n('<leader>gn', "<Cmd>Neogit<CR>", "Neogit")
        end,
        config = function()
            -- it seems NeogitDiffAdd and Delete are set to color fg + no bg version of DiffAdd and DiffDelete (DiffAdd by default is the same but DiffDelete has bg and no fg).
            -- So: no need to do anything to coordinate colors.
            require "neogit".setup {
                -- go directly to insert mode after pressing cc for commit, if the commit message is empty
                disable_insert_on_commit = "auto",
                disable_context_highlighting = true,
                signs = {
                    -- { CLOSED, OPENED }
                    section = { "", "" },
                    item = { "", "" },
                    hunk = { "", "" },
                },
                integrations = { diffview = true }, -- adds integration with diffview.nvim
                sections = {
                    untracked = {
                        folded = true,
                    }
                },
            }
        end,
    },
    -- highlight git conflicts, jump with [x and ]x,
    -- resolve by keeping none (cn), theirs (ct), our (co), both (cb), or both reverse (cB)
    {
        "rhysd/conflict-marker.vim",
        event = "VeryLazy",
    },
    -- git decoration to the left
    {
        "lewis6991/gitsigns.nvim",
        event = "VeryLazy",
        dependencies = { 'nvim-lua/plenary.nvim' },
        init = function()
            local gs = require "gitsigns"

            map.n("<leader>ga", gs.stage_hunk, "Add/stage hunk")
            map.n("<leader>gA", gs.stage_buffer, "Add/stage buffer")
            map.n("<leader>gu", gs.undo_stage_hunk, "Unstage hunk")
            map.n("<leader>gb", gs.blame_line, "Blame line")
            map.n("<leader>gp", gs.preview_hunk, "Preview hunk")
            map.n("<leader>gi", gs.preview_hunk_inline, "Preview hunk inline")
            map.n("<leader>gr", gs.reset_hunk, "Reset hunk")
            map.n("<leader>gR", gs.reset_buffer, "Reset buffer")
            -- not as good as diffview
            -- map.n("<leader>gd", gs.diffthis       , "Diff buffer")
            map.n("<leader>gt", gs.toggle_deleted, "Toggle deleted")
            map.n("[h", function() gs.nav_hunk('prev') end, "Previous hunk")
            map.n("]h", function() gs.nav_hunk('next') end, "Next hunk")
            map.x("<leader>gs", function() gs.stage_hunk { vim.fn.line('.'), vim.fn.line('v') } end, "Stage hunks")
            map.x("<leader>gr", function() gs.reset_hunk { vim.fn.line('.'), vim.fn.line('v') } end, "Reset hunks")
            -- Text object
            map.ox('ih', ':<C-U>Gitsigns select_hunk<CR>', "Hunk")
        end,
        opts = {
            -- by default they are '~' to indicate that try to be compromise between bar and underscore.
            -- We can show underline instead as a better compromise. See lua/highlights
            -- We set delete to nothing, istead of default underscore, and use underline there as well for consistency.
            signs         = {
                changedelete = { text = '┃' },
                delete       = { text = ' ' },
            },
            signs_staged  = {
                changedelete = { text = '┃' },
                delete       = { text = ' ' },
            },
            numhl         = true, -- highlight line number
            watch_gitdir  = { interval = 100 },
            sign_priority = 5,
        },
    },
    -- :DiffviewOpen (<leader>gd) and other commands for seeing git diff and git history for files.
    {
        'sindrets/diffview.nvim',
        dependencies = { 'nvim-lua/plenary.nvim', 'nvim-tree/nvim-web-devicons' },
        -- cmd = {"DiffviewOpen", "DiffviewClose", "DiffviewFileHistory", "DiffviewFocusFiles", "DiffviewToggleFiles", "DiffviewLog", "DiffviewRefresh"},
        -- only some commands are relevant before first use
        cmd = { "DiffviewOpen", "DiffviewFileHistory" },
        init = function()
            -- hide untracked files with -uno.
            -- hide gitsigns' file explorer with DiffviewToggleFiles (unhide with <leader>e like NvimTreeToggle)
            -- Open during merge or rebase should show conflicts nicer automatically.
            map.n("<leader>gd", "<Cmd>DiffviewOpen -uno<CR>:DiffviewToggleFiles<CR>", "Diffview open")
            -- pretty cool: works on ranges
            map({ "n", "x" }, "<leader>gh", "<Cmd>DiffviewFileHistory %<CR>", { desc = "History" })
        end,
        config = function()
            local diffview = require "diffview"
            local actions = require "diffview.actions"
            -- defaults copied and changed from
            -- https://github.com/sindrets/diffview.nvim
            diffview.setup {
                -- more subtle coloring knowing the context is git diff between old and new and not just file diff
                -- See ':h diffview-config-enhanced_diff_hl'
                enhanced_diff_hl = true,
                keymaps = {
                    -- Customise them here:
                    -- - small change for <leader>e and <leader>b
                    -- - only set the <leader>gq here
                    disable_defaults = true,
                    view = {
                        -- The `view` bindings are active in the diff buffers, only when the current
                        -- tabpage is a Diffview.
                        ["<tab>"]          = actions.select_next_entry,         -- Open the diff for the next file
                        ["<s-tab>"]        = actions.select_prev_entry,         -- Open the diff for the previous file
                        ["gf"]             = actions.goto_file,                 -- Open the file in a new split in the previous tabpage
                        ["<C-w><C-f>"]     = actions.goto_file_split,           -- Open the file in a new split
                        ["<C-w>gf"]        = actions.goto_file_tab,             -- Open the file in a new tabpage
                        ["<leader>e"]      = actions.focus_files,               -- Bring focus to the file panel
                        -- ["<leader>e"]  = actions.toggle_files,              -- Toggle the file panel.
                        ["g<C-x>"]         = actions.cycle_layout,              -- Cycle through available layouts.
                        ["[x"]             = actions.prev_conflict,             -- In the merge_tool: jump to the previous conflict
                        ["]x"]             = actions.next_conflict,             -- In the merge_tool: jump to the next conflict
                        ["<LocalLeader>o"] = actions.conflict_choose("ours"),   -- Choose the OURS version of a conflict
                        ["<LocalLeader>t"] = actions.conflict_choose("theirs"), -- Choose the THEIRS version of a conflict
                        ["<LocalLeader>b"] = actions.conflict_choose("base"),   -- Choose the BASE version of a conflict
                        ["<LocalLeader>a"] = actions.conflict_choose("all"),    -- Choose all the versions of a conflict
                        ["dx"]             = actions.conflict_choose("none"),   -- Delete the conflict region
                        ["<LocalLeader>d"] = actions.conflict_choose("none"),   -- Use dx instead
                        ["<Leader>gq"]     = actions.close,                     -- Use :tabclose instead
                        ["<LocalLeader>q"] = actions.close,                     -- Use :tabclose instead
                        ["<Leader>bc"]     = "<Cmd>tabclose<CR>",
                    },
                    diff1 = { --[[ Mappings in single window diff layouts ]] },
                    diff2 = { --[[ Mappings in 2-way diff layouts ]] },
                    diff3 = {
                        -- Mappings in 3-way diff layouts
                        { { "n", "x" }, "2do", actions.diffget("ours") },   -- Obtain the diff hunk from the OURS version of the file
                        { { "n", "x" }, "3do", actions.diffget("theirs") }, -- Obtain the diff hunk from the THEIRS version of the file
                    },
                    diff4 = {
                        -- Mappings in 4-way diff layouts
                        { { "n", "x" }, "1do", actions.diffget("base") },   -- Obtain the diff hunk from the BASE version of the file
                        { { "n", "x" }, "2do", actions.diffget("ours") },   -- Obtain the diff hunk from the OURS version of the file
                        { { "n", "x" }, "3do", actions.diffget("theirs") }, -- Obtain the diff hunk from the THEIRS version of the file
                    },
                    file_panel = {
                        ["j"]              = actions.next_entry,   -- Bring the cursor to the next file entry
                        ["<down>"]         = actions.next_entry,
                        ["k"]              = actions.prev_entry,   -- Bring the cursor to the previous file entry.
                        ["<up>"]           = actions.prev_entry,
                        ["<cr>"]           = actions.select_entry, -- Open the diff for the selected entry.
                        ["o"]              = actions.select_entry,
                        ["<2-LeftMouse>"]  = actions.select_entry,
                        ["-"]              = actions.toggle_stage_entry, -- Stage / unstage the selected entry.
                        ["S"]              = actions.stage_all,          -- Stage all entries.
                        ["U"]              = actions.unstage_all,        -- Unstage all entries.
                        ["X"]              = actions.restore_entry,      -- Restore entry to the state on the left side.
                        ["L"]              = actions.open_commit_log,    -- Open the commit log panel.
                        ["<c-b>"]          = actions.scroll_view(-0.25), -- Scroll the view up
                        ["<c-f>"]          = actions.scroll_view(0.25),  -- Scroll the view down
                        ["<tab>"]          = actions.select_next_entry,
                        ["<s-tab>"]        = actions.select_prev_entry,
                        ["gf"]             = actions.goto_file,
                        ["<C-w><C-f>"]     = actions.goto_file_split,
                        ["<C-w>gf"]        = actions.goto_file_tab,
                        ["i"]              = actions.listing_style,       -- Toggle between 'list' and 'tree' views
                        ["f"]              = actions.toggle_flatten_dirs, -- Flatten empty subdirectories in tree listing style.
                        ["R"]              = actions.refresh_files,       -- Update stats and entries in the file list.
                        -- ["<leader>e"]     = actions.focus_files,
                        ["<leader>e"]      = actions.toggle_files,
                        ["g<C-x>"]         = actions.cycle_layout,
                        ["[x"]             = actions.prev_conflict,
                        ["]x"]             = actions.next_conflict,
                        ["<leader>gq"]     = "<Cmd>tabclose<CR>",
                        ["<localleader>q"] = "<Cmd>tabclose<CR>",
                        ["<leader>bc"]     = "<Cmd>tabclose<CR>",
                    },
                    file_history_panel = {
                        ["g!"]             = actions.options,          -- Open the option panel
                        ["<C-A-d>"]        = actions.open_in_diffview, -- Open the entry under the cursor in a diffview
                        ["y"]              = actions.copy_hash,        -- Copy the commit hash of the entry under the cursor
                        ["L"]              = actions.open_commit_log,
                        ["zR"]             = actions.open_all_folds,
                        ["zM"]             = actions.close_all_folds,
                        ["j"]              = actions.next_entry,
                        ["<down>"]         = actions.next_entry,
                        ["k"]              = actions.prev_entry,
                        ["<up>"]           = actions.prev_entry,
                        ["<cr>"]           = actions.select_entry,
                        ["o"]              = actions.select_entry,
                        ["<2-LeftMouse>"]  = actions.select_entry,
                        ["<c-b>"]          = actions.scroll_view(-0.25),
                        ["<c-f>"]          = actions.scroll_view(0.25),
                        ["<tab>"]          = actions.select_next_entry,
                        ["<s-tab>"]        = actions.select_prev_entry,
                        ["gf"]             = actions.goto_file,
                        ["<C-w><C-f>"]     = actions.goto_file_split,
                        ["<C-w>gf"]        = actions.goto_file_tab,
                        ["<leader>e"]      = actions.focus_files,
                        -- ["<leader>e"]     = actions.toggle_files,
                        ["g<C-x>"]         = actions.cycle_layout,
                        ["<leader>gq"]     = "<Cmd>tabclose<CR>",
                        ["<localleader>q"] = "<Cmd>tabclose<CR>",
                        ["<leader>bc"]     = "<Cmd>tabclose<CR>",
                    },
                    option_panel = {
                        ["<tab>"] = actions.select_entry,
                        ["q"]     = actions.close,
                    },
                },
            }
        end,
    },
}
