local g = vim.g

g.nvim_tree_indent_markers = 1
g.nvim_tree_git_hl = 1
-- some issue with highlighting shell file logo
-- g.nvim_tree_highlight_opened_files = 1
g.nvim_tree_allow_resize = 0

g.nvim_tree_show_icons = {git=1, folders=1, files=1}

g.nvim_tree_icons = {
    default = " ",
    symlink = " ",
    git = {
        unstaged = "✗",
        staged = "✓",
        unmerged = "",
        renamed = "➜",
        untracked = "★",
        deleted = "",
        ignored = "◌"
    },
    folder = {
        default = "",
        open = "",
        symlink = "",
        empty = "",
        empty_open = "",
        symlink_open = ""
    }
}

-- Mappings

vim.api.nvim_set_keymap("n", "<leader>o", ":NvimTreeToggle<CR>", {noremap = true, silent = true})


local tree_cb = require"nvim-tree.config".nvim_tree_callback

require'nvim-tree'.setup {
    ["actions.open_file"] = {
        quit_on_open = 1,
        -- when opening a file and a terminal is taking up one of the windows, ignore it for selection of where to put the new window.
        ["window_picker.exclude.buftype"] = {'terminal'}
    },
    -- closes neovim automatically when the tree is the last **WINDOW** in the view
    auto_close = true,
    -- opens the tree when changing/opening a new tab if the tree wasn't previously opened
    open_on_tab = true,
    -- update the focused file on `BufEnter`, un-collapses the folders recursively until it finds the file
    update_focused_file = {
        -- enables the feature
        enable = true,
    },
    view = { mappings = { list = {
        { key = {"<CR>", "<2-LeftMouse>"}, cb = tree_cb("edit") },  
        { key = "o", cb = tree_cb("system_open") },
    }}},
    filters = {
        dotfiles = true,
        custom = {".git", ".cache", "Icon\r", ".DS_Store"}
    }
}


