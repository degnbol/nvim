local g = vim.g

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

-- open file in default program with o
local tree_cb = require"nvim-tree.config".nvim_tree_callback

require'nvim-tree'.setup {
    actions = {open_file = {quit_on_open = true}},
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
    },
    renderer={indent_markers={enable = true}},
}


