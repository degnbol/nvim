local g = vim.g

vim.o.termguicolors = true

g.nvim_tree_side = "left"
g.nvim_tree_width = 25
g.nvim_tree_ignore = {".git", ".cache"}
g.nvim_tree_auto_open = 0
g.nvim_tree_auto_close = 1
g.nvim_tree_quit_on_open = 1
g.nvim_tree_follow = 1
g.nvim_tree_indent_markers = 1
g.nvim_tree_hide_dotfiles = 1
g.nvim_tree_git_hl = 1
g.nvim_tree_highlight_opened_files = 1
g.nvim_tree_tab_open = 1
g.nvim_tree_allow_resize = 0

g.nvim_tree_show_icons = {git = 1, folders = 1, files = 1}

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

-- when opening a file and a terminal is taking up one of the windows, ignore it for selection of where to put the new window.
g.nvim_tree_window_picker_exclude = {buftype={'terminal'}}


-- Mappings

vim.api.nvim_set_keymap("n", "<leader>t", ":NvimTreeToggle<CR>", {noremap = true, silent = true})


function NvimTreeOSOpen()
  local lib = require "nvim-tree.lib"
  local node = lib.get_node_at_cursor()
  if node then
    vim.fn.jobstart("open '" .. node.absolute_path .. "' &", {detach = true})
  end
end


local tree_cb = require "nvim-tree.config".nvim_tree_callback
 
g.nvim_tree_bindings = {
    { key = "o", cb = ":lua require'tree'NvimTreeOSOpen()<CR>" },
}
