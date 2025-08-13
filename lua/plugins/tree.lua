-- tree file explorer to the left.
return {
    -- alt with more features:
    -- { "ms-jpq/chadtree", build="python3 -m chadtree deps", }
    -- alt that adds icons to the default netrw which isn't a project drawer.
    -- Combine with https://github.com/tpope/vim-vinegar
    -- "prichrd/netrw.nvim"
    -- a version of nvim-tree/nvim-tree with material icons:
    -- or consider a lua version: https://github.com/miversen33/netman.nvim
    {
        "nvim-tree/nvim-tree.lua",
        cmd = "NvimTreeToggle",
        dependencies = { 'nvim-tree/nvim-web-devicons' },
        init = function()
            -- disable netrw
            vim.g.loaded_netrw = 1
            vim.g.loaded_netrwPlugin = 1

            -- some issue with highlighting shell file logo
            -- vim.g.nvim_tree_highlight_opened_files = 1
            vim.g.nvim_tree_allow_resize = 0

            vim.keymap.set("n", "<leader>e", "<Cmd>NvimTreeToggle<CR>", { desc = "Explorer" })
        end,
        opts = {
            actions = { open_file = { quit_on_open = true } },
            -- opens the tree when changing/opening a new tab if the tree wasn't previously opened
            open_on_tab = true,
            -- update the focused file on `BufEnter`, un-collapses the folders recursively until it finds the file
            update_focused_file = { enable = true, },
            filters = {
                dotfiles = true,
                custom = { ".git", ".cache", "Icon\r", ".DS_Store" }
            },
            renderer = {
                indent_markers = { enable = true },
                icons = {
                    show = { git = true, folder = true, file = true },
                    glyphs = {
                        default = "",
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
                    },
                },
                highlight_git = true,
            },
            on_attach = function(bufnr)
                local api = require 'nvim-tree.api'

                local function nmap(lhs, rhs, desc)
                    vim.keymap.set('n', lhs, rhs, { desc = desc, buffer = bufnr, nowait = true })
                end

                nmap('<C-]>', api.tree.change_root_to_node, 'CD')
                -- important: don't use this since <C-e> is for scrolling.
                -- nmap('<C-e>', api.node.open.replace_tree_buffer, 'Open: In Place')
                nmap('<C-k>', api.node.show_info_popup, 'Info')
                nmap('<C-r>', api.fs.rename_sub, 'Rename: Omit Filename')
                nmap('<C-t>', api.node.open.tab, 'Open: New Tab')
                nmap('<C-v>', api.node.open.vertical, 'Open: Vertical Split')
                nmap('<C-x>', api.node.open.horizontal, 'Open: Horizontal Split')
                nmap('<BS>', api.node.navigate.parent_close, 'Close Directory')
                nmap('<CR>', api.node.open.edit, 'Open')
                nmap('<Tab>', api.node.open.preview, 'Open Preview')
                nmap('>', api.node.navigate.sibling.next, 'Next Sibling')
                nmap('<', api.node.navigate.sibling.prev, 'Previous Sibling')
                nmap('.', api.node.run.cmd, 'Run Command')
                nmap('-', api.tree.change_root_to_parent, 'Up')
                nmap('a', api.fs.create, 'Create')
                nmap('bmv', api.marks.bulk.move, 'Move Bookmarked')
                nmap('B', api.tree.toggle_no_buffer_filter, 'Toggle No Buffer')
                nmap('c', api.fs.copy.node, 'Copy')
                nmap('C', api.tree.toggle_git_clean_filter, 'Toggle Git Clean')
                nmap('[c', api.node.navigate.git.prev, 'Prev Git')
                nmap(']c', api.node.navigate.git.next, 'Next Git')
                nmap('d', api.fs.remove, 'Delete')
                nmap('D', api.fs.trash, 'Trash')
                nmap('E', api.tree.expand_all, 'Expand All')
                nmap('e', api.fs.rename_basename, 'Rename: Basename')
                nmap(']e', api.node.navigate.diagnostics.next, 'Next Diagnostic')
                nmap('[e', api.node.navigate.diagnostics.prev, 'Prev Diagnostic')
                nmap('F', api.live_filter.clear, 'Clean Filter')
                nmap('f', api.live_filter.start, 'Filter')
                nmap('g?', api.tree.toggle_help, 'Help')
                nmap('gy', api.fs.copy.absolute_path, 'Copy Absolute Path')
                nmap('H', api.tree.toggle_hidden_filter, 'Toggle Dotfiles')
                nmap('I', api.tree.toggle_gitignore_filter, 'Toggle Git Ignore')
                nmap('J', api.node.navigate.sibling.last, 'Last Sibling')
                nmap('K', api.node.navigate.sibling.first, 'First Sibling')
                nmap('m', api.marks.toggle, 'Toggle Bookmark')
                nmap('o', api.node.open.edit, 'Open')
                nmap('O', api.node.open.no_window_picker, 'Open: No Window Picker')
                nmap('p', api.fs.paste, 'Paste')
                nmap('P', api.node.navigate.parent, 'Parent Directory')
                nmap('q', api.tree.close, 'Close')
                nmap('r', api.fs.rename, 'Rename')
                nmap('R', api.tree.reload, 'Refresh')
                nmap('s', api.node.run.system, 'Run System')
                nmap('S', api.tree.search_node, 'Search')
                nmap('U', api.tree.toggle_custom_filter, 'Toggle Hidden')
                nmap('W', api.tree.collapse_all, 'Collapse')
                nmap('x', api.fs.cut, 'Cut')
                nmap('y', api.fs.copy.filename, 'Copy Name')
                nmap('Y', api.fs.copy.relative_path, 'Copy Relative Path')
                nmap('<2-LeftMouse>', api.node.open.edit, 'Open')
                nmap('<2-RightMouse>', api.tree.change_root_to_node, 'CD')
                nmap('<CR>', api.node.open.edit, 'Open')
                nmap('<2-LeftMouse>', api.node.open.edit, 'Open')
                nmap('o', api.node.run.system, 'Run System')
            end,
        }
    }
}
