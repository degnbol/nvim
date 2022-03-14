local whichkey = require("which-key")


whichkey.setup {
    plugins = {
        marks = true, -- shows a list of your marks on ' and `
        registers = true, -- shows your registers on " in NORMAL or <C-r> in INSERT mode
        -- the presets plugin, adds help for a bunch of default keybindings in Neovim
        -- No actual key bindings are created
        spelling = {
            enabled = false, -- enabling this will show WhichKey when pressing z= to select spelling suggestions
            suggestions = 20 -- how many suggestions should be shown in the list?
        },
        presets = {
            operators = true, -- adds help for operators like d, y, ... and registers them for motion / text object completion
            motions = true, -- adds help for motions
            text_objects = true, -- help for text objects triggered after entering an operator
            windows = true, -- default bindings on <c-w>
            nav = true, -- misc bindings to work with windows
            z = true, -- bindings for folds, spelling and others prefixed with z
            g = true -- bindings for prefixed with g
        }
    },
    -- add operators that will trigger motion and text object completion
    -- to enable all native operators, set the preset / operators plugin above
    operators = {gc = "Comments"},
    icons = {
        breadcrumb = "»", -- symbol used in the command line area that shows your active key combo
        separator = "➜", -- symbol used between a key and it's label
        group = "+" -- symbol prepended to a group
    },
    window = {
        border = "none", -- none, single, double, shadow
        position = "bottom", -- bottom, top
        margin = {1, 0, 1, 0}, -- extra window margin [top, right, bottom, left]
        padding = {2, 2, 2, 2} -- extra window padding [top, right, bottom, left]
    },
    layout = {
        height = {min = 4, max = 25}, -- min and max height of the columns
        width = {min = 20, max = 50}, -- min and max width of the columns
        spacing = 3 -- spacing between columns
    },
    ignore_missing = false, -- enable this to hide mappings for which you didn't specify a label
    hidden = {"<silent>", "<cmd>", "<Cmd>", "<CR>", "call", "lua", "^:", "^ "}, -- hide mapping boilerplate
    show_help = true, -- show help message on the command line when the popup is visible
    triggers = "auto" -- automatically setup triggers
    -- triggers = {"<leader>"} -- or specifiy a list manually
}


-- name all kinds of shortcuts
whichkey.register{
    ["<leader>"] = {
        ["/"] = "(un)comment",
        ["1"] = "buffer 1 (bufferline)",
        ["2"] = "buffer 2 (bufferline)",
        ["3"] = "buffer 3 (bufferline)",
        ["4"] = "...",
        ["5"] = "which_key_ignore",
        ["6"] = "which_key_ignore",
        ["7"] = "which_key_ignore",
        ["8"] = "which_key_ignore",
        ["9"] = "which_key_ignore",
        ["?"] = "cheatsheet",
        ["<CR>"] = {":lua kittyWindow()<CR>", "kitty REPL"},
        a = "swap arg (treesitter textobjects)",
        A = "swap arg back (treesitter textobjects)",
        -- currently just for latex but can be depend on filetype
        b = {":VimtexCompile<CR>", "build"},
        d = "defintion? (treesitter textobjects)",
        D = "defintion? (treesitter textobjects)",
        e = {":NvimTreeToggle<CR>", "explorer"},
        -- telescope and dashboard mappings
        f = {
            name = "file",
            b = {":lua require('telescope.builtin').buffers()<CR>", "buffers (telescope)"},
            f = {":lua require('telescope.builtin').find_files()<CR>", "find (telescope)"},
            h = {":lua require('telescope.builtin').help_tags()<CR>", "help tags (telescope)"},
            n = {":DashboardNewFile<CR>", "new (dashboard)"},
            o = {":lua require('telescope.builtin').oldfiles()<CR>", "recent (telescope)"},
            w = {":Telescope live_grep<CR>", "live grep (telescope)"},
        },
        h = {
            name = "git(signs)",
            b = "blame line",
            p = "preview hunk",
            r = "reset hunk",
            s = "stage hunk",
            u = "undo stage hunk",
        },
        m = {
            m = {":DashboardJumpMarks<CR>", "jump marks (dashboard)"}
        },
        s = {
            name = "substitute OR session",
            l = {":SessionLoad<CR>", "load session"},
            s = {":SessionSave<CR>", "save session"},
        },
        t = {
            name = "terminal",
            t = {':let b:repl_id = input("window id: ")<CR>', "set REPL id"},
        },
        x = "delete buffer",
    },
    d = {
        name = "delete",
        s = {
            -- from the surround package
            name = "delete surround"
        },
    },
    g = {
        -- from UnconditionalPaste
        ["#"] = "commented paste",
        [","] = "comma paste",
        ["="] = "expression paste",
        [">"] = "indent paste",
        ["\\"] = "escape paste",
        ["["] = "indented paste",
        ["]"] = "indented paste",
        b = "block paste",
        B = "jagged paste",
        -- without UnconditionalPaste plugin it will (un)comment motion
        c = "char paste OR (un)comment",
        C = "char condensed paste",
        h = "combinatorial paste",
        H = "recombinatorial paste",
        l = "line paste",
        q = {
            name = "queried OR delimited paste",
            b = {
                name = "delimited paste",
                p = "after",
                P = "before",
            },
            g = {
                name = "queried joined paste",
                p = "after",
                P = "before",
            },
            p = "queried after",
            P = "queried before",
            q = "format line",
        },
        Q = "requeried OR redelimited paste",
        -- grep paste from UnconditionalPaste, refactor from treesitter
        r = "grep paste OR refactor",
        R = "regrep paste",
        s = {
            name = "spaced paste",
            p = "after",
            P = "before",
        },
        S = {
            name = "paragraphed paste",
            p = "after",
            P = "before",
        },
    }
}
