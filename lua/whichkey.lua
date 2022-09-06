local wk = require("which-key")


wk.setup {
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
        border = "none",
        position = "bottom",
        margin = {0, 0, 0, 0},
        padding = {0, 0, 0, 0}
    },
    layout = {
        height = {min = 4, max = 25}, -- min and max height of the columns
        width = {min = 20, max = 50}, -- min and max width of the columns
        spacing = 0 -- spacing between columns
    },
    ignore_missing = false, -- enable this to hide mappings for which you didn't specify a label
    hidden = {"<silent>", "<cmd>", "<Cmd>", "<CR>", "call", "lua", "^:", "^ "}, -- hide mapping boilerplate
    show_help = true, -- show help message on the command line when the popup is visible
    triggers = "auto" -- automatically setup triggers
    -- triggers = {"<leader>"} -- or specifiy a list manually
}

wk.register({
    ["<TAB>"] = {":BufferLineCycleNext<CR>", "next buffer"},
    ["<S-TAB>"] = {":BufferLineCyclePrev<CR>", "previous buffer"},
    ["<leader>"] = {
        ["<up>"] = "select increment (treesitter)",
        ["<down>"] = "select decrement (treesitter)",
        ["<S-up>"] = "select scope increment (treesitter)",
        ["/"] = {":CommentToggle<CR>", "(un)comment"}, -- see comment.lua
        ["1"] = {":BufferLineGoToBuffer 1<CR>", "buffer 1 (bufferline)"},
        ["2"] = {":BufferLineGoToBuffer 2<CR>", "buffer 2 (bufferline)"},
        ["3"] = {":BufferLineGoToBuffer 3<CR>", "buffer 3 (bufferline)"},
        ["4"] = {":BufferLineGoToBuffer 4<CR>", "..."},
        ["5"] = {":BufferLineGoToBuffer 5<CR>", "which_key_ignore"},
        ["6"] = {":BufferLineGoToBuffer 6<CR>", "which_key_ignore"},
        ["7"] = {":BufferLineGoToBuffer 7<CR>", "which_key_ignore"},
        ["8"] = {":BufferLineGoToBuffer 8<CR>", "which_key_ignore"},
        ["9"] = {":BufferLineGoToBuffer 9<CR>", "which_key_ignore"},
        ["?"] = "cheatsheet",
        ["<CR>"] = "kitty REPL",
        a = "swap arg (treesitter textobjects)",
        A = "swap arg back (treesitter textobjects)",
        -- currently just for latex but can be depend on filetype
        c = {
            name = "code (LSP)...",
            a = "code action (LSP)...",
        },
        d = {
            name = "defintion peek (treesitter textobjects + LSP)",
            f = "function",
            F = "class",
        },
        D = "type defintion (LSP)",
        e = {":NvimTreeToggle<CR>", "explorer"},
        E = "errors (LSP diagnostics)",
        -- telescope and dashboard mappings
        f = {
            name = "file",
            b = {":Telescope buffers<CR>", "buffers (telescope)"},
            f = {":Telescope find_files<CR>", "find (telescope)"},
            h = {":Telescope help_tags<CR>", "help tags (telescope)"},
            n = {":DashboardNewFile<CR>", "new (dashboard)"},
            o = {":Telescope oldfiles<CR>", "recent (telescope)"},
            w = {":Telescope live_grep<CR>", "find in files (telescope)"},
            W = {":Telescope grep_string<CR>", "find word under cursor (telescope)"},
        },
        h = {
            name = "git(signs)",
            b = "blame line",
            p = "preview hunk",
            r = "reset hunk",
            s = "stage hunk",
            u = "undo stage hunk",
        },
        -- <leader>K since K is regular defintion of word under cursor.
        K = {":DashWord<CR>", "Dash word"},
        l = {":noh<CR>", "clear highlights"},
        L = {
            name = "LaTeX",
            c = {":VimtexCompile<CR>", "compile"},
            e = {':lua require("nabla").popup()<CR>', "equation"},
        },
        m = {
            m = {":DashboardJumpMarks<CR>", "jump marks (dashboard)"}
        },
        p = {
            name = "paste after",
            ["#"] = {"<Plug>UnconditionalPasteCommentedAfter", "commented"},
            -- the `[ shouldn't be necessary but for some reason the cursor may 
            -- be moved to the first empty line of the paste, so we make sure 
            -- to move it to start of paste before autoindent.
            ["="] = {"<Plug>UnconditionalPasteLineAfter()`[=']", "autoindent"},
            b = {"<Plug>UnconditionalPasteBlockAfter", "blockwise"},
            c = {"<Plug>UnconditionalPasteCharAfter", "charwise"},
            C = {"<Plug>UnconditionalPasteCharCondensedAfter", "charwise condensed"},
            i = {"<Plug>UnconditionalPasteIndentedAfter", "indented"},
            j = {"<Plug>UnconditionalPasteJustJoinedAfter", "just joined"},
            l = {"<Plug>UnconditionalPasteLineAfter", "linewise"},
            n = {"<Plug>UnconditionalPasteInlinedAfter", "inlined"},
            -- switch to paste mode even though it is reported in :h nopaste as 
            -- obsolete. p=paste after. gq=format a motion. '[ and '] marks for 
            -- start and end of last edit (the paste). Depending on setting for 
            -- p the cursor will be at one or the other end after paste. 
            -- Default vim settings are cursor unchanged after paste. I here 
            -- assume it is moved after pasted text.
            q = {":set paste<CR>p:set nopaste<CR>`[gq']", "format"},
            -- ... there are many more to consider https://github.com/inkarkat/vim-UnconditionalPaste
        },
        P = {
            name = "paste before",
            ["#"] = {"<Plug>UnconditionalPasteCommentedBefore", "commented"},
            ["="] = {"<Plug>UnconditionalPasteLineBefore()`[=']", "autoindent"},
            b = {"<Plug>UnconditionalPasteBlockBefore", "blockwise"},
            c = {"<Plug>UnconditionalPasteCharBefore", "charwise"},
            C = {"<Plug>UnconditionalPasteCharCondensedBefore", "charwise condensed"},
            i = {"<Plug>UnconditionalPasteIndentedBefore", "indented"},
            j = {"<Plug>UnconditionalPasteJustJoinedBefore", "just joined"},
            l = {"<Plug>UnconditionalPasteLineBefore", "linewise"},
            n = {"<Plug>UnconditionalPasteInlinedBefore", "inlined"},
            q = {":set paste<CR>P:set nopaste<CR>`[gq']", "format"},
            -- ... there are many more to consider https://github.com/inkarkat/vim-UnconditionalPaste
        },
        r = {
            -- as long as there is only one function under r
            name = "rename...",
            n = "rename",
        },
        s = {
            name = "substitute OR session",
            -- substitute is an optional feature enabled from the substitute package where
            -- I can substitute e.g. all occurrences of a word in a paragraph with some new text by writing <leader>swip then the replacement text.
            l = {":SessionLoad<CR>", "load session"},
            s = {":SessionSave<CR>", "save session"},
        },
        t = {
            name = "toggle OR terminal",
            -- switch to/from Danish æøå and to insert mode, which is convenient.
            d = {'i<C-^>', "Danish (ctrl+^)"},
            l = {':silent HlSearchLensToggle<CR>', "HlSearchLens"},
            m = {':MarkdownToggle', "markdown preview"},
            s = {':ScrollbarToggle<CR>', "scrollbar"},
            t = {':let b:repl_id = input("window id: ")<CR>', "set REPL id"},
        },
        x = {':BufDel<CR>', "delete buffer"}, -- ojroques BufDel
    },
    c = {
        name = "change...",
        s = "surround...",
    },
    d = {
        name = "delete...",
        -- mini package
        s = "surround...",
    },
    g = {
        ["*"] = "search word under cursor flexibly", -- flexibly=ignore case and whole word
        ["#"] = "search word under cursor flexibly",
        ["/"] = "highlight last search",
        c = {
            name = "(un)comment motion",
            c = "line",
        },
        j = {function() require"trevj".format_at_cursor() end, "unjoin (trevj)"},
        -- uses % from andymass/vim-matchup which first jumps to container start, then visual, then container end, then core vim Join.
        -- Couldn't get it to work with whichkey probably since we need to remap %, so it is mapped in keymap.vim
        J = "join container",
        -- normally gp is p where cursor is moved at end. Since we do that by default, we can use it for whitepaste
        p = "whitepaste after",
        P = "whitepaste before",
        q = {
            name = "format motion",
            q = "line",
        },
        r = "references (LSP)",
    }
}, {mode='n'})

-- visual
wk.register({
    ["<leader>"] = {
        ["/"] = {":CommentToggle<CR>", "(un)comment"}, -- see comment.lua
    },
    ["<ScrollWheelUp>"] = "which_key_ignore",
    ["<ScrollWheelDown>"] = "which_key_ignore",
    ["."] = "increment",
    [","] = "decrement",
}, {mode='v'})

