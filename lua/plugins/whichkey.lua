-- pop-up to help with keybindings that have been started
return {
    {
    "folke/which-key.nvim",
    -- doesn't dependend on legendary but we want to call it first
    dependencies={ "mrjones2014/legendary.nvim", },
    config=function()
--
require('legendary').setup({ which_key = { auto_register = true } })


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
            g = true, -- bindings for prefixed with g
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
        spacing = 2 -- spacing between columns
    },
    ignore_missing = false, -- enable this to hide mappings for which you didn't specify a label
    hidden = {"<silent>", "<cmd>", "<Cmd>", "<CR>", "call", "lua", "^:", "^ "}, -- hide mapping boilerplate
    show_help = true, -- show help message on the command line when the popup is visible
    -- triggers = "auto", -- automatically setup triggers
    triggers = {"<leader>", "c", "cr", "g"}, -- or specifiy a list manually
    triggers_nowait = { "cr", }, -- doesn't work
}

wk.register({
    ["<S-CR>"] = {":lua focus_repl()<CR>", "Focus REPL"},
    ["<TAB>"] = {":BufferLineCycleNext<CR>", "next buffer"},
    ["<S-TAB>"] = {":BufferLineCyclePrev<CR>", "previous buffer"},
    ["<leader>"] = {
        ["<leader>"] = {
            s = {"<cmd>source ~/.config/nvim/after/plugin/luasnip.lua<CR>", "reload snippets"},
        },
        ["<up>"] = "select increment (treesitter)",
        ["<down>"] = "select decrement (treesitter)",
        ["<S-up>"] = "select scope increment (treesitter)",
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
        b = {
            name = "bibliography (papis)",
            o = {"<>", "open"},
        },
        c = {
            -- works with both since Diffview only overwrites keybindings for their own buffer types
            -- For choosing none, use dx (delete conflict)
            name = "code (LSP) || choose (Diffview)...",
            a = "code action (LSP)... || all (diffview)",
            b = "base (diffview)",
            d = {"dx", "none/delete (diffview). Use dx"},
            o = "ours (diffview)",
            t = "theirs (diffview)",
        },
        C = {"ga", "Character under cursor"}, -- set here because ga is replaced with https://github.com/junegunn/vim-easy-align 
        d = {
            name = "line diagnostics || definition peek",
            -- treesitter textobjects + LSP
            f = "function",
            F = "class",
            d = {":lua vim.diagnostic.open_float()<CR>", "Diagnostic on this line"},
        },
        D = "type defintion (LSP)",
        -- is overwritten by Diffview for its own buffer types so
        -- the description is correct here even though the command only indicates NvimTree
        e = {":NvimTreeToggle<CR>", "explorer (NvimTree || Diffview)"},
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
        h = {":noh<CR>", "clear highlights"},
        -- <leader>K since K is regular defintion of word under cursor.
        j = {
            name = "to multiline (revj)",
            j = "line",
        },
        K = {":DashWord<CR>", "Dash word"},
        l = {
            name = "LaTeX (vimtex || telescope-bibtex || nabla)",
            a = "context menu",
            -- align table see ftplugin/tex.lua. ma ... `a to not move cursor.
            A = {"ma<plug>AlignTable<CR>`a", "Align table"},
            b = {":Telescope bibtex<CR>", "insert citation (bibtex)"},
            c = "clean",
            C = "clean full",
            e = "errors",
            E = {':lua require("nabla").popup()<CR>', "equation"},
            g = "status",
            G = "status all",
            i = "info",
            -- j for jump
            j = {"<plug>TexJumpPre", "goto/from preamble (table)"},
            I = "info full",
            k = "stop",
            K = "stop all",
            l = "compile",
            L = "compile selected",
            m = "imaps list",
            o = "compile output",
            -- p and P will be for pasting and reformatting table
            q = "log",
            s = "toggle main",
            t = "TOC open",
            T = "TOC toggle",
            u = {"<plug>Latex2Unicode", "latex2unicode"},
            U = {"<plug>Unicode2Latex", "unicode2latex"},
            v = "view",
            x = "reload",
            X = "reload state",
            y = {"<plug>YankTable", "yank table as tsv"},
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
            name = "substitute || session",
            -- substitute is an optional feature enabled from the substitute package where
            -- I can substitute e.g. all occurrences of a word in a paragraph with some new text by writing <leader>swip then the replacement text.
            l = {":SessionLoad<CR>", "load session"},
            s = {":SessionSave<CR>", "save session"},
        },
        t = {
            name = "toggle || terminal",
            -- switch to/from Danish æøå and to insert mode, which is convenient.
            d = {'i<C-^>', "Danish (ctrl+^)"},
            l = {':silent HlSearchLensToggle<CR>', "HlSearchLens"},
            m = {':MarkdownToggle', "markdown preview"},
            s = {':ScrollbarToggle<CR>', "scrollbar"},
            T = {':let b:repl_id = input("window id: ")<CR>', "set REPL winid"},
            t = {':lua set_repl_last()<CR>', "set last window as REPL"},
        },
        x = {':BufDel<CR>', "delete buffer"}, -- ojroques BufDel
    },
    c = {
        name = "change...",
        s = {
            name = "surrounding...",
            c = "command (vimtex)",
            e = "environment (vimtex)",
            -- e.g. change $ $ to equation environment by typing "equation" at prompt
            ['$'] = "math (vimtex)"
        },
        a = {
            name = "a(round)...",
            -- e.g. parenthesis in math
            d = "delimiter (vimtex)",
            ['$'] = "math (vimtex)",
            P = "section (vimtex)",
            i = "item (vimtex)", -- defined in ftplugin
            m = "math (vimtex)", -- changed in ftplugin
        },
        i = {
            name = "in(side)...",
            -- e.g. parenthesis in math
            d = "delimiter (vimtex)",
            ['$'] = "math (vimtex)",
            P = "section (vimtex)",
            i = "item (vimtex)", -- defined in ftplugin
            m = "math (vimtex)", -- changed in ftplugin
        },
    },
    cr = {
            name = "coerce (casing)",
            c = "camelCase",
            s = "snake_case",
            m = "MixedCase",
            ["<space>"] = "mixed case",
        },
    d = {
        name = "delete...",
        -- mini package
        s = {
            name = "surrounding...",
            -- see :h vimtex-default-mappings
            c = "command (vimtex)",
            -- includes \left and \right if connected to e.g. ( and )
            d = "delimiter (vimtex)",
            e = "environment (vimtex)",
        },
    },
    g = {
        ["*"] = "search word under cursor flexibly", -- flexibly=ignore case and whole word
        ["#"] = "search word under cursor flexibly",
        ["/"] = "highlight last search",
        a = "show char info",
        c = {
            name = "(un)comment motion",
            c = "line",
            A = "new EOL",
            a = {
                name = "a(round)...",
                c = "comment line",
                C = "comments",
            },
            i = {
                name = "in(side)...",
                c = "comments",
            },
            o = "new under",
            O = "new above",
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
        w = {
            -- use the default vim formatter instead of whichever is set for a language.
            name = "vim format motion",
            w = "line",
        },
        r = "references (LSP)",
    },
    t = {
        name = "toggle... (vimtex)",
        s = {
            name = "style...",
            c = "command", -- e.g. section
            d = "delimiter", -- e.g. with(out) \left
            -- same as d, but looks through g:vimtex_delim_toggle_mod_list in reverse
            D = "delimiter reverse",
            e = "environment",
            f = "fraction", -- toggle / <-> \frac
            ['$'] = "equation", -- inline vs display etc
        },
    },
    ['['] = {
        name = "Previous...",
        d = "diagnostic (LSP)",
        h = {"<Cmd>Gitsigns prev_hunk<CR>", "hunk (gitsigns)"},
        y = "Change paste (Yoink)",
        x = "conflict (Diffview)",
        ['4'] = {"<Plug>(vimtex-[n)", "equation (vimtex)"}, -- without shift
        ['$'] = {"<Plug>(vimtex-[N)", "equation (vimtex)"},
        m = {"<Plug>(vimtex-[n)", "equation (vimtex)"}, -- without shift
        M = {"<Plug>(vimtex-[N)", "equation (vimtex)"},
    },
    [']'] = {
        name = "Next...",
        d = "diagnostic (LSP)",
        h = {"<Cmd>Gitsigns next_hunk<CR>", "hunk (gitsigns)"},
        y = "Change paste (Yoink)",
        x = "conflict (Diffview)",
        ['4'] = {"<Plug>(vimtex-]n)", "equation (vimtex)"}, -- without shift
        ['$'] = {"<Plug>(vimtex-]N)", "equation (vimtex)"},
        m = {"<Plug>(vimtex-]n)", "equation (vimtex)"}, -- without shift
        M = {"<Plug>(vimtex-]N)", "equation (vimtex)"},
    },
}, {mode='n'})

wk.register({
    ["<leader>"] = {
        l = {
            name = "latex",
            u = {"<plug>Latex2Unicode_visual", "latex2unicode"},
            U = {"<plug>Unicode2Latex_visual", "unicode2latex"},
        },
    },
    ["<ScrollWheelUp>"] = "which_key_ignore",
    ["<ScrollWheelDown>"] = "which_key_ignore",
    ["<C-space>"] = "increment (TS)",
    ["<C-backspace>"] = "decrement (TS)",
    ["<C-s>"] = "increment scope (TS)",
    -- see lua/plugin/treesitter.lua
    -- . and , kinda works like incr, decrement but we already have c-space and c-backspace
    ["."] = "smart (TS subj)",
    [","] = "prev (TS subj)",
    ["a;"] = "a container",
    ["i;"] = "in container",
}, {mode='v'})

wk.register({
    ["<leader>"] = {
        g = {
            name = "git",
            b = {"<Cmd>Gitsigns blame_line<CR>", "blame line (gitsigns)"},
            -- hide untracked files with -uno.
            -- hide gitsigns' file explorer with DiffviewToggleFiles (unhide with <leader>e like NvimTreeToggle)
            -- Open during merge or rebase should show conflicts nicer automatically.
            -- Isn't allowed for ranges (shows error), but keep it here to have all git stuff one place and to make it explicit if typed for whatever reason.
            d = {"<Cmd>DiffviewOpen -uno<CR>:DiffviewToggleFiles<CR>", "open diffview"},
            q = {"<Cmd>DiffviewClose<CR>", "quit (Diffview)"},
            -- pretty cool: history of a specific range of code
            h = {"<Cmd>DiffviewFileHistory<CR>", "history (Diffview)"},
            p = {"<Cmd>Gitsigns preview_hunk<CR>", "preview hunk (gitsigns)"},
            r = {"<Cmd>Gitsigns reset_hunk<CR>", "reset hunk (gitsigns)"},
            s = {"<Cmd>Gitsigns stage_hunk<CR>", "stage hunk (gitsigns)"},
            u = {"<Cmd>Gitsigns undo_stage_hunk<CR>", "undo stage hunk (gitsigns)"},
        },
    },
    ['\\'] = {
        f = "forward to (leap)",
        F = "backward to (leap)",
        t = "forward till (leap)",
        T = "backward till (leap)",
    },
}, {mode={'n', 'v'}})

wk.register({
    ["<C-l>"] = {":lua require'telescope.builtin'.keymaps(require('telescope.themes').get_ivy())<CR>", "Keymaps"},
})

end}


}
