-- pop-up to help with keybindings that have been started
return {
    {
        "mrjones2014/legendary.nvim",
        opts = { which_key = { auto_register = true } },
    },
    {
    "folke/which-key.nvim",
    event = "VeryLazy",
    -- doesn't dependend on legendary but we want to call it first
    dependencies={ "mrjones2014/legendary.nvim", },
    config=function()
--
local wk = require "which-key"

wk.setup {
    plugins = {
        marks = true, -- shows a list of your marks on ' and `
        registers = true, -- shows your registers on " in NORMAL or <C-r> in INSERT mode
        -- the presets plugin, adds help for a bunch of default keybindings in Neovim
        -- No actual key bindings are created
        spelling = {
            enabled = true, -- enabling this will show WhichKey when pressing z= to select spelling suggestions
            suggestions = 20 -- how many suggestions should be shown in the list?
        },
        presets = {
            -- adds help for operators like d, y, ... and registers them for motion / text object completion
            -- disabled since gU etc gives a window without wait which I don't need
            operators = false,
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
        height = {min = 3, max = 18}, -- min and max height of the columns
        width = {min = 15, max = 60}, -- min and max width of the columns
        spacing = 2 -- spacing between columns
    },
    ignore_missing = false, -- hide unlabeled?
    -- hide boilerplate
    hidden = {"<silent>", "<cmd>", "<Cmd>", "<CR>", "call", "lua", "^:", "^ "},
    show_help = true, -- show line at bottom with whichkey keybindings
    -- triggers = "auto", -- automatically setup triggers
    triggers = {"<leader>", "g", "y"}, -- or specify a list manually
    triggers_nowait = {}, 
}

wk.register({
    ["<S-CR>"] = {":lua focus_repl()<CR>", "Focus REPL"},
    ["<TAB>"] = {"<Cmd>BufferLineCycleNext<CR>", "next buffer"},
    ["<S-TAB>"] = {"<Cmd>BufferLineCyclePrev<CR>", "previous buffer"},
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
        ["4"] = {":BufferLineGoToBuffer 4<CR>", "…"},
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
            -- TODO: maybe only list Diffview options for diffview types?
            name = "code (LSP)…|choose (Diffview)…|color…",
            a = "code action (LSP)…|all (diffview)",
            b = "base (diffview)",
            d = {"dx", "none/delete (diffview). Use dx"},
            o = "ours (diffview)",
            s = {function()
                -- Load colorschemes. Using lazy keys didn't work
                -- https://www.reddit.com/r/neovim/comments/12tcx0b/attempt_at_adding_color_schemes_to_list_of/
                vim.api.nvim_exec_autocmds("User", { pattern = "ColorSchemeLoad" })
                -- simple attempt at only showing dark themes if we are in dark-mode
                if vim.o.background == "dark" then
                    vim.api.nvim_exec_autocmds("User", { pattern = "ColorSchemeLoadDark" })
                end
                require("telescope.builtin").colorscheme()
            end, "colorscheme"},
            t = "theirs (diffview)",
        },
        C = { "ga", "Character under cursor" }, -- set here because ga is replaced with https://github.com/junegunn/vim-easy-align 
        d = {
            name = "line diagnostics|definition peek",
            -- treesitter textobjects + LSP
            f = "function",
            F = "class",
            d = {":lua vim.diagnostic.open_float()<CR>", "Diagnostic on this line"},
        },
        D = "type defintion (LSP)",
        -- is overwritten by Diffview for its own buffer types so
        -- the description is correct here even though the command only indicates NvimTree
        e = {"<Cmd>NvimTreeToggle<CR>", "explorer (NvimTree|Diffview)"},
        E = "errors (LSP diagnostics)",
        -- telescope and dashboard mappings
        f = {
            name = "file",
            b = {"<Cmd>Telescope buffers<CR>", "buffers (telescope)"},
            f = {"<Cmd>Telescope find_files<CR>", "find (telescope)"},
            h = {"<Cmd>Telescope help_tags<CR>", "help tags (telescope)"},
            n = {"<Cmd>enew<CR>", "new file"},
            o = {"<Cmd>Telescope oldfiles<CR>", "recent (telescope)"},
            w = {"<Cmd>Telescope live_grep<CR>", "find in files (telescope)"},
            W = {"<Cmd>Telescope grep_string<CR>", "find word under cursor (telescope)"},
        },
        h = {":noh<CR>", "clear highlights"},
        j = { "<Plug>Join", "join (see splitjoin.lua)" },
        -- <leader>K since K is regular definition of word under cursor.
        -- K = {":DashWord<CR>", "Dash word"},
        l = {
            name = "LaTeX (vimtex|telescope-bibtex|nabla)",
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
            p = {"<Plug>UnconditionalPasteParagraphedAfter", "paragraphed"},
            q = {":set paste<CR>p:set nopaste<CR>`[gq']", "format"},
            s = {"<Plug>UnconditionalPasteSpacedAfter", "spaced"},
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
            p = {"<Plug>UnconditionalPasteParagraphedBefore", "paragraphed"},
            q = {":set paste<CR>P:set nopaste<CR>`[gq']", "format"},
            s = {"<Plug>UnconditionalPasteSpacedBefore", "spaced"},
            -- ... there are many more to consider https://github.com/inkarkat/vim-UnconditionalPaste
        },
        r = {
            -- as long as there is only one function under r
            name = "rename…", n = "rename",
        },
        s = {"<Plug>Split", "split (see splitjoin.lua)"},
        t = {
            name = "toggle|terminal",
            -- switch to/from Danish æøå and to insert mode, which is convenient.
            d = {'i<C-^>', "Danish (ctrl+^)"},
            j = {"<Plug>JoinToggle", "splitjoin"},
            l = {':silent HlSearchLensToggle<CR>', "HlSearchLens"},
            m = {':MarkdownToggle', "markdown preview"},
            s = {':ScrollbarToggle<CR>', "scrollbar"},
            T = {':let b:repl_id = input("window id: ")<CR>', "set REPL winid"},
            t = {':lua set_repl_last()<CR>', "set last window as REPL"},
        },
        x = { ':BufDel<CR>', "delete buffer" }, -- ojroques BufDel
    },
    c = {
        name = "change…",
        s = {
            name = "surrounding…",
            c = "command (vimtex)",
            e = "environment (vimtex)",
            -- e.g. change $ $ to equation environment by typing "equation" at prompt
            ['$'] = "math (vimtex)"
        },
        a = {
            name = "a(round)…",
            -- e.g. parenthesis in math
            d = "delimiter (vimtex)",
            ['$'] = "math (vimtex)",
            P = "section (vimtex)",
            i = "item (vimtex)", -- defined in ftplugin
            m = "math (vimtex)", -- changed in ftplugin
        },
        i = {
            name = "in(side)…",
            -- e.g. parenthesis in math
            d = "delimiter (vimtex)",
            ['$'] = "math (vimtex)",
            P = "section (vimtex)",
            i = "item (vimtex)", -- defined in ftplugin
            m = "math (vimtex)", -- changed in ftplugin
        },
    },
    d = {
        name = "delete…",
        -- small hack to remove excess whitespace.
        -- iw also captures whitespace under cursor.
        ['i '] = {"ciw <Esc>", "delete excess whitespace"},
        -- mini package
        s = {
            name = "surrounding…",
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
                name = "a(round)…",
                c = "comment line",
                C = "comments",
            },
            i = {
                name = "in(side)…",
                c = "comments",
            },
            o = "new under",
            O = "new above",
        },
        j = "go down line",
        J = "join simply",
        k = "go up line",
        -- normally gp is p where cursor is moved at end. Since we do that by default, we can use it for whitepaste.
        -- the plugin uses ,p and ,P by default but that slows down using , for ,; moving between f/t searching.
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
        -- see below and in plugin/keymaps.vim for more subversive
        ['ss'] = {"<plug>(SubversiveSubstituteWordRange)", "substitute word under cursor"}
    },
    t = {
        name = "toggle… (vimtex)",
        s = {
            name = "style…",
            c = "command", -- e.g. section
            d = "delimiter", -- e.g. with(out) \left
            -- same as d, but looks through g:vimtex_delim_toggle_mod_list in reverse
            D = "delimiter reverse",
            e = "environment",
            f = "fraction", -- toggle / <-> \frac
            ['$'] = "equation", -- inline vs display etc
        },
    },
    y = {
        name = "you…",
        o = {
            name = "toggle option (unimpaired)",
            h = "hlsearch",
            i = "ignorecase",
            l = "list",
            n = "number",
            r = "relativenumber",
            s = "spell",
            w = "wrap",
        }
    },
    ['['] = {
        name = "Previous…",
        d = "diagnostic (LSP)",
        h = {"<Cmd>Gitsigns prev_hunk<CR>", "hunk (gitsigns)"},
        y = "Change paste (Yoink)",
        x = "conflict (Diffview)",
        ['4'] = {"<Plug>(vimtex-[n)", "equation (vimtex)"}, -- without shift
        ['$'] = {"<Plug>(vimtex-[N)", "equation (vimtex)"},
        m = {"<Plug>(vimtex-[n)", "equation (vimtex)"}, -- without shift
        M = {"<Plug>(vimtex-[N)", "equation (vimtex)"},
        p = "put above, same indent (unimpaired)",
        P = "put above, same indent (unimpaired)",
    },
    [']'] = {
        name = "Next…",
        d = "diagnostic (LSP)",
        h = {"<Cmd>Gitsigns next_hunk<CR>", "hunk (gitsigns)"},
        y = "Change paste (Yoink)",
        x = "conflict (Diffview)",
        ['4'] = {"<Plug>(vimtex-]n)", "equation (vimtex)"}, -- without shift
        ['$'] = {"<Plug>(vimtex-]N)", "equation (vimtex)"},
        m = {"<Plug>(vimtex-]n)", "equation (vimtex)"}, -- without shift
        M = {"<Plug>(vimtex-]N)", "equation (vimtex)"},
        p = "put below, same indent (unimpaired)",
        P = "put below, same indent (unimpaired)",
    },
    ['='] = {
        p = "paste formatted after (unimpaired)",
        P = "paste formatted before (unimpaired)",
        s = {
            name = "toggle setting (unimpaired)",
            h = "hlsearch",
            i = "ignorecase",
            l = "list",
            n = "number",
            r = "relativenumber",
            s = "spell",
            w = "wrap",
        },
    },
}, {mode='n'})

wk.register({
    ["<leader>"] = {
        j = "join (revj)",
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
    g = {
        -- substitute is an optional feature enabled from the substitute package where
        -- I can substitute e.g. all occurrences of a word in a paragraph with some new text by writing <leader>Swip then the replacement text.
        -- example: gsiwip to replace all instances of the current word under the cursor that exist within the paragraph under the cursor. 
        -- example: gsl_ to replace all instances of the character under the cursor on the current line.
        -- example: gssip to replace the word under cursor in the current paragraph. Matches complete words so is different from <leader>siwip
        -- See normal gss mapping above.
        s = { "<plug>(SubversiveSubstituteRange)", "substitute motion in motion" },
    },
}, {mode={'n', 'v'}})

wk.register({
    ['\\'] = {
        f = {'<Plug>(leap-forward-to)', "forward to (leap)"},
        F = {'<Plug>(leap-backward-to)', "backward to (leap)"},
        t = {'<Plug>(leap-forward-till)', "forward till (leap)"},
        T = {'<Plug>(leap-backward-till)', "backward till (leap)"},
    },
}, {mode={'n', 'x', 'o'}})

wk.register({
    ["<C-l>"] = {":lua require'telescope.builtin'.keymaps(require('telescope.themes').get_ivy())<CR>", "Keymaps"},
})

end}


}
