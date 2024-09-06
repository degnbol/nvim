-- TODO: p = default
-- gp to paste after.

return {
    -- add substitution functions to e.g. replace a word with clipboard content by writing siw
    {
        "svermeulen/vim-subversive",
        init = function ()
            vim.keymap.set({"n", "x"}, "s", "<Plug>(SubversiveSubstitute)", { desc="Substitute" })
            vim.keymap.set("n", "ss", "<Plug>(SubversiveSubstituteLine)", { desc="Substitute line" })
            vim.keymap.set("n", "S", "<Plug>(SubversiveSubstituteToEndOfLine)", { desc="Substitute to EOL" })
            -- replace default useless "Sleep"
            -- substitute is an optional feature enabled from the substitute package where
            -- I can substitute e.g. all occurrences of a word in a paragraph with some new text by writing <leader>Swip then the replacement text.
            -- example: gsiwip to replace all instances of the current word under the cursor that exist within the paragraph under the cursor. 
            -- example: gsl_ to replace all instances of the character under the cursor on the current line.
            -- example: gssip to replace the word under cursor in the current paragraph. Matches complete words so is different from <leader>siwip
            -- See normal gss mapping above.
            vim.keymap.set({"n", "x"}, "gs", "<Plug>(SubversiveSubstituteRange)", { desc="Substitute motion in motion" })
            vim.keymap.set({"n", "x"}, "gss", "<Plug>(SubversiveSubstituteWordRange)", { desc="Substitute word under cursor" })
            -- Is integrated with yoink to track clipboard history
            vim.keymap.set("x", "p", "<plug>(SubversiveSubstitute)", { desc="Substitute" })
            vim.keymap.set("x", "P", "<plug>(SubversiveSubstitute)", { desc="Substitute" })
        end,
    },
    -- yank history that you can cycle. I chose non-default [p (and ]p)
    {
        "svermeulen/vim-yoink",
        config = function ()
            -- yoink integration with cutlass
            vim.g.yoinkIncludeDeleteOperations = 1
            -- add yanks to numbered register
            vim.g.yoinkSyncNumberedRegisters = 1
            -- move cursor to end instead of start of multi-line paste
            vim.g.yoinkMoveCursorToEndOfPaste = 0
            -- preserve yank between neovim sessions
            vim.g.yoinkSavePersistently = 1

            vim.keymap.set('n', 'p', "<plug>(YoinkPaste_p)", { desc="Paste below with history" })
            vim.keymap.set('n', 'P', "<plug>(YoinkPaste_P)", { desc="Paste above with history" })
            -- whitepaste gets preference
            -- vim.keymap.set('n', 'gp', "<plug>(YoinkPaste_gp)", { desc="Paste below with history" })
            -- vim.keymap.set('n', 'gP', "<plug>(YoinkPaste_gP)", { desc="Paste above with history" })
            vim.keymap.set({'n', 'x'}, 'y', "<plug>(YoinkYankPreserveCursorPosition)", { desc="Yoink" })
            vim.keymap.set("n", "[p", "<Plug>(YoinkPostPasteSwapBack)", { desc="Swap paste" })
            vim.keymap.set("n", "]p", "<Plug>(YoinkPostPasteSwapForward)", { desc="Swap paste forward" })
            vim.keymap.set("n", "[y", "<Plug>(YoinkRotateBack)", { desc="Swap yank" })
            vim.keymap.set("n", "]y", "<Plug>(YoinkRotateForward)", { desc="Swap yank forward" })
            vim.keymap.set("n", "<C-=>", "<Plug>(YoinkPostPasteToggleFormat)", { desc="Toggle pasted indent" })
            -- esc, reindent back til start of paste, jump back to where we 
            -- just jumped from (`] has been changed by the reindent), then 
            -- insert mode, then <C-F> to reindent current line without leaving 
            -- insert mode (F for flush with line above). This should mean we 
            -- can handle both a charwise and a linewise paste.
            vim.keymap.set("i", "<C-=>", [=[<Esc>=`[``a<C-f>]=], { desc="Toggle pasted indent" })
        end,
    },
    -- c(hange), d(elete) no longer copies, remapped in keymapping file so x will cut. Since we have added backspace and delete button support in normal mode there is no need for default x behavior
    {"gbprod/cutlass.nvim",
    opts={
        cut_key='x',
        override_del = true, -- true -> del key will put in blackhole register}
    }},
    -- yank in tmux and over ssh
    {
        'ibhagwan/smartyank.nvim',
        enabled = false, -- breaks blockwise paste
        opts = {
            -- same as
            -- api.nvim_create_autocmd({"TextYankPost"}, { callback=function() vim.highlight.on_yank{on_visual=false} end })
            highlight = {
                -- enabled = true,         -- highlight yanked text
                -- higroup = "IncSearch",  -- highlight group of yanked text
                timeout = 100,         -- timeout for clearing the highlight
            },
            -- clipboard = {
            --   enabled = true
            -- },
            -- tmux = {
            --   enabled = true,
            --   -- remove `-w` to disable copy to host client's clipboard
            --   cmd = { 'tmux', 'set-buffer', '-w' }
            -- },
            -- osc52 = {
            --   enabled = true,
            --   ssh_only = true,        -- false to OSC52 yank also in local sessions
            --   silent = false,         -- true to disable the "n chars copied" echo
            --   echo_hl = "Directory",  -- highlight group of the OSC52 echo message
            -- }
        }
    },
    -- lots of ways to paste
    {"inkarkat/vim-UnconditionalPaste",
    event = "VeryLazy",
    dependencies={'inkarkat/vim-ingo-library'}, init=function()
        -- use <leader>p ... and <leader>P ... instead of the harder to remember gbp etc.
        vim.g.UnconditionalPaste_no_mappings = true
        function mapp(char, desc)
            vim.keymap.set("n", "<leader>P" .. char, "<Plug>UnconditionalPaste" .. desc .. "Before", { desc=desc })
            vim.keymap.set("n", "<leader>p" .. char, "<Plug>UnconditionalPaste" .. desc .. "After", { desc=desc })
        end
        mapp("#", "Commented")
        mapp("b", "Block")
        mapp("c", "Char")
        mapp("C", "CharCondensed")
        mapp("i", "Indented")
        mapp("j", "JustJoined")
        mapp("l", "Line")
        mapp("n", "Inlined")
        mapp("p", "Paragraphed")
        mapp("s", "Spaced")
        -- ... there are many more to consider https://github.com/inkarkat/vim-UnconditionalPaste

        -- special additions:
        -- the `[ shouldn't be necessary but for some reason the cursor may 
        -- be moved to the first empty line of the paste, so we make sure 
        -- to move it to start of paste before autoindent.
        vim.keymap.set("n", "<leader>P=", "<Plug>UnconditionalPasteLineBefore()`[=`]", { desc="Autoindent" })
        vim.keymap.set("n", "<leader>p=", "<Plug>UnconditionalPasteLineAfter()`[=`]", { desc="Autoindent" })
        -- Next four maps are not really unimpaired specific except they take up some of the same letters.
        -- switch to paste mode temp to disable various disturbances, see :h 
        -- nopaste. It is obsolete for regular cmd+v pasting.
        -- gq=format a motion. `[ and `] marks for 
        -- start and end of last edit (the paste).
        vim.keymap.set("n", "<leader>Pq", ":set paste<CR>P:set nopaste<CR>`[gq`]", { desc="Format" })
        vim.keymap.set("n", "<leader>pq", ":set paste<CR>p:set nopaste<CR>`[gq`]", { desc="Format" })
        -- default vim formatting, ignoring any formatfunc option
        vim.keymap.set("n", "<leader>Pw", ":set paste<CR>P:set nopaste<CR>`[gw`]", { desc="Format default" })
        vim.keymap.set("n", "<leader>pw", ":set paste<CR>p:set nopaste<CR>`[gw`]", { desc="Format default" })
    end},
    -- paste with multiple empty lines around contents reduced to single empty lines.
    {"AndrewRadev/whitespaste.vim",
    init = function ()
        -- normally gp is p where cursor is moved at end. Since we do that by default, we can use it for whitepaste.
        -- the plugin uses ,p and ,P by default but that slows down using , for ,; moving between f/t searching.
        vim.api.nvim_set_var("whitespaste_before_mapping", 'gP')
        vim.api.nvim_set_var("whitespaste_after_mapping",  'gp')

        vim.schedule(function ()
            require "utils/keymap"
            set_keymap_desc('n', 'gP', "Put before + trim empty")
            set_keymap_desc('n', 'gp', "Put after + trim empty")
        end)
    end,},
}

