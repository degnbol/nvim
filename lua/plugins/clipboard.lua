return {
    -- yank history that you can cycle with c-n and c-p
    {
        "svermeulen/vim-yoink",
    },
    -- c(hange), d(elete) no longer copies, remapped in keymapping file so x will cut. Since we have added backspace and delete button support in normal mode there is no need for default x behavior
    {"gbprod/cutlass.nvim", opts={
        cut_key='x',
        override_del = true, -- true -> del key will put in blackhole register}
    }},
    -- yank in tmux and over ssh
    {
        'ibhagwan/smartyank.nvim',
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
    {"AndrewRadev/whitespaste.vim", init = function ()
        -- normally gp is p where cursor is moved at end. Since we do that by default, we can use it for whitepaste.
        -- the plugin uses ,p and ,P by default but that slows down using , for ,; moving between f/t searching.
        vim.api.nvim_set_var("whitespaste_before_mapping", 'gP')
        vim.api.nvim_set_var("whitespaste_after_mapping",  'gp')

        vim.defer_fn(function ()
            require "utils/keymap"
            set_keymap_desc('n', 'gP', "Put before + trim empty")
            set_keymap_desc('n', 'gp', "Put after + trim empty")
        end, 0)
    end,},
}

