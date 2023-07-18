return {
    -- yank history that you can cycle with c-n and c-p
    "svermeulen/vim-yoink",
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
        -- set keybindings myself in whichkey.lua
        vim.g.UnconditionalPaste_no_mappings = true
    end},
    -- paste with multiple empty lines around contents reduced to single empty lines.
    -- TODO: the cmd isn't slow, it's just that there is a wait since other keys may be pressed.
    {"AndrewRadev/whitespaste.vim", init = function ()
        vim.cmd [[
        let g:whitespaste_before_mapping = 'gP'
        let g:whitespaste_after_mapping  = 'gp'
        ]]
    end,},
}

