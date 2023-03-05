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
            highlight = {
                -- enabled = true,         -- highlight yanked text
                -- higroup = "IncSearch",  -- highlight group of yanked text
                timeout = 200,         -- timeout for clearing the highlight
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
    -- lots of ways to paste using g{c,C,l,b}{,i}{p,P} and many others
    {"inkarkat/vim-UnconditionalPaste", dependencies={'inkarkat/vim-ingo-library'}, config=function()
        -- set keybindings myself in whichkey.lua
        vim.g.UnconditionalPaste_no_mappings = true
    end},
    -- paste without empty newlines
    {"AndrewRadev/whitespaste.vim", config=function()
        -- so p and P will paste whatever is in register, while gp and gP will paste 
        -- with fewer empty lines.
        vim.g.whitespaste_after_mapping  = 'gp'
        vim.g.whitespaste_before_mapping = 'gP'
    end},
}
