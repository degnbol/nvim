return {
    -- more languages supported but sometimes does nothing when trevj does the split
    {
        "Wansmer/treesj",
        dependencies = {
            'nvim-treesitter/nvim-treesitter',
            -- make sure is loaded for the init to work
            'AckslD/nvim-trevJ.lua',
            'AckslD/nvim-revJ.lua',
            'AndrewRadev/splitjoin.vim',
            'andymass/vim-matchup',
        },
        opts = {
            -- set only for supported filetypes
            use_default_keymaps = false
        },
        init = function ()
            -- but always map toggle since the other plugins doesn't implement it
            vim.keymap.set("n", "<Plug>JoinToggle", "<Cmd>TSJToggle<CR>")

            -- use alt plugins for specific filetypes
            local fts_split = { splitjoin = { }, trevj = { } }
            local fts_join  = { splitjoin = {}, matchup = {}, }
            local grp = vim.api.nvim_create_augroup("splitjoin", {clear=true})

            -- call nmap buffer thru vimscript since nowait isn't implemented in lua API
            function nowait(lhs, rhs)
                vim.cmd('nnoremap <buffer><nowait> ' .. lhs .. ' ' .. rhs)
            end
            function nowaitre(lhs, rhs)
                vim.cmd('nmap <buffer><nowait> ' .. lhs .. ' ' .. rhs)
            end
            
            vim.api.nvim_create_autocmd("Filetype", {
                pattern = "*", group = grp,
                callback = function ()
                    if fts_split['splitjoin'][vim.bo.filetype] then
                        nowait("<Plug>Split", "<Cmd>SplitjoinSplit<CR>")
                    elseif fts_split['trevj'][vim.bo.filetype] then
                        nowait("<Plug>Split", ":lua require'trevj'.format_at_cursor()<CR>")
                    else
                        nowait("<Plug>Split", "<Cmd>TSJSplit<CR>")
                    end
                    
                    if fts_join['splitjoin'][vim.bo.filetype] then
                        nowait("<Plug>Join", "<Cmd>SplitjoinJoin<CR>")
                    elseif fts_join['matchup'][vim.bo.filetype] then
                        -- Uses a% from https://github.com/andymass/vim-matchup which first jumps to container 
                        nowaitre("<Plug>Join", "va%J")
                    else
                        nowait("<Plug>Join", "<Cmd>TSJJoin<CR>")
                    end
                end
            })
        end
    },
    {"AckslD/nvim-trevJ.lua",
    lazy = true,
    config=function()
        local make_default_opts = function()
            return { final_separator = ",", final_end_line = true, skip = { } }
        end

        local make_no_final_sep_opts = function()
            return { final_separator = false, final_end_line = true, }
        end

        require'trevj'.setup {
            containers = {
                julia = {
                    parameters = make_default_opts(),
                    argument_list = make_default_opts(),
                    list = make_default_opts(),
                    tuple = make_default_opts(),
                    dictionary = make_default_opts(),
                    set = make_default_opts(),
                    list_comprehension = make_no_final_sep_opts(),
                    generator_expression = make_no_final_sep_opts(),
                    dictionary_comprehension = make_no_final_sep_opts(),
                }
            }
        }
    end},
    {
        "AndrewRadev/splitjoin.vim",
        lazy = true,
        init = function ()
            -- remove default mappings
            vim.g.splitjoin_split_mapping = ''
            vim.g.splitjoin_join_mapping = ''
        end
    },
    -- still using it since it's the only one with visual version and motion
    {"AckslD/nvim-revJ.lua", lazy = true, dependencies={'kana/vim-textobj-user', 'sgur/vim-textobj-parameter'}, opts={
        brackets = {first = '([{<', last = ')]}>'}, -- brackets to consider surrounding arguments
        new_line_before_last_bracket = true, -- add new line between last argument and last bracket (only if no last seperator)
        add_seperator_for_last_parameter = true, -- if a seperator should be added if not present after last parameter
        enable_default_keymaps = false,
        keymaps = {
            operator = 'gS', -- for operator (+motion)
            line = 'gSS', -- for formatting current line
            visual = '<leader>s', -- for formatting visual selection
        },
        parameter_mapping = ',', -- specifies what text object selects an arguments (ie a, and i, by default)
        -- if you're using `vim-textobj-parameter` you can also set this to `vim.g.vim_textobj_parameter_mapping`
    }},
}
