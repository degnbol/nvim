return {
    -- more languages supported but sometimes does nothing when trevj does the split
    {
        "Wansmer/treesj",
        dependencies = {
            'nvim-treesitter/nvim-treesitter',
            -- make sure is loaded for the init to work
            'AckslD/nvim-trevJ.lua',
            'AndrewRadev/splitjoin.vim',
            'andymass/vim-matchup',
        },
        -- opts = {
        --     -- set only for supported filetypes
        --     use_default_keymaps = false,
        -- },
        config = function ()
            local treesj = require"treesj"
            -- local default_preset = require"treesj.langs.default_preset"
            -- default_preset['both']['no_format_with'] = {}
            -- require"treesj.langs".default_preset = default_preset
            treesj.setup {
                use_default_keymaps = false, langs = { lua =  { both = { no_format_with = nil }}},
            }
        end,
        init = function ()
            -- but always map toggle since the other plugins doesn't implement it
            -- alt use "<Plug>JoinToggle" if you want to make binding somewhere else
            vim.keymap.set("n", "<leader>tj", "<Cmd>TSJToggle<CR>", { desc="Split/join" })

            -- use alt plugins for specific filetypes
            -- it is also a possibility to define them for treesj
            local fts_split = { splitjoin = { tex=true, }, trevj = { julia=true, lua=true, } }
            local fts_join  = { splitjoin = { tex=true, }, matchup = { julia=true,}, }
            local grp = vim.api.nvim_create_augroup("splitjoin", {clear=true})

            -- call nmap buffer thru vimscript since nowait isn't implemented in lua API
            -- keep here for a momemt too see if there is any issues without 
            -- nowait, which would be better to avoid since we can then use the 
            -- lua mapping with desc field.
            function nowait(lhs, rhs)
                vim.cmd('nnoremap <buffer><nowait> ' .. lhs .. ' ' .. rhs)
            end
            function nowaitre(lhs, rhs)
                vim.cmd('nmap <buffer><nowait> ' .. lhs .. ' ' .. rhs)
            end

            -- alt use "<Plug>Split" if you wish to choose in another file
            local key_join = "<leader>j"
            local key_split = "<leader>s"
            
            vim.api.nvim_create_autocmd("Filetype", {
                pattern = "*", group = grp,
                callback = function ()
                    if fts_split['splitjoin'][vim.bo.filetype] then
                        -- nowait(key_join, "<Cmd>SplitjoinSplit<CR>")
                        vim.keymap.set("n", key_join, "<Cmd>SplitjoinSplit<CR>", { desc="Join" })
                    elseif fts_split['trevj'][vim.bo.filetype] then
                        -- nowait(key_split, ":lua require'trevj'.format_at_cursor()<CR>")
                        vim.keymap.set("n", key_split, function ()
                            return require"trevj".format_at_cursor()
                        end, { desc="Split" })
                    else
                        -- nowait(key_split, "<Cmd>TSJSplit<CR>")
                        vim.keymap.set("n", key_split, "<Cmd>TSJSplit<CR>", { desc="Split" })
                    end
                    
                    if fts_join['splitjoin'][vim.bo.filetype] then
                        -- nowait(key_join, "<Cmd>SplitjoinJoin<CR>")
                        vim.keymap.set("n", key_join, "<Cmd>SplitjoinJoin<CR>", { desc="Join" })
                    elseif fts_join['matchup'][vim.bo.filetype] then
                        -- Uses a% from https://github.com/andymass/vim-matchup which first jumps to container 
                        -- nowaitre(key_join, "va%J")
                        vim.keymap.set("n", key_join, "va%J", { desc="Join" })
                    else
                        -- nowait(key_join, "<Cmd>TSJJoin<CR>")
                        vim.keymap.set("n", key_join, "<Cmd>TSJJoin<CR>", { desc="Join" })
                    end
                end
            })
        end
    },
    {"AckslD/nvim-trevJ.lua",
    lazy = true,
    config=function()
        local make_default_opts = function()
            return {
                -- was using "," but it gets added to a final comment as well
                final_separator = "",
                final_end_line = true,
                skip = { },
            }
        end

        local make_no_final_sep_opts = function()
            return { final_separator = false, final_end_line = true }
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
    -- unmaintained but using it (my fork fixing linewise visual) since it's the only one with visual and motion
    {
        dir = "$XDG_CONFIG_HOME/nvim/nvim-revJ.lua",
        dev = true,
        dependencies={'kana/vim-textobj-user', 'sgur/vim-textobj-parameter'},
        opts={
            brackets = {first = '([{<', last = ')]}>'}, -- brackets to consider surrounding arguments
            new_line_before_last_bracket = true, -- add new line between last argument and last bracket (only if no last seperator)
            add_seperator_for_last_parameter = true, -- if a seperator should be added if not present after last parameter
            enable_default_keymaps = false,
            keymaps = {
                operator = 'gS', -- for operator (+motion)
                line = 'gSS', -- for formatting current line
                -- only works for visual char mode not line, fix it in init keymap
                visual = '<leader>s', -- for formatting visual selection
            },
            parameter_mapping = ',', -- specifies what text object selects an arguments (ie a, and i, by default)
            -- if you're using `vim-textobj-parameter` you can also set this to `vim.g.vim_textobj_parameter_mapping`
        },
    },
}
