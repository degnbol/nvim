return {
    -- more languages supported but sometimes does nothing when trevj does the split
    -- there is also
    -- https://github.com/CKolkey/ts-node-action/
    -- and mini splitjoin which doesn't use treesitter
    -- https://www.reddit.com/r/neovim/comments/11mtm3m/minisplitjoin_split_and_join_arguments/
    -- I tried it out and it joins comments by uncommenting them. It has 
    -- ability to be customized with hooks, e.g. potentially to turn line 
    -- comments into comment regions, but I think we will never actually want 
    -- to do this.
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
                use_default_keymaps = false,
                langs = { lua =  { both = { no_format_with = nil }}},
            }
        end,
        init = function ()
            -- but always map toggle since the other plugins doesn't implement it
            -- alt use "<Plug>JoinToggle" if you want to make binding somewhere else
            vim.keymap.set("n", "<leader>J", "<Cmd>TSJToggle<CR>", { desc="Toggle split/join" })

            -- use alt plugins for specific filetypes
            -- it is also a possibility to define them for treesj
            local fts_split = {
                splitjoin = { tex=true, },
                trevj = { julia=true, lua=false, },
            }
            local fts_join  = {
                splitjoin = { tex=true, },
                matchup = { julia=true,},
            }
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
                        vim.keymap.set("n", key_split, "<Cmd>SplitjoinSplit<CR>", { buffer=true, desc="Split (splitjoin)" })
                    elseif fts_split['trevj'][vim.bo.filetype] then
                        vim.keymap.set("n", key_split, function ()
                            return require"trevj".format_at_cursor()
                        end, { buffer=true, desc="Split (trevj)" })
                    else
                        vim.keymap.set("n", key_split, "<Cmd>TSJSplit<CR>", { buffer=true, desc="Split (TSJ)" })
                    end
                    
                    if fts_join['splitjoin'][vim.bo.filetype] then
                        vim.keymap.set("n", key_join, "<Cmd>SplitjoinJoin<CR>", {buffer=true, desc="Join (splitjoin)" })
                    elseif fts_join['matchup'][vim.bo.filetype] then
                        -- Uses a% from https://github.com/andymass/vim-matchup which first jumps to container 
                        vim.keymap.set("n", key_join, "va%J", { buffer=true, desc="Join (matchup)" })
                    else
                        vim.keymap.set("n", key_join, "<Cmd>TSJJoin<CR>", { buffer=true, desc="Join (TSJ)" })
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
    -- delete continuation characters.
    -- Just like formatoptions+=j will remove leading comment chars, conjoin will remove trailing continuation markers.
    {
        'flwyd/vim-conjoin',
        -- depend on splitjoin so splitjoin is loaded first, which makes conjoin work with it.
        -- https://www.reddit.com/r/vim/comments/g71wyq/delete_continuation_characters_when_joining_lines/
        dependencies = { "AndrewRadev/splitjoin.vim" },
        init = function ()
            -- add custom language support
            vim.g.conjoin_filetypes = {
                asciidoc = {trailing='+$'},
            }
        end,
    },
}
