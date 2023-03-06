return {
    -- if it fails (gj), try revj (<leader>j[j] or motion INSIDE brackets)
    {"AckslD/nvim-trevJ.lua", config=function()
        local make_default_opts = function()
            return {
                final_separator = ",",
                final_end_line = true,
                skip = {},
            }
        end

        local make_no_final_sep_opts = function()
            return {
                final_separator = false,
                final_end_line = true,
            }
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
    {"AckslD/nvim-revJ.lua", dependencies={'kana/vim-textobj-user', 'sgur/vim-textobj-parameter'}, opts={
        brackets = {first = '([{<', last = ')]}>'}, -- brackets to consider surrounding arguments
        new_line_before_last_bracket = true, -- add new line between last argument and last bracket (only if no last seperator)
        add_seperator_for_last_parameter = true, -- if a seperator should be added if not present after last parameter
        enable_default_keymaps = false, -- enables default keymaps without having to set them below
        keymaps = {
            operator = '<Leader>j', -- for operator (+motion)
            line = '<Leader>jj', -- for formatting current line
            visual = '<Leader>j', -- for formatting visual selection
        },
        parameter_mapping = ',', -- specifies what text object selects an arguments (ie a, and i, by default)
        -- if you're using `vim-textobj-parameter` you can also set this to `vim.g.vim_textobj_parameter_mapping`
    }},
}
