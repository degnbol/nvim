#!/usr/bin/env lua
return {
    {
        'echasnovski/mini.nvim',
        version = false,
        priority = 100, -- higher than default 50 to allow mini parts loaded in other files.
        config = function ()

            require('mini.bracketed').setup {
                -- ]i to go to more indented region.
                indent = { suffix = 'i', options = { change_type="more" } },
                -- compliment with using [[, ]], [], ][ to jump to less indented region and to next region
            }

            -- mini_indentscope = require('mini.indentscope')
            -- mini_indentscope.setup {
            --     draw = { animation = mini_indentscope.gen_animation.none() },
            --     symbol = "▏",
            -- }


            require'mini.surround'.setup {
                mappings = {
                    add = 'ys', -- Add surrounding in Normal. Not Visual modes, see below.
                    delete = 'ds', -- Delete surrounding
                    find = '', -- Find surrounding (to the right)
                    find_left = '', -- Find surrounding (to the left)
                    highlight = '', -- Highlight surrounding
                    replace = 'cs', -- Replace surrounding
                    update_n_lines = '', -- Update `n_lines`
                },
                custom_surroundings = {
                    -- [d]elete, e.g. only as cstd to delete surround """
                    ['d'] = { output = { left = '', right = '' } },
                    -- normally used for tags such as <div> </div> which is cool but I don't use them so t for triple
                    ['t'] = {
                        -- figured out these patterns from the [[ ]] example in :h MiniSurround.config
                        input = { '%"%"%"().-()%"%"%"' },
                        output = { left = '"""', right = '"""' }
                    },
                    -- Taken from `:h MiniSurround.config` for use in lua.
                    -- I never really use the open versions of brackets that adds space so we override '['.
                    ['['] = {
                        input = { '%[%[().-()%]%]' },
                        output = { left = '[[', right = ']]' },
                    },
                    -- for latex ``...'', ~ is chosen since it is shift+`
                    ['~'] = {
                        input = { "``().-()''" },
                        output = { left = '``', right = "''" },
                    },
                },
                -- Number of lines within which surrounding is searched
                n_lines = 100,
            }

            -- Remap adding surrounding to Visual mode selection
            vim.keymap.del('x', 'ys') -- del to avoid y(ank) waiting
            vim.keymap.set('x', 'S', [[:<C-u>lua MiniSurround.add('visual')<CR>]], { noremap = true, desc="Surround" })
            -- Make special mapping for "add surrounding for line"
            vim.keymap.set('n', 'yss', 'ys_', { remap = true, desc="Surround line" })

            -- highlight hex colors and todos, notes etc
            -- consider this or https://github.com/folke/paint.nvim if you want to highlight custom things.
            local hipatterns = require'mini.hipatterns'
            hipatterns.setup {
                highlighters = {
                    -- Highlight standalone 'FIXME', 'HACK', 'TODO', 'NOTE'
                    fixme = { pattern = '%f[%w]()FIXME()%f[%W]', group = 'MiniHipatternsFixme' },
                    hack  = { pattern = '%f[%w]()HACK()%f[%W]',  group = 'MiniHipatternsHack'  },
                    todo  = { pattern = '%f[%w]()TODO()%f[%W]',  group = 'MiniHipatternsTodo'  },
                    note  = { pattern = '%f[%w]()NOTE()%f[%W]',  group = 'MiniHipatternsNote'  },

                    -- Highlight hex color strings (`#rrggbb`) using that color
                    hex_color = hipatterns.gen_highlighter.hex_color(),
                }
            }


            local clue = require'mini.clue'

            clue.setup({
                triggers = {
                    { mode = 'n', keys = '<Leader>' },
                    { mode = 'x', keys = '<Leader>' },
                    { mode = 'n', keys = 'g' },
                    { mode = 'x', keys = 'g' },
                    { mode = 'n', keys = 'y' },
                    { mode = 'x', keys = 'y' },
                    { mode = 'n', keys = 'z' },
                    { mode = 'x', keys = 'z' },
                    { mode = 'n', keys = ']' },
                    { mode = 'n', keys = '[' },
                    { mode = 'n', keys = '=' },
                    { mode = 'n', keys = '>' },
                    { mode = 'n', keys = '<' },
                    -- Built-in completion
                    { mode = 'i', keys = '<C-x>' },
                    -- Marks
                    { mode = 'n', keys = "'" },
                    { mode = 'n', keys = '`' },
                    { mode = 'x', keys = "'" },
                    { mode = 'x', keys = '`' },
                    -- Registers
                    { mode = 'n', keys = '"' },
                    { mode = 'x', keys = '"' },
                    { mode = 'i', keys = '<C-r>' },
                    { mode = 'c', keys = '<C-r>' },
                    -- Window
                    { mode = 'n', keys = '<C-w>' },
                    -- vimtex default map
                    { mode = 'n', keys = 'ts' },
                    -- list text objects
                    { mode = 'x', keys = 'i' },
                    { mode = 'x', keys = 'a' },
                },

                clues = {
                    -- use e.g. postkeys='<C-w>' to make a submode. 
                    { mode = 'n', keys = '<leader>a', desc = "Argument" },
                    { mode = 'n', keys = '<leader>b', desc = "Buffer/Tab" },
                    { mode = 'n', keys = '<leader>c', desc = "Compile|Code|Color" },
                    { mode = 'n', keys = '<leader>d', desc = "Diagnostic|Definition" },
                    { mode = 'n', keys = '<leader>f', desc = "Find|File" },
                    { mode = 'n', keys = '<leader>g', desc = "Git" },
                    { mode = 'n', keys = '<leader>l', desc = "LSP|Lang|Completion" },
                    { mode = 'n', keys = '<leader>m', desc = "Multicursor" },
                    -- UnconditionalPaste
                    { mode = 'n', keys = '<leader>p', desc = "Paste after" },
                    { mode = 'n', keys = '<leader>P', desc = "Paste before" },

                    { mode = 'n', keys = '<leader>r', desc = "Re|REPL" },
                    { mode = 'n', keys = '<leader>t', desc = "Tree-sitter" },
                    { mode = 'n', keys = '<leader>w', desc = "Workspace" },
                    { mode = 'n', keys = '<leader>q', desc = "Quickfix" },
                    { mode = 'n', keys = '<leader>Q', desc = "Locationlist" },
                    { mode = 'n', keys = ']s', desc = "Spell" },
                    { mode = 'n', keys = '[s', desc = "Spell" },
                    -- vim unimpaired
                    { mode = 'n', keys = 'yob', desc = "background" },
                    { mode = 'n', keys = 'yoh', desc = "hlsearch" },
                    { mode = 'n', keys = 'yoi', desc = "ignorecase" },
                    { mode = 'n', keys = 'yol', desc = "list" },
                    { mode = 'n', keys = 'yon', desc = "number" },
                    { mode = 'n', keys = 'yor', desc = "relativenumber" },
                    { mode = 'n', keys = 'yos', desc = "spell" },
                    { mode = 'n', keys = 'yow', desc = "wrap" },
                    { mode = 'n', keys = 'yo-', desc = "cursorline" },
                    { mode = 'n', keys = 'yo_', desc = "cursorline" },
                    { mode = 'n', keys = 'yox', desc = "cursorcolumn" },
                    { mode = 'n', keys = '=sh', desc = "hlsearch" },
                    { mode = 'n', keys = '=si', desc = "ignorecase" },
                    { mode = 'n', keys = '=sl', desc = "list" },
                    { mode = 'n', keys = '=sn', desc = "number" },
                    { mode = 'n', keys = '=sr', desc = "relativenumber" },
                    -- replaced by Substitute+reindent. 
                    -- Other motion could also be candidates for being replaced.
                    -- { mode = 'n', keys = '=ss', desc = "spell" },
                    { mode = 'n', keys = '=sw', desc = "wrap" },
                    { mode = 'n', keys = '<sl', desc = "list" },
                    { mode = 'n', keys = '>sl', desc = "nolist" },
                    { mode = 'n', keys = '<sn', desc = "number" },
                    { mode = 'n', keys = '>sn', desc = "nonumber" },
                    { mode = 'n', keys = '<sr', desc = "relativenumber" },
                    { mode = 'n', keys = '>sr', desc = "norelativenumber" },
                    { mode = 'n', keys = '<ss', desc = "spell" },
                    { mode = 'n', keys = '>ss', desc = "nospell" },
                    { mode = 'n', keys = '<sw', desc = "wrap" },
                    { mode = 'n', keys = '>sw', desc = "nowrap" },

                    -- Enhance this by adding descriptions for <Leader> mapping groups
                    clue.gen_clues.builtin_completion(),
                    clue.gen_clues.g(),
                    clue.gen_clues.marks(),
                    clue.gen_clues.registers(),
                    clue.gen_clues.windows {
                        -- making submodes of the window stuff means we can increase hight with <C-w>++ instead of <C-w>+<C-w>+
                        submode_move = true,
                        submode_navigate = true,
                        submode_resize = true,
                    },
                    clue.gen_clues.z(),
                },

                window = {
                    config = {
                        width = "auto",
                    },
                },
            })

            vim.schedule(function ()
                require "utils/keymap"
                set_keymap_desc('n', 'gc', "Comment")
                set_keymap_desc('n', 'gcc', "Line")
                set_keymap_desc('n', 'g/', "Last search")
                set_keymap_desc('n', '>>', "Indent line")
                set_keymap_desc('n', '<<', "Unindent line")
            end)

            require('mini.notify').setup()

            local MiniFiles = require 'mini.files'
            local minifiles_toggle = function(...)
                if not MiniFiles.close() then MiniFiles.open(...) end
            end
            vim.keymap.set('n', '<leader>e', minifiles_toggle, { desc="Toggle mini.files" })
            -- :h MiniFiles-examples
            local show_dotfiles = false
            local filter_show = function(fs_entry) return true end
            local filter_hide = function(fs_entry)
                return not vim.startswith(fs_entry.name, '.')
            end
            local toggle_dotfiles = function()
                show_dotfiles = not show_dotfiles
                local new_filter = show_dotfiles and filter_show or filter_hide
                MiniFiles.refresh { content = { filter = new_filter } }
            end

            MiniFiles.setup {
                content = {
                    filter = filter_hide
                },
                windows = {
                    preview = true,
                    -- Width of focused window
                    width_focus = 30,
                    -- width_nofocus = 15,
                    width_preview = 60,
                }
            }

            vim.api.nvim_create_autocmd('User', {
                pattern = 'MiniFilesBufferCreate',
                callback = function(args)
                    local buf_id = args.data.buf_id
                    -- overide - so an oil buffer isn't opened within the mini files buffer
                    vim.keymap.set('n', '-', 'h', { buffer = buf_id, remap=true })
                    vim.keymap.set('n', '<leader>bd', MiniFiles.close, { buffer = buf_id, })
                    -- There is no option for force closing without prompt
                    vim.keymap.set('n', '<leader>bD', MiniFiles.close, { buffer = buf_id, })
                    -- use l to open file without closing explorer
                    vim.keymap.set('n', '<Enter>', function () MiniFiles.go_in({close_on_file=true}) end, {
                        buffer = buf_id, desc="Open file and close explorer",
                    })
                    vim.keymap.set('n', 'g.', toggle_dotfiles, { buffer = buf_id })
                    -- NOTE: couldn't get a function working to toggle preview.
                end,
            })

        end
    },
}
