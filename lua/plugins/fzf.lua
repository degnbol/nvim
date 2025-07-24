return {
    {
        "folke/snacks.nvim",
        enabled = false,
        ---@type snacks.Config
        opts = {
            picker = {
            },
            explorer = {}
        },
        keys = {
            -- TODO
        },
    },
    -- We prefer fzf-lua over telescope for the following reasons:
    -- - There are claims that it is faster.
    -- - Telescope is not as actively developed/maintained.
    -- - Telescope (at least by default) needs double ESC to exit.
    --   This is to allow for normal mode editing of search but the idea of fzf
    --   searching is quick and messy search so we prefer single ESC exit.
    -- - The preview is worse (at least by default), i.e. the line highlight doesn't span the whole pane,
    --   and the relevant col is not indicated.
    -- - There is a bit more builtins from fzf lua.
    {
        "ibhagwan/fzf-lua",
        lazy = true,
        cmd = "FzfLua",
        -- optional for icon support
        dependencies = { "nvim-tree/nvim-web-devicons" },
        init = function()
            -- Simple picker with TAB (by default) for toggling preview, where fzf lua has preview by default.
            require 'mini.pick'.setup()
            -- contains the explorer picker that allows moving through folders.
            local MiniExtra = require 'mini.extra'
            MiniExtra.setup()

            local function nmap(keys, func, desc, leader)
                leader = leader or "<leader>"
                vim.keymap.set('n', leader .. keys, func, { desc = desc })
            end
            local function map(keys, func, desc, args, leader)
                -- use anon funcs so require only gets called when fzf is needed, i.e. lazy load.
                nmap(keys, function()
                    require 'fzf-lua'[func](args or {})
                end, desc, leader)
            end

            -- General starting point
            map("f<leader>", "builtin", "Builtin", {})
            map("ff", "resume", "Resume")

            -- Buffers and Files

            -- starts in folder of current buffer (like Oil)
            nmap("fp", function() MiniExtra.pickers.explorer { cwd = vim.fn.expand("%:h") } end, "Path explorer")
            -- differs by starting in CWD, which may be different from dir for current buffer
            map("fP", "files", "PWD files")
            map("fD", "files", "~/dotfiles/", { cwd = "~/dotfiles/" })
            map("fN", "files", "~/nvim/", { cwd = "~/nvim/" })
            map("fb", "buffers", "Open buffers")
            map("fo", "oldfiles", "Opened files history")
            -- map("ft", "tabs", "Tabs")
            map("fq", "quickfix", "Quickfix list")
            map("fQ", "quickfix_stack", "Quickfix stack")
            -- location list is window local quickfix list, see :h location-list
            map("fl", "loclist", "Location list")
            map("fL", "loclist_stack", "Location stack")
            map("f$", "args", "Argument list")
            map("f_", "grep_curbuf", "Grep current buffer lines")
            -- map("f_", "lines", "Lines") -- redudant with the grep_curbuf
            map("f-", "grep_project", "Grep project")
            map("fx", "tmux_buffers", "List tmux paste buffers")

            map("ts", "treesitter", "Symbols", { prompt = "Symbols❯ " })

            -- Regex pattern search and file content grep refinement

            map("fg", "grep", "Grep", { input_prompt = "Grep❯ " })
            -- we already have resume with <leader>ff so not so important:
            map("fG", "grep_last", "Grep last")
            map("fw", "grep_cword", "Grep cword")
            map("fW", "grep_cWORD", "Grep cWORD")
            vim.keymap.set('v', '<leader>fg', function() require "fzf-lua".grep_visual {} end, { desc = "Grep" })
            -- seems redundant:
            -- map("/?", "blines", "Current buffer lines")
            -- Won't map grep_quickfix and lgrep_quickfix
            -- since they didn't work in latex

            -- Tags

            map("]]", "btags", "Tags in buffer")
            map("f]", "tags", "Tags")
            map("]g", "tags_grep", "Tags regex")
            map("]w", "tags_grep_cword", "Tag grep cword")
            map("]W", "tags_grep_cWORD", "Tag grep cWORD")

            -- Git
            map("gf", "git_files", "Files")
            map("gc", "git_bcommits", "Commits (current buf)")
            map("gC", "git_commits", "Commits")
            map("gB", "git_branches", "Branches")
            map("gs", "git_status", "Status")

            -- Misc

            map("fh", "highlights", "highlight groups")
            map("fK", "helptags", "Help tags")
            map("fm", "marks", ":marks")
            map("fM", "manpages", "Manual pages")
            map("fj", "jumps", ":jumps")
            map("fc", "colorschemes", "Colorschemes")
            map("fC", "awesome_colorschemes", "Awesome colorschemes")
            map("f:", "commands", "Ex commands")
            map("f;", "command_history", "Ex command history")
            map("f/", "search_history", "Search history")
            map("fr", "registers", ":registers")
            map("fa", "autocmds", ":autocmd")
            map("fk", "keymaps", "key mappings")
            map(":t", "filetypes", "Set filetype")
            map(":]", "tagstack", ":tags")
            map(":c", "changes", ":changes")
            map(":p", "packadd", ":packadd")

            -- Editing. Prefix with x like builtin completion with ctrl-x ctrl-s etc.
            map("xs", "spell_suggest", "Spelling suggestions")
            map("xf", "complete_file", "Complete file under cursor (excl dirs)")
            map("xF", "complete_path", "Complete path under cursor (incl dirs)")
            map("x_", "complete_bline", "Complete line (current buffer only)")


            -- LSP in LSP.lua

            -- Rare. access from builtins instead:
            -- profiles
            -- menus
        end,
        config = function()
            local fzflua = require "fzf-lua"
            fzflua.setup {
                winopts = {
                    preview = {
                        -- Reduce lag from default 20 ms.
                        -- It is there for a purpose relating to fast scrolling.
                        -- https://github.com/ibhagwan/fzf-lua
                        -- Reducing it gives noticeable improvement in preview update speed.
                        delay = 10,
                    },
                },
                lsp = {
                    -- if there is a single LSP result only, as is common for goto def,
                    -- use it directly.
                    jump1 = true,
                    code_actions = {
                        async_or_timeout = 5000,
                        -- Requires git-delta for prettier code action preview
                        previewer        = "codeaction_native",
                    },
                },
                files = {
                    -- Same as default but adding build/ for exclusion.
                    find_opts = [[-type f \! -path '*/.git/*' -and \! -path '*/build/*']],
                    rg_opts   = [[--color=never --hidden --files -g '!.git' -g '!build']],
                    fd_opts   = [[--color=never --hidden --type f --type l --exclude .git --exclude build]],
                },
            }
        end,
    },
    -- recommended compiled fuzzy finder for telescope. Cannot be opt=true when needed by tzachar/cmp-fuzzy-path
    {
        "nvim-telescope/telescope-fzf-native.nvim",
        build = 'make',
        lazy = true, -- load as (fake) dependency of telescope
        config = function()
            -- load the native fzf as recommended
            -- require'telescope'.load_extension('fzf')
            -- pcall is protected call, i.e. doesn't make a big deal out of errors
            -- taken from https://github.com/nvim-lua/kickstart.nvim/blob/master/init.lua
            pcall(require('telescope').load_extension, 'fzf')
        end
    },
    {
        "nvim-telescope/telescope.nvim",
        cmd = "Telescope",
        dependencies = {
            "nvim-lua/plenary.nvim",
            "nvim-telescope/telescope-fzf-native.nvim",
        },
        opts = {
            defaults = {
                layout_strategy = "vertical",
                prompt_prefix = '',
                selection_caret = '',
                entry_prefix = '',
                multi_icon = '',
            },
            pickers = {
                colorscheme = {
                    -- window has to be big enough to show the preview window
                    enable_preview = true,
                }
            },
            extensions = {
                dash = {
                    file_type_keywords = {
                        python = { "python", "numpy", "scipy", "pandas" }
                    }
                }
            }
        }
    },
    {
        'sudormrfbin/cheatsheet.nvim',
        init = function()
            vim.keymap.set("n", "<leader>f?", "<Cmd>Cheatsheet<CR>", { desc = "Cheatsheet" })
        end,
        cmd = "Cheatsheet",
        dependencies = { 'nvim-telescope/telescope.nvim' }
    },
    -- TODO we don't actually use this yet
    -- search stackoverflow quicker
    { "lalitmee/browse.nvim", enabled = false, dependencies = { "nvim-telescope/telescope.nvim" } },
    {
        "nvim-telescope/telescope-bibtex.nvim",
        dependencies = { 'nvim-telescope/telescope.nvim' },
        ft = "tex",
        config = function()
            -- https://github.com/nvim-telescope/telescope-bibtex.nvim
            local telescope = require "telescope"
            telescope.load_extension("bibtex")
            -- the following is used to detect *.bib file used
            -- by looking for \bibliography and \addbibresource.
            telescope.setup { context = true, }
        end
    },
    -- consider papis as well.
    -- bibliography references, mostly relevant for citations in .tex documents.
    -- { "jghauser/papis.nvim",
    --     dependencies = { "kkharji/sqlite.lua", "nvim-lua/plenary.nvim", "MunifTanjim/nui.nvim", "nvim-treesitter/nvim-treesitter", "nvim-telescope/telescope.nvim", "hrsh7th/nvim-cmp" },
    --     rocks="lyaml", config = function() require("papis").setup() end,
    -- },
    {
        "debugloop/telescope-undo.nvim",
        dependencies = { -- note how they're inverted to above example
            {
                "nvim-telescope/telescope.nvim",
                dependencies = { "nvim-lua/plenary.nvim" },
            },
        },
        keys = {
            { -- lazy style key map
                "<leader>fu",
                "<cmd>Telescope undo<CR>",
                desc = "undo history",
            },
        },
        opts = {
            -- don't use `defaults = { }` here, do this in the main telescope spec
            extensions = {
                undo = {
                    -- telescope-undo.nvim config, see below
                },
                -- no other extensions here, they can have their own spec too
            },
        },
        config = function(_, opts)
            -- Calling telescope's setup from multiple specs does not hurt, it will happily merge the
            -- configs for us. We won't use data, as everything is in it's own namespace (telescope
            -- defaults, as well as each extension).
            require("telescope").setup(opts)
            require("telescope").load_extension("undo")

            vim.keymap.set("n", "<leader>fu", "<cmd>Telescope undo<CR>")
        end,
    },
}
