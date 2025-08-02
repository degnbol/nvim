local hi = require "utils/highlights"
local map = require "utils/keymap"

return {
    -- Still needs some polishing.
    {
        "dmtrKovalenko/fff.nvim",
        enabled = false,
        build = "cargo build --release",
        -- Config and opts has no effect.
        opts = {
            width = 1.0,
            height = 1.0,
            prompt = '',
        },
        keys = {
            {
                "<leader>fF",
                function()
                    require("fff").toggle()
                end,
                desc = "FFF cwd files",
            },
            {
                "<leader>fg",
                function()
                    require "fff".find_in_git_root()
                end,
                desc = "FFF git files",
            },
        },
    },
    {
        "folke/snacks.nvim",
        priority = 1000,
        ---@type snacks.Config
        opts = {
            picker = { enabled = true, layout = { fullscreen = true } },
            explorer = { enabled = true },
            quickfile = { enabled = false }, -- doesn't seem to make a difference.
        },
        keys = {
            { "\\",         function() require "snacks".picker.smart() end,      desc = "Smart Find Files" },
            { "<leader>fh", function() require "snacks".picker.highlights() end, desc = "Highlights" },
        },
        init = function()
            -- Instead of default float.
            hi.def("SnacksPicker", "Normal")
            hi.def("SnacksPickerPreview", "Normal")
            -- Instead of to special.
            hi.def("SnacksPickerMatch", "IncSearch")
        end,
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
            require 'mini.pick'.setup {
                options = {
                    -- Whether to show content from bottom to top
                    content_from_bottom = true,
                    -- Whether to cache matches (more speed and memory on repeated prompts)
                    use_cache = true,
                },
                window = {
                    config = { width = 100 },
                }
            }
            -- contains the explorer picker that allows moving through folders.
            local MiniExtra = require 'mini.extra'
            MiniExtra.setup()

            local function map_fzf(keys, func, desc, args)
                map.n("<leader>" .. keys, function()
                    require 'fzf-lua'[func](args or {})
                end, desc)
            end

            -- General starting point
            map_fzf("f<leader>", "builtin", "Builtin")
            map_fzf("ff", "resume", "Resume")

            -- Buffers and Files

            -- starts in folder of current buffer (like Oil)
            map.n("<leader>fp", function() MiniExtra.pickers.explorer { cwd = vim.fn.expand("%:h") } end, "Path explorer")
            -- differs by starting in CWD, which may be different from dir for current buffer
            map_fzf("fP", "files", "PWD files")
            map_fzf("fD", "files", "~/dotfiles/", { cwd = "~/dotfiles/" })
            map_fzf("fN", "files", "~/nvim/", { cwd = "~/nvim/" })
            map_fzf("fb", "buffers", "Open buffers")
            map_fzf("fo", "oldfiles", "Opened files history")
            -- map("ft", "tabs", "Tabs")
            map_fzf("fq", "quickfix", "Quickfix list")
            map_fzf("fQ", "quickfix_stack", "Quickfix stack")
            -- location list is window local quickfix list, see :h location-list
            map_fzf("fl", "loclist", "Location list")
            map_fzf("fL", "loclist_stack", "Location stack")
            map_fzf("f$", "args", "Argument list")
            map_fzf("f_", "grep_curbuf", "Grep current buffer lines")
            -- map("f_", "lines", "Lines") -- redudant with the grep_curbuf
            map_fzf("f-", "grep_project", "Grep project")
            map_fzf("fx", "tmux_buffers", "List tmux paste buffers")

            map_fzf("ts", "treesitter", "Symbols", { prompt = "Symbols❯ " })

            -- Regex pattern search and file content grep refinement

            map_fzf("fg", "grep", "Grep", { input_prompt = "Grep❯ " })
            -- we already have resume with <leader>ff so not so important:
            map_fzf("fG", "grep_last", "Grep last")
            map_fzf("fw", "grep_cword", "Grep cword")
            map_fzf("fW", "grep_cWORD", "Grep cWORD")
            map.x('<leader>fg', function() require "fzf-lua".grep_visual {} end, "Grep")
            -- seems redundant:
            -- map("/?", "blines", "Current buffer lines")
            -- Won't map grep_quickfix and lgrep_quickfix
            -- since they didn't work in latex

            -- Tags

            map_fzf("]]", "btags", "Tags in buffer")
            map_fzf("f]", "tags", "Tags")
            map_fzf("]g", "tags_grep", "Tags regex")
            map_fzf("]w", "tags_grep_cword", "Tag grep cword")
            map_fzf("]W", "tags_grep_cWORD", "Tag grep cWORD")

            -- Git
            map_fzf("gf", "git_files", "Files")
            map_fzf("gc", "git_bcommits", "Commits (current buf)")
            map_fzf("gC", "git_commits", "Commits")
            map_fzf("gB", "git_branches", "Branches")
            map_fzf("gs", "git_status", "Status")

            -- Misc

            -- The preview from snacks is better.
            -- This one previews other hl groups that are basically unrelated.
            -- map_fzf("fh", "highlights", "highlight groups")
            map_fzf("fK", "helptags", "Help tags")
            map_fzf("fm", "marks", ":marks")
            map_fzf("fM", "manpages", "Manual pages")
            map_fzf("fj", "jumps", ":jumps")
            map_fzf("fc", "colorschemes", "Colorschemes")
            map_fzf("fC", "awesome_colorschemes", "Awesome colorschemes")
            map_fzf("f:", "commands", "Ex commands")
            map_fzf("f;", "command_history", "Ex command history")
            map_fzf("f/", "search_history", "Search history")
            map_fzf("fr", "registers", ":registers")
            map_fzf("fa", "autocmds", ":autocmd")
            map_fzf("fk", "keymaps", "key mappings")
            map_fzf(":t", "filetypes", "Set filetype")
            map_fzf(":]", "tagstack", ":tags")
            map_fzf(":c", "changes", ":changes")
            map_fzf(":p", "packadd", ":packadd")

            -- Editing. Prefix with x like builtin completion with ctrl-x ctrl-s etc.
            map_fzf("xs", "spell_suggest", "Spelling suggestions")
            map_fzf("xf", "complete_file", "Complete file under cursor (excl dirs)")
            map_fzf("xF", "complete_path", "Complete path under cursor (incl dirs)")
            map_fzf("x_", "complete_bline", "Complete line (current buffer only)")


            -- LSP in LSP.lua

            -- Rare. access from builtins instead:
            -- profiles
            -- menus
        end,
        opts = {
            winopts = {
                -- Don't dim other windows.
                backdrop = 100,
                -- Fullscreen is actually nice for max focus.
                fullscreen = true,
                preview = {
                    -- Reduce lag from default 20 ms.
                    -- It is there for a purpose relating to fast scrolling.
                    -- https://github.com/ibhagwan/fzf-lua
                    -- Reducing it gives noticeable improvement in preview update speed.
                    delay = 10,
                    -- Doesn't move with the preview anyways.
                    scrollbar = false,
                },
            },
            keymap = {
                builtin = {
                    ["<C-f>"]   = "preview-page-down",
                    ["<C-b>"]   = "preview-page-up",
                    -- Single line movement for finer control.
                    ["<C-S-f>"] = "preview-down",
                    ["<C-S-b>"] = "preview-up",
                }
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
        },
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
            map.n("<leader>f?", "<Cmd>Cheatsheet<CR>", "Cheatsheet")
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
}
