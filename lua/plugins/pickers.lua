local hi = require "utils/highlights"
local map = require "utils/keymap"

-- Toggle between pickers by changing this flag.
-- Options: "fzf-lua", "snacks"
-- local picker = "fzf-lua"
local picker = "snacks"

-- Shared keymap definitions. Each entry: { keys, fzf_func, snacks_func, desc, opts }
-- opts is picker-specific and passed through.
-- snacks_func can be a string (Snacks.picker method name) or nil (no snacks equivalent).
local keymaps = {
    -- General
    { "Fl",  "builtin",      "pickers",          "All pickers" },
    { "f<leader>", "resume", "resume",            "Resume" },

    -- Buffers and Files
    { "fP",  "files",         "files",            "PWD files" },
    { "fD",  "files",         "files",            "~/dotfiles/",       fzf = { cwd = "~/dotfiles/" }, snacks = { dirs = { "~/dotfiles/" } } },
    { "fN",  "files",         "files",            "~/nvim/",           fzf = { cwd = "~/nvim/" },     snacks = { dirs = { "~/nvim/" } } },
    { "fb",  "buffers",       "buffers",          "Open buffers" },
    { "fo",  "oldfiles",      "recent",           "Opened files history" },
    { "fq",  "quickfix",      "qflist",           "Quickfix list" },
    { "fQ",  "quickfix_stack", nil,               "Quickfix stack" },
    { "fl",  "loclist",       "loclist",          "Location list" },
    { "fL",  "loclist_stack", nil,                "Location stack" },
    { "f$",  "args",          nil,                "Argument list" },
    { "f_",  "grep_curbuf",   "lines",           "Grep current buffer lines" },
    { "f-",  "grep_project",  "grep",            "Grep project" },
    { "fx",  "tmux_buffers",  nil,                "List tmux paste buffers" },

    { "ts",  "treesitter",    "treesitter",       "Symbols" },

    -- Grep
    { "fg",  "grep",          "grep",             "Grep" },
    { "fG",  "grep_last",     nil,                "Grep last" },
    { "fw",  "grep_cword",    "grep_word",        "Grep cword" },
    { "fW",  "grep_cWORD",    "grep_word",        "Grep cWORD",
        snacks = { search = function() return vim.fn.expand("<cWORD>") end } },

    -- Tags
    { "]]",  "btags",         "tags",             "Tags in buffer" },
    { "f]",  "tags",          "tags",             "Tags" },
    { "]g",  "tags_grep",     nil,                "Tags regex" },
    { "]w",  "tags_grep_cword", nil,              "Tag grep cword" },
    { "]W",  "tags_grep_cWORD", nil,              "Tag grep cWORD" },

    -- Git
    { "gf",  "git_files",     "git_files",        "Files" },
    { "gc",  "git_bcommits",  "git_log_file",     "Commits (current buf)" },
    { "gC",  "git_commits",   "git_log",          "Commits" },
    { "gB",  "git_branches",  "git_branches",     "Branches" },
    { "gs",  "git_status",    "git_status",       "Status" },

    -- Misc
    { "fh",  "highlights",    "highlights",       "Highlight groups" },
    { "fK",  "helptags",      "help",             "Help tags" },
    { "fm",  "marks",         "marks",            ":marks" },
    { "fM",  "manpages",      "man",              "Manual pages" },
    { "fj",  "jumps",         "jumps",            ":jumps" },
    { "fc",  "colorschemes",  "colorschemes",     "Colorschemes" },
    { "fC",  "awesome_colorschemes", nil,          "Awesome colorschemes" },
    { "f:",  "commands",      "commands",          "Ex commands" },
    { "f;",  "command_history", "command_history", "Ex command history" },
    { "f/",  "search_history", "search_history",  "Search history" },
    { "fr",  "registers",    "registers",         ":registers" },
    { "fa",  "autocmds",     "autocmds",          ":autocmd" },
    { "fk",  "keymaps",      "keymaps",           "Key mappings" },
    { ":t",  "filetypes",    nil,                 "Set filetype" },
    { ":]",  "tagstack",     nil,                 ":tags" },
    { ":c",  "changes",      nil,                 ":changes" },
    { ":p",  "packadd",      nil,                 ":packadd" },

    -- Editing
    { "xs",  "spell_suggest", "spelling",         "Spelling suggestions" },
    { "xf",  "complete_file", nil,                "Complete file under cursor (excl dirs)" },
    { "xF",  "complete_path", nil,                "Complete path under cursor (incl dirs)" },
    { "x_",  "complete_bline", nil,               "Complete line (current buffer only)" },
}

--- Ensure the active picker plugin is loaded (packadd + lz.n before/after).
local function load_picker()
    if picker == "fzf-lua" then
        require("lz.n").trigger_load("fzf-lua")
    elseif picker == "snacks" then
        require("lz.n").trigger_load("snacks.nvim")
    end
end

--- Register all shared keymaps for the active picker.
local function register_keymaps()
    for _, km in ipairs(keymaps) do
        local keys, fzf_func, snacks_func, desc = km[1], km[2], km[3], km[4]

        if picker == "fzf-lua" then
            map.n("<leader>" .. keys, function()
                load_picker()
                require("fzf-lua")[fzf_func](km.fzf or {})
            end, desc)
        elseif picker == "snacks" and snacks_func then
            map.n("<leader>" .. keys, function()
                load_picker()
                require("snacks").picker[snacks_func](km.snacks or {})
            end, desc)
        end
    end

    -- Visual grep (not expressible in the shared table)
    if picker == "fzf-lua" then
        map.x("<leader>fg", function()
            load_picker()
            require("fzf-lua").grep_visual {}
        end, "Grep")
    elseif picker == "snacks" then
        map.x("<leader>fg", function()
            load_picker()
            require("snacks").picker.grep_word()
        end, "Grep")
    end
end

-- Register keymaps and mini.pick eagerly (before any picker loads).
-- The keymap callbacks require("fzf-lua")/require("snacks") lazily on invocation.
register_keymaps()

-- mini.pick path explorer (independent of picker choice).
-- Setup deferred to first invocation since mini.pick is an opt package.
local mini_pick_ready = false
map.n("<leader>fp", function()
    if not mini_pick_ready then
        vim.cmd.packadd("mini.nvim")
        require("mini.pick").setup {
            options = { content_from_bottom = true, use_cache = true },
            window = { config = { width = 100 } },
        }
        require("mini.extra").setup()
        mini_pick_ready = true
    end
    require("mini.extra").pickers.explorer { cwd = vim.fn.expand("%:h") }
end, "Path explorer")

return {
    {
        "snacks.nvim",
        ---@type snacks.Config
        after = function()
            require("snacks").setup {
                picker = {
                    enabled = true,
                    layout = { fullscreen = true },
                },
                explorer = { enabled = true },
                image = {
                    enabled = true,
                    -- Keep struts visible: the inline-math override below wraps
                    -- expressions in a transparent \strut, which snacks' default
                    -- `-trim` would crop straight back off.
                    convert = {
                        magick = {
                            math = { "-density", 192, "{src}[{page}]" }, -- snacks default appends "-trim"
                        },
                    },
                    math = {
                        latex = {
                            font_size = "normalsize", -- default Large
                        }
                    }
                },
            }

            -- Render inline `$...$` (and `\(...\)`) math flowing within the line,
            -- while `$$...$$` / `\[...\]` stay display blocks. snacks discards the
            -- delimiter and wraps every expression as display `\[...\]`, whose
            -- vertical glue makes even a single glyph ~3 cells tall — over
            -- placement.lua's `height<=2` inline-collapse threshold, so any inline
            -- expression followed by text renders below the line with an icon. We
            -- rewrite inline math to `\begin{math}<strut>...\end{math}` before
            -- snacks' transform: math mode drops the glue, and the strut floors
            -- every expression to one line-height box (uniform 2 cells) so it
            -- collapses to a single inline row without magnifying short glyphs.
            -- The strut is asymmetric (not \strut's 70/30): snacks fills the cell
            -- with no baseline awareness and maps the box bottom to the line
            -- bottom, so the height/depth split must mirror the editor font's
            -- baseline ratio or descender-less glyphs hover. Two knobs, both in
            -- \baselineskip units: their *ratio* (0.15 : 1.02 ~ 15 : 85) is the
            -- baseline split — 0.15 depth is the floor, any shallower and a
            -- subscript (k_{cat}) drops below the strut and grows the box so
            -- `$k$` and `$k_{cat}$` stop matching. Their *sum* is the size knob:
            -- the glyph ink is fixed, so a taller box shrinks the glyph once the
            -- collapse squeezes the box into one row. Scale both together to
            -- resize without disturbing the baseline.
            local strut = "\\rule[-0.18\\baselineskip]{0pt}{1.02\\baselineskip}"
            local doc = require("snacks").image.doc
            local latex = doc.transforms.latex
            doc.transforms.latex = function(img, ctx)
                if img.content and img.ext == "math.tex" then
                    local raw = vim.trim(img.content)
                    if raw:match("^%$[^$]") or raw:match("^\\%(") then
                        local inner = raw:gsub("^%$+`?", ""):gsub("`?%$+$", "")
                            :gsub("^\\%(", ""):gsub("\\%)$", "")
                        img.content = ("\\begin{math}%s%s\\end{math}"):format(strut, inner)
                    end
                end
                return latex(img, ctx)
            end

            -- Make collapsed inline math occupy exactly its source width.
            -- placement.lua:489 sizes the image by pixel aspect + 2 padding
            -- cells, but it conceals the source `$...$` and replaces it inline —
            -- so any mismatch reflows the line. Width = source display width keeps
            -- the glyph where the LaTeX was (left-aligned, no shift when conceal
            -- toggles under the cursor). At conceallevel 1 the concealed `$...$`
            -- leaves one residual space cell after the image (level 2 hides it
            -- fully), so drop a cell there to land back on the source width.
            local placement = require("snacks").image.placement
            local state = placement.state
            placement.state = function(self)
                local st = state(self)
                local r = self.opts.range
                if self.opts.type == "math" and st.loc.height == 1 and r and r[1] == r[3] then
                    local src = vim.api.nvim_buf_get_text(self.buf, r[1] - 1, r[2], r[3] - 1, r[4], {})[1]
                    local win = st.wins[1]
                    if src and win then
                        local residual = vim.wo[win].conceallevel == 1 and 1 or 0
                        st.loc.width = math.max(1, vim.api.nvim_strwidth(src) - residual)
                    end
                end
                return st
            end

            -- Persist overrides across colorscheme changes. Snacks has a
            -- managed ColorScheme autocmd that re-applies its defaults after
            -- :highlight clear — onColorScheme fires after it, so our force
            -- links win.
            hi.onColorScheme(function()
                hi.link("SnacksPicker", "Normal")
                hi.link("SnacksPickerPreview", "Normal")
                hi.link("SnacksPickerMatch", "IncSearch")
                hi.link("SnacksPickerFile", "@string.special.path")
                hi.link("SnacksPickerDir", "@string.special.path")
                hi.link("SnacksPickerLink", "NonText")
                hi.link("SnacksPickerBorder", "FloatBorder")
                hi.link("SnacksPickerInputBorder", "FloatBorder")
                hi.link("SnacksPickerPreviewBorder", "FloatBorder")
                hi.link("SnacksPickerListBorder", "FloatBorder")
                hi.link("SnacksPickerTitle", "Title")
                hi.link("SnacksPickerInputTitle", "Title")
                -- Don't color math at all, it's clear from the symbols it's 
                -- different from regular text and normal colors makes it blend 
                -- better into the flow of reading.
                hi.link("SnacksImageMath", "Normal")
            end)
        end,
    },
    {
        "fzf-lua",
        enabled = picker == "fzf-lua",
        lazy = true,
        cmd = "FzfLua",
        before = function()
            -- hi.def runs before packadd, so these win over fzf-lua's own
            -- default=true in setup_highlights(). fzf-lua has no ColorScheme
            -- autocmd, so onColorScheme re-applies after :highlight clear.
            hi.onColorScheme(function()
                hi.link("FzfLuaBorder", "FloatBorder")
                hi.link("FzfLuaTitle", "Title")
                hi.link("FzfLuaDirPart", "@string.special.path")
                hi.link("FzfLuaFilePart", "@string.special.path")
            end)
        end,
        after = function()
            require("fzf-lua").setup {
                -- Use FuzzyMatch colour (colour 9) for fzf match highlighting.
                -- `true` enables auto-generation from highlight groups (reads FzfLuaFzfMatch).
                -- hl/hl+ override the treesitter default of "-1:reverse".
                fzf_colors = {
                    true,
                    ["hl"]  = { "fg", "FuzzyMatch" },
                    ["hl+"] = { "fg", "FuzzyMatch" },
                },
                winopts = {
                    backdrop = 100,
                    fullscreen = true,
                    treesitter = {
                        fzf_colors = {
                            ["hl"]  = { "fg", "FuzzyMatch" },
                            ["hl+"] = { "fg", "FuzzyMatch" },
                        },
                    },
                    preview = {
                        delay = 10,
                        scrollbar = false,
                        horizontal = "right:50%",
                        vertical = "down:50%",
                    },
                },
                keymap = {
                    builtin = {
                        ["<C-f>"]   = "preview-page-down",
                        ["<C-b>"]   = "preview-page-up",
                        ["<C-S-f>"] = "preview-down",
                        ["<C-S-b>"] = "preview-up",
                    },
                },
                lsp = {
                    jump1 = true,
                    code_actions = {
                        async_or_timeout = 5000,
                        previewer        = "codeaction_native",
                    },
                },
                -- Match ripgrep config (~/.config/ripgrep/config) colours explicitly
                -- so headless subprocesses don't rely on RIPGREP_CONFIG_PATH.
                grep = {
                    rg_opts = "--column --line-number --no-heading --color=always --smart-case "
                        .. "--max-columns=4096 "
                        .. "--colors=line:fg:3 --colors=column:fg:3 --colors=match:fg:9 --colors=match:style:bold "
                        .. "--colors=path:none --colors=path:style:underline -e",
                },
                files = {
                    find_opts = [[-type f \! -path '*/.git/*' -and \! -path '*/build/*']],
                    rg_opts   = [[--color=never --hidden --files -g '!.git' -g '!build']],
                    fd_opts   = [[--color=never --hidden --type f --type l --exclude .git --exclude build]],
                },
                oldfiles = {
                    include_current_session = true,
                },
            }
        end,
    },
    {
        "cheatsheet.nvim",
        before = function()
            map.n("<leader>f?", "<Cmd>Cheatsheet<CR>", "Cheatsheet")
        end,
        cmd = "Cheatsheet",
    },
    {
        "snacks-bibtex.nvim",
        ft = { "tex", "typst" },
        after = function()
            require("snacks-bibtex").setup {
                mappings = {
                    ["<C-p>"] = false,
                },
            }
        end,
        keys = {
            {
                "<leader>fc",
                function() require("snacks-bibtex").bibtex() end,
                desc = "BibTeX citations",
            },
        },
    },
}
