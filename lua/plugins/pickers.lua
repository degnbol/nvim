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
        after = function()
            ---@type snacks.Config
            require("snacks").setup {
                -- https://github.com/folke/snacks.nvim/blob/main/docs/input.md
                input = {
                    enabled = true,
                    icon = "",
                },
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
                            -- Snippets insert real unicode (Δ, ∑, ∂, →, ≤, …;
                            -- luasnippets/tex/unicode.lua), so compile with
                            -- unicode-math under tectonic's XeTeX (reads
                            -- codepoints natively). unicode-math supersedes and
                            -- clashes with amssymb/amsfonts, so drop both from
                            -- the snacks default. The template sorts packages
                            -- alphabetically (doc.lua), so unicode-math lands
                            -- after amsmath/mathtools — the required load order.
                            -- Latin Modern Math (the default) covers the glyphs;
                            -- no \setmathfont needed.
                            packages = { "amsmath", "amscd", "mathtools", "unicode-math" },
                        }
                    }
                },
            }

            -- Workaround for folke/snacks.nvim#2539 (open; fix PR #2871 unmerged):
            -- the picker `select` source (backing `vim.ui.select`) sizes its list
            -- box to `math.min(#items, vim.o.lines * 0.8 - 10)`, which is fractional
            -- when `lines` isn't a multiple of 5 and #items exceeds it. snacks' win.lua
            -- passes that straight to nvim_win_set_config, which rejects non-integral
            -- dims (E5108). Hits large `vim.ui.select` menus (e.g. agentic's mode/
            -- config-option selector) intermittently, depending on window height.
            -- The buggy size is injected via the picker's own `layout.config`/
            -- `on_update_pre`, both of which the picker hard-overrides, so it can't be
            -- reached from config. Floor at the single chokepoint instead: `win:dim()`,
            -- whose return `win_opts()` hands to nvim_win_set_config. Remove once fixed
            -- upstream.
            local Win = require("snacks.win")
            if not Win._dim_floored then
                Win._dim_floored = true
                local dim = Win.dim
                function Win.dim(self, parent)
                    local d = dim(self, parent)
                    d.height, d.width = math.floor(d.height), math.floor(d.width)
                    d.row, d.col = math.floor(d.row), math.floor(d.col)
                    return d
                end
            end

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

            -- Prefetch margin: render math within one screenful above/below the
            -- viewport so it's ready on arrival instead of converting only once
            -- scrolled on-screen (conversions are content-hash cached, so the
            -- cost is paid once, ahead of time). find_visible (discovery) and
            -- inline.visible (keep-alive/reconcile) MUST use the same range — a
            -- prefetched off-screen placement that find_visible converts but
            -- visible() doesn't match gets recreated as a duplicate every tick.
            local function prefetch_range(info)
                local margin = info.botline - info.topline
                return math.max(info.topline - 1 - margin, 1), info.botline + margin
            end

            -- The line range the cursor "occupies" for conceal purposes: the
            -- visual selection (or just the cursor line), or nil when
            -- concealcursor keeps this mode concealed (image stays shown).
            -- Mirrors snacks' own conceal() check so the freeze below and the
            -- hide agree on which line is "being edited". Current window only.
            local function cursor_lines()
                local mode = vim.fn.mode():sub(1, 1):lower()
                if vim.wo.concealcursor:find(mode, 1, true) then
                    return nil
                end
                local a, b = vim.fn.line("v"), vim.api.nvim_win_get_cursor(0)[1]
                return math.min(a, b), math.max(a, b)
            end

            -- Whether a match/placement is inline math overlapping [cf, ct].
            -- range is {start_row, start_col, end_row, end_col}, 1-indexed rows.
            local function math_on_cursor(kind, range, cf, ct)
                return cf and kind == "math" and range and range[1] <= ct and range[3] >= cf
            end

            -- Reimplement find_visible (range is hardcoded inside it, so a
            -- wrapper can't inject the margin) — near-copy of doc.find_visible
            -- with three drops from the found set:
            --   * widen the query range by the prefetch margin (above).
            --   * conceallevel 0 shows the literal `$...$`, so drop math matches
            --     in any window at cl=0; reconcile then :close()s the placement,
            --     revealing the source. Non-math doc images pass through.
            --   * drop math on the cursor line (current buffer) so the
            --     actively-edited (often syntactically incomplete) expression
            --     isn't reconverted on every keystroke. Paired with the same drop
            --     in inline.visible below, the existing placement is left
            --     untouched — frozen alive (hidden by conceal), not reconverted,
            --     not closed — until the cursor leaves and conceal() re-renders.
            doc.find_visible = function(buf, cb)
                local cf, ct
                if vim.api.nvim_get_current_buf() == buf then
                    cf, ct = cursor_lines()
                end
                local ret = {}
                local wins = vim.fn.win_findbuf(buf)
                local count = #wins
                for _, win in ipairs(wins) do
                    local info = vim.fn.getwininfo(win)[1]
                    local cl0 = vim.wo[win].conceallevel == 0
                    local from, to = prefetch_range(info)
                    -- doc.find's iter_matches end_row is end-inclusive (covers
                    -- [from..to+1]) — the same off-by-one as stock get. With get
                    -- now fixed to [from..to], pass to-1 here so discovery and
                    -- keep-alive (visible→get) stay aligned on [from..to];
                    -- otherwise bottom-prefetch math is found but not kept alive
                    -- and gets recreated every tick. Cost: one line less bottom
                    -- prefetch (negligible slack margin).
                    doc.find(buf, function(matches)
                        for _, i in ipairs(matches) do
                            if not ((cl0 and i.type == "math") or math_on_cursor(i.type, i.range, cf, ct)) then
                                ret[i.id] = i
                            end
                        end
                        count = count - 1
                        if count == 0 and cb then
                            cb(vim.tbl_values(ret))
                        end
                    end, { from = from, to = to })
                end
            end

            -- Match find_visible's prefetch margin (the coupling: a prefetched
            -- off-screen placement that visible() doesn't match is recreated as a
            -- duplicate every tick) and its cursor-line drop (so the frozen
            -- placement is neither matched nor closed by reconcile). Reimplemented
            -- (not wrapped) because the range is hardcoded inside the method.
            local inline = require("snacks").image.inline

            -- Reimplement M:get with the correct end bound. Stock get
            -- (inline.lua) passes a 1-indexed-inclusive `to` straight into
            -- nvim_buf_get_extmarks' 0-indexed end row {to, -1}, so it covers
            -- 1-indexed lines [from..to+1] — one line too many. conceal's
            -- self:get(cf, ct) then hides the cursor line *plus* the line below
            -- it (that line shows its source = the cursorline bug). The 0-indexed
            -- inclusive end for 1-indexed `to` is to-1, so query {to-1, -1}.
            -- conceal/visible call this and inherit the fix unchanged.
            local ns = require("snacks").image.placement.ns
            inline.get = function(self, from, to)
                local ret = {}
                local marks = vim.api.nvim_buf_get_extmarks(
                    self.buf, ns, { from - 1, 0 }, { to - 1, -1 }, { overlap = true, hl_name = false })
                for _, m in ipairs(marks) do
                    local p = self.idx[m[1]]
                    if p and not self.imgs[p.id] then
                        self.idx[m[1]] = nil
                        p = nil
                    end
                    if p then
                        ret[p.id] = p
                    end
                end
                return ret
            end

            inline.visible = function(self)
                local cf, ct
                if vim.api.nvim_get_current_buf() == self.buf then
                    cf, ct = cursor_lines()
                end
                local ret = {}
                for _, win in ipairs(vim.fn.win_findbuf(self.buf)) do
                    local info = vim.fn.getwininfo(win)[1]
                    local from, to = prefetch_range(info)
                    for k, v in pairs(self:get(from, to)) do
                        if not math_on_cursor(v.opts.type, v.opts.range, cf, ct) then
                            ret[k] = v
                        end
                    end
                end
                return ret
            end

            -- Drive cursor-line conceal as one synchronous state. Stock conceal()
            -- (inline.lua) img:show()s ALL images, THEN hides the cursor-line
            -- ones — the show-then-hide on that line is the flash on entry. Build
            -- the hide-set first and only show() what isn't hidden, so the cursor
            -- line is never momentarily re-shown over its source. The cursor-line
            -- placement stays alive (frozen by the find/visible drops above) so
            -- leaving re-shows it instantly. Editing it leaves it stale; nothing
            -- else fires update() on cursor-exit, so re-render here — debounced,
            -- only on line change so horizontal moves within an expression don't
            -- thrash.
            -- A placement's live extmark row (1-indexed) — what nvim tracked the
            -- image to under edits, vs the stored opts.range which only refreshes
            -- on reconcile.
            local function live_row(p)
                local eid = p.eids and p.eids[1]
                local m = eid and vim.api.nvim_buf_get_extmark_by_id(p.buf, ns, eid, {})
                return m and m[1] and m[1] + 1 or nil -- nil if the mark is gone
            end

            -- placement._render re-pins every extmark to opts.range[1]-1
            -- (placement.lua:386), discarding the row the extmark tracked to.
            -- After a structural edit (line :move via ]e/[e, dd, …) a frozen
            -- cursor-line placement's opts.range is stale, so the next render
            -- snaps its image back to the old row — divorced from its text — and
            -- the reconcile then can't match it and spawns duplicates. Sync
            -- opts.range to the live extmark span before rendering, but only when
            -- the row actually moved (in-place edits keep the frozen range so a
            -- half-typed expression isn't re-measured mid-edit).
            local function sync_range(p)
                local eid = p.eids and p.eids[1]
                local m = eid and vim.api.nvim_buf_get_extmark_by_id(p.buf, ns, eid, { details = true })
                local r = p.opts.range
                if m and m[1] and r and m[1] + 1 ~= r[1] then
                    p.opts.range = { m[1] + 1, m[2], (m[3].end_row or m[1]) + 1, m[3].end_col or m[2] }
                end
            end

            inline.conceal = function(self)
                local cf, ct = cursor_lines()
                local hide = cf and self:get(cf, ct) or {}
                for id, img in pairs(self.imgs) do
                    sync_range(img)
                    if hide[id] and img.opts.conceal then
                        -- img:hide() debounces its re-render 10ms (placement.lua),
                        -- so on entry the image lingers ~10ms over the source vim
                        -- has already revealed — the flash. Force the real (un-
                        -- debounced) update so the image vanishes the same redraw.
                        img:hide()
                        require("snacks").image.placement.update(img)
                    else
                        img:show()
                    end
                end
                local key = cf and (cf .. ":" .. ct) or ""
                if key ~= self._cl_key then
                    self._cl_key = key
                    self._refresh = self._refresh or require("snacks").util.debounce(function()
                        self:update()
                    end, { ms = 100 })
                    self._refresh()
                end
            end

            -- Reconcile binds each fresh match to a live placement by content
            -- hash (img.src = the cache PNG path, so the LaTeX→PNG compile is
            -- already deduped + disk-cached — identical exprs never recompile).
            -- But identical adjacent `$…$` share one src, so stock's
            -- first-src-wins (inline.lua:100) picks arbitrarily among them and
            -- can hand a match the placement sitting on a *different* line —
            -- swapping their stored opts.range and recreating both extmarks
            -- (churn + flicker on edit). Extmarks track text edits live, so the
            -- placement whose live extmark row already equals the match's row IS
            -- the right one: prefer it, fall back to first-src-wins only when no
            -- live extmark lands on the row (mark invalidated by an edit across
            -- the expression, or a brand-new match).
            -- ponytail: row-only match; two identical `$…$` on one line still
            -- bind by pairs() order (no row relocation to flicker, so harmless).
            local Snacks = require("snacks")
            inline.update = function(self)
                -- conceal config is a bool or a (lang, type) -> bool predicate.
                local conceal_cfg = Snacks.image.config.doc.conceal
                doc.find_visible(self.buf, function(imgs)
                    local visible = self:visible()
                    for _, i in ipairs(imgs) do
                        local img
                        local tr = i.range and i.range[1]
                        -- pass 1: position-stable (same src AND same live row);
                        -- pass 2: first-src-wins fallback.
                        for pass = 1, 2 do
                            if img then break end
                            for v, o in pairs(visible) do
                                if o.img.src == i.src and (pass == 2 or live_row(o) == tr) then
                                    img, visible[v] = o, nil
                                    break
                                end
                            end
                        end
                        if not img then
                            img = Snacks.image.placement.new(
                                self.buf,
                                i.src,
                                Snacks.config.merge({}, Snacks.image.config.doc, {
                                    pos = i.pos,
                                    range = i.range,
                                    inline = true,
                                    conceal = vim.b[self.buf].snacks_image_conceal
                                        or (type(conceal_cfg) == "function" and conceal_cfg(i.lang, i.type) or conceal_cfg),
                                    type = i.type,
                                    on_update = function(p)
                                        for _, eid in ipairs(p.eids) do
                                            self.idx[eid] = p
                                        end
                                    end,
                                })
                            )
                            for _, eid in ipairs(img.eids) do
                                self.idx[eid] = img
                            end
                            self.imgs[img.id] = img
                        else
                            img.opts.pos = i.pos
                            img.opts.range = i.range
                            img:update()
                        end
                    end
                    for _, img in pairs(visible) do
                        img:close()
                        self.imgs[img.id] = nil
                    end
                end)
            end

            -- Size collapsed inline math by how hard the source is concealed.
            -- The image stretch-fills loc.width × 1 cells, so the right width
            -- depends on whether the source footprint is still on screen:
            --   1: `$...$` collapses to one residual space cell — match source
            --      width minus that cell so surrounding text doesn't move.
            --   2+: `$...$` fully hidden — size to the glyph's true width (pixel
            --      aspect at one cell tall, rounded to whole cells) so every
            --      glyph renders at a consistent size with no surrounding
            --      whitespace. snacks' own ceil(w/h)+2 (placement.lua:490) can't
            --      be reused: pixels_to_cells already ceils both axes, so the
            --      ratio double-rounds and squashes wider expressions (`$k_{cat}$`
            --      lands at 2 cells when its true width is ~2.8).
            -- (conceallevel 0 never reaches here; find_visible filters it out.)
            local placement = require("snacks").image.placement
            local state = placement.state
            placement.state = function(self)
                local st = state(self)
                local r = self.opts.range
                if self.opts.type == "math" and st.loc.height == 1 and r and r[1] == r[3] then
                    local win = st.wins[1]
                    if win and vim.wo[win].conceallevel == 1 then
                        local src = vim.api.nvim_buf_get_text(self.buf, r[1] - 1, r[2], r[3] - 1, r[4], {})[1]
                        if src then
                            st.loc.width = math.max(1, vim.api.nvim_strwidth(src) - 1)
                        end
                    else
                        local sz = self.img.info and self.img.info.size
                        if sz then
                            local cell = require("snacks").image.terminal.size()
                            st.loc.width = require("utils.inline_math").cell_width(
                                sz.width, sz.height, cell.cell_width, cell.cell_height)
                        end
                    end
                end
                return st
            end

            -- Inline-math width depends on conceallevel (see above), but snacks
            -- only recomputes placements on scroll/edit/enter. Re-fire the inline
            -- manager's own BufWinEnter handler (scoped to its augroup, so no
            -- other BufWinEnter autocmds run) when conceallevel changes.
            vim.api.nvim_create_autocmd("OptionSet", {
                pattern = "conceallevel",
                callback = function()
                    local buf = vim.api.nvim_get_current_buf()
                    if vim.b[buf].snacks_image_attached then
                        pcall(vim.api.nvim_exec_autocmds, "BufWinEnter", {
                            buffer = buf,
                            group = "snacks.image.inline." .. buf,
                        })
                    end
                end,
            })

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
