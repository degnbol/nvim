local util = require "utils/init"
local read_lines = require("utils.file").read_lines
local map = require "utils/keymap"
local lsp = require "utils.lsp"
local paths = require "utils.paths"
local ts = require "utils/treesitter"

require "keymaps/options"
require "keymaps/danglish"
require "keymaps/surround"
require "keymaps/blockim"
require "keymaps/comments"

local hscroll = require("utils/hscroll").hscroll

for _, i in ipairs({ "", "2-", "3-", "4-" }) do
    for _, d in ipairs({ "Left", "Right" }) do
        local rhs = hscroll(d == "Left" and "h" or "l")
        local k = i .. "ScrollWheel" .. d
        map({ 'n', 'v', 'o', 'i' }, '<' .. k .. '>', rhs)
        map({ 'n', 'v', 'o', 'i' }, '<S-' .. k .. '>', rhs)
    end
    for _, d in ipairs({ "Up", "Down" }) do
        local k = i .. "ScrollWheel" .. d
        -- Shift+vertical scroll = horizontal scroll; multi-click variants just pass through
        if i == "" then
            map({ 'n', 'v', 'o', 'i' }, '<S-' .. k .. '>', hscroll(d == "Up" and "h" or "l"))
        else
            map({ 'n', 'v', 'o', 'i' }, '<S-' .. k .. '>', '<' .. k .. '>', { remap = true })
        end
    end
end

map({ 'n', 'v' }, '<Space>', '<Nop>', { silent = true })

-- Visual za: close all folds in selection if any are open, else open all.
-- Plain zo/zc are built-in in visual mode but neither handles the mixed case.
map.x('za', function()
    local s, _, e, _ = util.get_visual_range()
    local any_open = false
    for lnum = s, e do
        if vim.fn.foldlevel(lnum) > 0 and vim.fn.foldclosed(lnum) == -1 then
            any_open = true
            break
        end
    end
    vim.cmd(("%d,%dfold%s!"):format(s, e, any_open and "close" or "open"))
end, "Toggle folds in selection")

map.n('<leader>cc', function()
	local file = vim.api.nvim_buf_get_name(0)
	if file == '' then
		vim.notify("No file to run", vim.log.levels.WARN)
		return
	end
	local first_line = vim.api.nvim_buf_get_lines(0, 0, 1, false)[1] or ''
	local has_shebang = first_line:match('^#!') ~= nil
	if not has_shebang and vim.b.interpreter == nil then
		vim.notify("No shebang and no b:interpreter for filetype: " .. vim.bo.filetype,
			vim.log.levels.WARN)
		return
	end
	vim.cmd.write({ mods = { silent = true } })
	-- A shebang wins over b:interpreter: it may name a wrapper such as `uv run`
	-- that supplies deps the bare interpreter can't see.
	local interpreter = nil
	if has_shebang then
		if not vim.uv.fs_access(file, 'X') then
			-- Mask off the file-type bits; chmod is only specified for 07777.
			local mode = bit.band(vim.uv.fs_stat(file).mode, tonumber('7777', 8))
			local ok, err = vim.uv.fs_chmod(file, bit.bor(mode, tonumber('111', 8)))
			if not ok then
				vim.notify("chmod +x failed: " .. err, vim.log.levels.ERROR)
				return
			end
		end
	else
		interpreter = vim.b.interpreter
	end
	vim.cmd('!' .. util.script_command(file, interpreter))
end, "Run script")


-- in the terminal map escape to changing from terminal mode (insert mode) to
-- normal terminal mode <C-\><C-n> then change window to the left assuming that
-- term is on the right
map('t', "<Esc>", [[<C-\><C-n><C-w>h]])

-- Remap for dealing with word wrap
map.i('<down>', [[v:count == 0 ? '<C-\><C-O>gj' : '<down>']], "gj", { expr = true, silent = true })
map.i('<up>', [[v:count == 0 ? '<C-\><C-O>gk' : '<up>']], "gk", { expr = true, silent = true })
map.n('<down>', [[v:count == 0 ? 'gj' : '<down>']], "gj", { expr = true, silent = true })
map.n('<up>', [[v:count == 0 ? 'gk' : '<up>']], "gk", { expr = true, silent = true })

-- Typos. I don't use command window much but I often press q: or q; when I mean :q
map.n('q:', ':q')
map.n('q;', ':q')
map.n('<C-;>', 'q:')

-- Typos.
vim.api.nvim_create_user_command("Q", "q", {})
vim.api.nvim_create_user_command("Qa", "qa<bang>", { bang = true })
vim.api.nvim_create_user_command("X", "x", {})
vim.api.nvim_create_user_command("WQ", "wq", {})
vim.api.nvim_create_user_command("Wq", "wq", {})
vim.api.nvim_create_user_command("Lw", "w", {})
-- abbrev instead of command since command has to start with uppercase
vim.cmd [[cnoreabbrev qq q]]

-- small hack to remove excess whitespace possible since iw also captures
-- whitespace under cursor.
map.n("di ", "ciw <Esc>", "Delete excess whitespace")

map.n('[e', "<Cmd>move -2<CR>==", "Exchange line up", { silent = true })
map.n(']e', "<Cmd>move +1<CR>==", "Exchange line down", { silent = true })

map.n('>p', "p'[V']>", "Put below + incr indent")
map.n('>P', "P'[V']>", "Put above + incr indent")
map.n('<p', "p'[V']<", "Put below + decr indent")
map.n('<P', "P'[V']<", "Put above + decr indent")
map.n('=p', "p'[V']=", "Put below + reindent")
map.n('=P', "P'[V']=", "Put above + reindent")

-- like =p but for substitution
map.n('=ss', 'ss=`]', "Substitute+reindent", { remap = true, silent = true })
map.n("<leader>Sr", "<cmd>source $XDG_CONFIG_HOME/nvim/after/plugin/luasnip.lua<CR>", "Reload snippets")

-- Snippet management (replaces nvim-scissors)
map.x("<leader>xa", function() require("luasnippets/add").add_from_visual() end, "Snippet: Add")
map.n("<leader>xe", function() require("luasnippets/add").edit() end, "Snippet: Edit")

-- Built-in vim.snippet jump (same keys as blink.cmp snippet_forward/backward).
-- Blink's mappings take priority when its snippet session is active; these handle
-- vim.snippet sessions (e.g. from BufNewFile templates).
-- Built-in vim.snippet jump (same keys as blink.cmp snippet_forward/backward).
-- Blink's mappings take priority when its snippet session is active; these handle
-- vim.snippet sessions (e.g. from BufNewFile templates).
-- Uses expr + <Cmd> (like neovim's default Tab mapping) so select mode is preserved.
local function snippet_jump_expr(direction)
    return function()
        if vim.snippet.active({ direction = direction }) then
            return string.format('<Cmd>lua vim.snippet.jump(%d)<CR>', direction)
        end
        return ''
    end
end
vim.keymap.set({ 'i', 's' }, '<C-.>', snippet_jump_expr(1), { expr = true, silent = true, desc = "snippet_forward" })
vim.keymap.set({ 'i', 's' }, '<C-,>', snippet_jump_expr(-1), { expr = true, silent = true, desc = "snippet_backward" })
map.n('<C-.>', function()
    if vim.snippet.active({ direction = 1 }) then vim.snippet.jump(1) end
end, "snippet_forward")
map.n('<C-,>', function()
    if vim.snippet.active({ direction = -1 }) then vim.snippet.jump(-1) end
end, "snippet_backward")

-- use the following two commands to enable spelling
-- setlocal spell
-- set spelllang=en_us
-- ctrl+s to fix last misspelled word
-- credit: https://castel.dev/post/lecture-notes-1/
-- with git: https://github.com/gillescastel/latex-snippets
-- "i_ctrl-g u" = something about undo history, probably to keep the edit part of the insert edit?
-- [s = goto last misspelled, 1z= = replace with first spell suggestion.
local function spell_correct_closest()
    -- enable spell temporarily to spell correct
    if not vim.opt_local.spell:get() then
        vim.opt_local.spell = true
        vim.schedule(function()
            vim.opt_local.spell = false
        end)
    end
    return [[<c-g>u<Esc>[s1z=ea<c-g>u]]
end
-- using <C-s> for spell suggestion completion, see blink.cmp config
map.i('<C-S-s>', spell_correct_closest, "Spell correct closest", { expr = true })
map.n('<C-S-s>', [=[:set spell<CR>[s1z=e]=], "Spell correct closest")

-- not the most elegant but it works.
-- LeftMouse to move cursor to pressed location.
-- Then set @/ to the current word (\< and \> are to search strictly).
-- Then enable hlsearch. This is all a way to search without going to the next
-- match (if we just pressed * for instance)
map.n('<RightMouse>',
    [[<LeftMouse>:let @/='\<'.expand('<cword>').'\>'|set hlsearch<CR>]],
    "Search pressed word", { silent = true }
)

-- fallback search replace if both treesitter and LSP are not attached.
map.n('<leader>rn', [[:%s/<C-r><C-w>/]], "Search/replace cword")
-- use a selection that isn't a perfect cword, or just to use the simple search/replace when LSP is attached etc.
map.x('<leader>rn', [["ry:%s/<C-r>r/]], "Search/replace")

map.n('<Esc>', function()
    vim.cmd.nohlsearch()
    require('mini.notify').clear()
end, "Clear highlights and notifications", { silent = true })

-- <sa = my keybind for enable setting autoformat
-- gww = autoformat line
-- 0   = goto column 0, so we scroll all the way back to the right
-- gi  = go to last insert location and enter insert mode. Works even with the change to the line.
map.i('<C-S-A>', "<Esc><sagww0gi", "Enable autoformat and apply it", { remap = true })

-- map.n("<leader>qo", "<Cmd>copen<CR>", "Open")
-- map.n("<leader>qq", "<Cmd>cclose<CR>",  "Close") -- q for quit and is fast
-- Trying out "quicker.nvim" alt to stock quickfix
map.n("<leader>qo", function() require("quicker").open() end, "Open quicker")
map.n("<leader>qq", function() require("quicker").close() end, "Close quicker")
map.n("<leader>Qo", function() require("quicker").open({ loclist = true }) end, "Open quicker loclist")
map.n("<leader>Qq", function() require("quicker").close({ loclist = true }) end, "Close quicker loclist")
map.n("<leader>q1", "<Cmd>cc 1<CR>", "Entry 1")
map.n("<leader>q2", "<Cmd>cc 2<CR>", "Entry 2")
map.n("<leader>q3", "<Cmd>cc 3<CR>", "Entry 3")
-- we don't map :cnext etc here since we have ]q etc

map.n("<leader>QO", "<Cmd>lopen<CR>", "Open (stock)")
map.n("<leader>QQ", "<Cmd>lclose<CR>", "Close") -- Q for quit and is fast
map.n("<leader>Q1", "<Cmd>ll 1<CR>", "Entry 1")
map.n("<leader>Q2", "<Cmd>ll 2<CR>", "Entry 2")
map.n("<leader>Q3", "<Cmd>ll 3<CR>", "Entry 3")
-- we don't map :lnext etc here since we have ]l etc

map.n('<leader>bd', function() require "mini.bufremove".delete(vim.v.count) end, "bdel without win close")
map.n('<leader>bD', function() require "mini.bufremove".delete(vim.v.count, true) end, "bdel! without win close")
map.n('<leader>bw', function() require "mini.bufremove".wipeout(vim.v.count) end, "bwipe without win close")
map.n('<leader>bW', function() require "mini.bufremove".wipeout(vim.v.count, true) end, "bwipe! without win close")
map.n('<leader>bu', function() require "mini.bufremove".unshow(vim.v.count) end, "unshow without win close")
map.n('<leader>bn', "<Cmd>enew<CR>", "New")
map.n('<leader>bc', "<Cmd>tabclose<CR>", "tabclose")

-- poor fix for cmp replacing capitlisation of buffer words.
-- This long form is used over e.g. b~ea since the shorter form doesn't work for 1 or 2 char long words.
map.i("<C-`>", "<Esc>viwo<Esc>~gvo<Esc>a", "Capitalise last word")
map.n("<C-`>", "viwo<Esc>~gvo<Esc>", "Capitalise last word")

-- tired for accidentally jumping really far when pressing shift+down
-- Getting mapped by multicursor instead
map.x("<S-down>", "<down>")
map.x("<S-up>", "<up>")
-- map.n("<S-down>", "v<down>")
-- map.n("<S-up>", "v<up>")
-- shift+up and down jumping way to far for anything that would make sense in insert mode.
-- Changed to start/end of line but could do other things too.
map.i("<S-up>", "<C-o>^")
map.i("<S-down>", "<C-o>$")

map.i("<C-l>", "<right>", "Right")

-- Mac bindings.
map.c('<A-left>', "<s-left>", "move back one word")
map.c('<A-right>', "<s-right>", "move forward one word")
map({ 'n', 'v', 'c', 'i' }, '<D-left>', "<home>", { desc = "Start of line" })
map({ 'n', 'v', 'c', 'i' }, '<D-right>', "<end>", { desc = "Start of line" })
map.c('<A-BS>', "<C-w>", "Delete back one word")
-- Doesn't work well since we go one WORD to the right but only delete one word back.
map.c('<A-delete>', "<S-right><C-w>", "Delete next word")
-- this ignored if kitty handles it.
map.c('<D-BS>', "<C-u>", "Delete to beginning of line")

---Cmdline map typing `pattern` in a `/` or `?` search, and `other` elsewhere.
---@param lhs string
---@param pattern string
---@param other string
---@param desc string
local function search_key(lhs, pattern, other, desc)
    map.c(lhs, function()
        local cmdtype = vim.fn.getcmdtype()
        return (cmdtype == "/" or cmdtype == "?") and pattern or other
    end, desc, { expr = true, replace_keycodes = false })
end
-- Multi-line search: words separated by any whitespace, line breaks included.
search_key('<M-Space>', [[\_s\+]], " ", "Search: whitespace incl. line breaks")
search_key('<M-S-Space>', [[\_W\+]], " ", "Search: non-word run incl. line breaks")
-- Left Option+Space types U+00A0 in kitty (only right Option is Alt).
search_key('<Char-160>', [[\_s\+]], "\194\160", "Search: whitespace incl. line breaks")

map.x('<D-c>', '"+y', "Copy selection to clipboard")

map.n('<D-v>', 'p<C-=>', "Paste after, auto-indent, place cursor after", { remap = true })
map.n('<S-D-v>', 'P<C-=>', "Paste before, auto-indent, place cursor after", { remap = true })
map.i('<D-v>', function()
    local clipboard = vim.fn.getreg('+')
    -- whether part of line vs one or more whole lines
    local linewise = clipboard:match('\n')
    -- disregard blank lines and split
    local lines = vim.split(clipboard:gsub('\n*$', ''), '\n')
    -- insert charwise at cursor location placing cursor after
    vim.api.nvim_put(lines, 'c', false, true)
    -- auto-indent pasted lines if whole lines were pasted
    if linewise then
        -- TODO: when autoindent increases indent the c is off
        -- remember cursor location since I don't know a way to auto indent without moving cursor (in all cases)
        local r, c = unpack(vim.api.nvim_win_get_cursor(0))
        vim.cmd.normal "=`["
        vim.api.nvim_win_set_cursor(0, { r, c })
    end
end, "Paste, auto-indent, place cursor after")
map.c('<D-v>', '<C-r>+', "Paste, place cursor after")

-- treesitter mappings.
map.n('<leader>th', function()
    if vim.treesitter.highlighter.active[vim.api.nvim_get_current_buf()] then
        vim.treesitter.stop()
    else
        vim.treesitter.start()
    end
end, "Toggle highlight")
map.n('<leader>ti', vim.show_pos, "Inspect") -- Same as :Inspect
map.n('<leader>tI', "<Cmd>checkhealth nvim-treesitter<CR>", "Treesitter health")
map.n('<leader>tt', "<Cmd>InspectTree<CR>", "Inspect tree")
map.n('<leader>tn', function()
    local node = vim.treesitter.get_node()
    if node == nil then
        print("No node found")
    else
        local text = vim.treesitter.get_node_text(node, 0)
        local type = node:type()
        print(text, "type =", type)
    end
end, "node")
map.n('<leader>tN', function()
    local node = vim.treesitter.get_node():parent()
    if node == nil then
        print("No parent node found")
    else
        local text = vim.treesitter.get_node_text(node, 0)
        local type = node:type()
        print(text, "type =", type)
    end
end, "parent")

-- window layout.
vim.api.nvim_create_autocmd("FileType", {
    pattern = "*",
    group = vim.api.nvim_create_augroup("WindowLayout", { clear = true }),
    callback = function()
        local ftapp = {
            tex = "skim",
            -- pymol
            python = "/opt/homebrew/Caskroom/miniforge/base/envs/pymol/bin/python",
        }
        local rectangle = function(layouts)
            local cmd = { "rectangle" }
            for layout, app in pairs(layouts) do
                if app ~= nil then
                    table.insert(cmd, app)
                    table.insert(cmd, layout)
                end
            end
            return function()
                vim.system(cmd, { timeout = 1500 }, function(obj)
                    if obj.code ~= 0 then
                        util.schedule_notify(obj)
                    end
                end)
            end
        end
        local this = "kitty"
        local other = ftapp[vim.bo.filetype]
        map.n('<LocalLeader>1', rectangle { maximize = this },
            "Whole screen layout", { buffer = true })
        map.n('<LocalLeader>2', rectangle { ["right-half"] = other, ["left-half"] = this },
            "Half screen layout", { buffer = true })
        map.n('<LocalLeader>3', rectangle { ["last-third"] = other, ["first-two-thirds"] = this },
            "Two-thirds screen layout", { buffer = true })
    end
})

-- Like doing new | r!<CMD> except the scratch buffer is wiped when hidden.
map.n("<leader>:!", function()
    vim.ui.input({}, function(cmd)
        if cmd and cmd ~= "" then
            vim.cmd("noswapfile new")
            vim.bo.buftype = "nofile"
            vim.bo.bufhidden = "wipe"
            vim.api.nvim_buf_set_lines(0, 0, -1, false, vim.fn.systemlist(cmd))
        end
    end)
end, "new|r!<CMD> with bh=wipe")


-- LSP and completion status, overall conf etc.
map.n("<leader>li", "<Cmd>checkhealth vim.lsp<CR>", "Info")
-- can't use backspace since it is hardcoded by mini.clue for up one level
map.n("<leader>l0", "<Cmd>lsp stop<CR>", "Stop")
-- Not `:lsp enable`: with no name it enables every config for the filetype.
map.n("<leader>l1", function()
    local configs = vim.lsp.get_configs { enabled = true, filetype = vim.bo.filetype }
    vim.lsp.enable(vim.tbl_map(function(config) return config.name end, configs))
end, "Enable")
map.n("<leader>l!", "<Cmd>lsp restart<CR>", "Restart")
map.n("<leader>lL", function()
    vim.cmd.edit(vim.lsp.log.get_filename())
end, "Log")
-- match other completion related entries under x
map.n("<leader>xS", "<Cmd>CmpStatus<CR>", "Cmp status")

-- See `:help vim.diagnostic.*` for documentation on any of the below functions
map.n('<leader>dd', vim.diagnostic.open_float, "Line diagnostic")
map.n('<leader>dv', function() vim.diagnostic.config { virtual_lines = true } end, "Enable virtual line diagnostics")
map.n('<leader>dV', function() vim.diagnostic.config { virtual_lines = false } end, "Enable virtual line diagnostics")
-- can't use backspace since it is hardcoded by mini.clue for up one level
map.n('<leader>d0', function() vim.diagnostic.enable(false) end, "Disable diagnostics")
map.n('<leader>d1', vim.diagnostic.enable, "Enable diagnostics")
-- Simplify: remove "forward/backward".
map.desc('n', '[d', "Diagnostic")
map.desc('n', ']d', "Diagnostic")
map.n('<leader>dl', vim.diagnostic.setloclist, "Loclist diagnostics")

map.n('grr', lsp.references, "References")
map.desc('n', 'gra', "Code actions")
map.desc('n', 'gri', "Implementations")
-- grn (rename) is handled by live-rename.nvim's lz.n keys spec.
map.desc('n', 'grt', "Type definitions")
map.n('grd', lsp.definition, "Definition")

map.i('<C-s>', function()
    -- TODO: modify float to remove empty lines at top and bottom.
    -- TODO: update signature help when pressing comma or deleting a comma or moving cursor.
    -- Also decide if repeated <C-s> should cycle the signatures, as is default.
    local call_expression = ts.ancestor('call_expression')
    if call_expression == nil then return end
    local r, c = unpack(vim.api.nvim_win_get_cursor(0))
    local start_row, start_col, end_row, end_col = vim.treesitter.get_node_range(call_expression)
    vim.lsp.buf.signature_help { title = nil, offset_x = start_col - c, close_events = { "WinScrolled", "ModeChanged" } }
end, "Signature help")

-- LSP hover/signature markdown renders in a conceallevel=3 float (see
-- ftplugin/markdown.lua), which fully hides conceal-replaced text — including
-- the HTML entities nvim's markdown query resolves to a glyph at level 2.
-- Decode them in the LSP markup before it becomes float lines so they survive:
-- a rendered hover should show `<`, not nothing. (peek_file below passes raw
-- file lines straight to open_floating_preview, bypassing this — file peeks
-- stay verbatim.)
local html_entities = {
    ["&lt;"] = "<", ["&gt;"] = ">", ["&amp;"] = "&", ["&quot;"] = '"',
    ["&nbsp;"] = " ", ["&ensp;"] = " ", ["&emsp;"] = " ",
}
local convert_markdown = vim.lsp.util.convert_input_to_markdown_lines
---@diagnostic disable-next-line: duplicate-set-field
vim.lsp.util.convert_input_to_markdown_lines = function(input, contents)
    local lines = convert_markdown(input, contents)
    for i, line in ipairs(lines) do
        lines[i] = line:gsub("&%a+;", html_entities)
    end
    return lines
end

local PEEK_LINES = 500

-- Marks the peek float for open_floating_preview, which uses it to focus a
-- float that is already up instead of building a second one — the KK round
-- trip that :h vim.lsp.buf.hover gets from focus_id = "textDocument/hover".
local PEEK_FOCUS = "peek_file"

--- Peek a file in a hover-style float, opening on `lnum` and showing
--- PEEK_LINES lines from there — more than fills the float. An out-of-reach
--- `lnum`, whether past EOF (a stale grep hit) or too deep to read, warns and
--- falls back to the top. Highlighting is the matched filetype as the float's
--- 'syntax'. ponytail: regex syntax is enough — add vim.treesitter.start if it
--- ever matters.
---@param path string
---@param lnum integer|nil 1-based line to show first, default the top
local function peek_file(path, lnum)
    local lines = read_lines(path, lnum or 1, PEEK_LINES)
    local missed = lnum ~= nil and #lines == 0
    if missed then lines = read_lines(path, 1, PEEK_LINES) end
    -- A zero-byte file: open_floating_preview throws "Invalid 'width'" on it.
    if #lines == 0 then
        vim.notify(("%s is empty"):format(vim.fs.basename(path)), vim.log.levels.WARN)
        return
    end
    if missed then
        vim.notify(("%s: line %d is past the end of the file or the read cap")
            :format(vim.fs.basename(path), lnum), vim.log.levels.WARN)
    end
    -- The compression is not the content: a .tsv.gz highlights as a tsv.
    local ft = vim.filetype.match({ filename = (path:gsub("%.gz$", "")) }) or ""
    local buf, win = vim.lsp.util.open_floating_preview(lines, ft, { focus_id = PEEK_FOCUS })
    -- The float was already up and open_floating_preview moved the cursor into
    -- it rather than filling a new buffer; its lines are these ones already.
    if win == vim.api.nvim_get_current_win() then return end
    -- open_floating_preview rewrites its contents — trimempty for a plain
    -- syntax, _normalize_markdown for markdown — either of which would shift
    -- `lnum` off the top of the float.
    vim.bo[buf].modifiable = true
    vim.api.nvim_buf_set_lines(buf, 0, -1, false, lines)
    vim.bo[buf].modifiable = false
end

local function lsp_hover_capable(bufnr)
    for _, client in ipairs(vim.lsp.get_clients({ bufnr = bufnr })) do
        if client:supports_method("textDocument/hover") then return true end
    end
    return false
end

-- Overload K: if the cursor is on a path resolving to a readable file, peek it;
-- else LSP hover; else built-in keywordprg. Global (not buffer-local on
-- LspAttach) so it works in no-LSP buffers — the primary use case — and so the
-- LspAttach default sees a K mapping and skips installing its own.
map.n("K", function()
    -- Inside the peek float, K hands focus back to the window it was opened
    -- from, closing the KK round trip. open_floating_preview does this itself,
    -- but only on its way to a float it would then refill with the wrong lines.
    if vim.w[0][PEEK_FOCUS] then
        vim.cmd("wincmd p")
        return
    end
    local path, lnum = paths.resolve_location_under_cursor(0)
    if path and paths.is_file(path) then -- dirs (case 1) fall to hover
        peek_file(path, lnum)
        return
    end
    if lsp_hover_capable(0) then
        local cword_start = util.cword_start_col()
        local c = vim.api.nvim_win_get_cursor(0)[2]
        vim.lsp.buf.hover({ offset_x = cword_start - c - 1 }) -- align float to word start
        return
    end
    -- n = noremap (bypass this + buffer-local maps → built-in K), x = execute now.
    vim.api.nvim_feedkeys(vim.v.count1 .. "K", "nx", false)
end, "Peek file / hover / keywordprg")

local GX_DESC = "Smart URL opener"
-- The runtime gx: LSP document links, url extmarks, treesitter url metadata
-- and <cfile>. Checked by desc since re-sourcing this file would otherwise
-- capture the map below.
local default_gx = vim.fn.maparg("gx", "n", false, true)
default_gx = default_gx.desc ~= GX_DESC and default_gx.callback or nil
if not default_gx then
    vim.notify("gx: no runtime default gx to fall back to", vim.log.levels.WARN)
end

---Open a URL with the system handler, notifying if that fails.
---@param url string
local function open_url(url)
    local _, err = vim.ui.open(url)
    if err then vim.notify(err, vim.log.levels.ERROR) end
end

-- custom gx function that supports more website links.
map.n("gx", function()
    -- Go to github for plugins easily.
    if vim.bo.filetype == "lua" then
        -- First check if we are editing a plugin spec file
        local config = vim.fn.stdpath("config")
        local filepath = vim.api.nvim_buf_get_name(0)
        if filepath:match(config .. '/lua/plugins/') then
            -- get first string on the line, assumes we don't list multiple plugins on one line.
            local line = vim.api.nvim_get_current_line()
            local repo = line:match([["([%w._-]+)"]]) or line:match([['([%w._-]+)']])
            if repo then
                if repo:match("http") then
                    return open_url(repo)
                elseif repo:match("/") then
                    -- already an account/repo slug
                    return open_url("https://github.com/" .. repo)
                end
                -- Post vim.pack migration lz.n specs carry only the bare repo
                -- name. Resolve the source from pack_specs.lua.
                local specs = vim.fn.readfile(config .. "/lua/pack_specs.lua")
                local url = require("utils.pluginspec").resolve(repo, specs)
                if url then return open_url(url) end
                -- unresolved: fall through to the general URL resolvers below
            end
        end
    end

    -- for latex packages...
    if vim.bo.filetype == "tex" then
        local line = vim.api.nvim_get_current_line()
        local pac = line:match("\\usepackage.*{([%w_-]+)}")
        if pac ~= nil then
            local ctan = "https://ctan.org/pkg/" .. pac .. "?lang=en"
            return open_url(ctan)
        end
    end

    -- DOI links (e.g. 10.1234/abc, doi:10.1234/abc)
    local line = vim.api.nvim_get_current_line()
    local doi = line:match("10%.%d%d%d%d+/[%w%.%-_/:]+[%w]")
    if doi then
        return open_url("https://doi.org/" .. doi)
    end

    local row, col = unpack(vim.api.nvim_win_get_cursor(0))
    local link = lsp.document_link_at(0, row - 1, col)
    if link then return open_url(link) end

    -- Using various-textobjs url finder, we will detect url further ahead,
    -- e.g. useful when being lazy and the url is right there on the line.
    -- visually select URL
    if default_gx then
        -- The default takes over on a miss, so silence the textobj's miss message.
        local notify = require("various-textobjs.config.config").config.notify
        local when_not_found = notify.whenObjectNotFound
        notify.whenObjectNotFound = false
        local ok, err = pcall(require("various-textobjs").url)
        notify.whenObjectNotFound = when_not_found
        if not ok then error(err, 0) end
    else
        require("various-textobjs").url()
    end
    -- plugin only switches to visual mode when a URL is found (and notifies on miss)
    if not util.get_mode():find("^v") then
        if default_gx then default_gx() end
        return
    end
    -- retrieve URL with the z-register as intermediary
    vim.cmd.normal { '"zy', bang = true }
    local url = vim.fn.getreg("z")
    open_url(url)
end, GX_DESC)

-- Clear the builtin C-leftclick which is goto tag def.
-- It can't be removed with vim.keymap.del since it's builtin.
-- We always just use C-] for explicit tag lookup and C-leftmouse is mapped from Cmd-leftmouse in arch.
map.n("<C-leftmouse>", "")


vim.api.nvim_create_user_command('MessagesCopy', function()
    local output = vim.api.nvim_exec2('messages', { output = true }).output
    output = vim.trim(output)
    if output == '' then
        vim.notify('No messages to copy', vim.log.levels.WARN)
        return
    end
    vim.fn.setreg('+', output)
    vim.notify('Messages copied (' .. #vim.split(output, '\n') .. ' lines)')
end, {})
