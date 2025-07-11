#!/usr/bin/env lua
local util = require "utils/init"
local map = vim.keymap.set

-- shift should have no effect on scroll
local counts = { "", "2-", "3-", "4-" }
local directions = { "Up", "Down", "Left", "Right" }
for _, i in ipairs(counts) do
    for _, d in ipairs(directions) do
        local k = i .. "ScrollWheel" .. d
        map({ 'n', 'v', 'o', 'i' }, '<S-' .. k .. '>', '<' .. k .. '>', { remap = true })
    end
end

-- Keymaps for better default experience
map({ 'n', 'v' }, '<Space>', '<Nop>', { silent = true })

-- in the terminal map escape to changing from terminal mode (insert mode) to
-- normal terminal mode <C-\><C-n> then change window to the left assuming that
-- term is on the right
map('t', "<Esc>", [[<C-\><C-n><C-w>h]])

-- Danglish support. For when Danglish keyboard is selected.
-- Generally you should instead stay in code keyboard and use iminsert=2
-- This can also be done with langmap but since these are
-- never used in normal mode then it doesn't hurt, and I tried and it didn't
-- seem to map to [ with remapping even though langremap was on so I don't know
map({ 'n', 'x' }, "æ", ";", { remap = true })
map({ 'n', 'x' }, "Æ", ":", { remap = true })
map({ 'n', 'x' }, "ø", "'", { remap = true })
map({ 'n', 'x' }, "Ø", '"', { remap = true })
map({ 'n', 'x' }, "å", "[", { remap = true })
map({ 'n', 'x' }, "Å", "{", { remap = true })
-- In Danglish I moved : and ; to the }] button
-- But this messes with things when I'm not in Danglish.

-- Mapping for function keys available on mechanical keyboard
-- TODO: too far away for being useful for something like regular completion,
-- tab is more convenient. Find another use, e.g. snippet
-- completion and navigation? Then you can free up ctrl+j/k/l although we're
-- not always on mech keyboard, and I don't mind using <s-c-K> for the rare
-- case of needing special chars.
map('i', "<F13>", "<C-j>", { remap = true })
map('i', "<F14>", "<C-k>", { remap = true })
map('i', "<F15>", "<C-l>", { remap = true })

-- When language is set to danish we can get ;:'" with alt the same way as we
-- do in danglish keyboard input, however that is with the left alt which is
-- deeper in the OS and the mapping below is for when esc+key (^[) is detected which
-- is what is sent to the terminal, e.g. with kitty's setting 'macos_option_as_alt right'.
-- Note that they are also available with the ]} and \| keys.
map('i', "<A-;>", ";")
map('i', "<A-S-;>", ":")
map('i', "<A-'>", "'")
map('i', "<A-S-'>", '"')
-- when writing text with Danish we might try to write : but the key is mapped to Æ.
-- : is written at ends of word where we would never write capital Æ, so we can check if we are at end of word.
-- The only exception would be if the entire word is uppercase. Currently choosing to ignore that edge case.
map('i', "Æ", function()
    local char = util.get_current_char()
    local put = char:match('[A-Åa-å0-9.,!?]') and ':' or 'Æ'
    util.put_char(put)
end)
-- similarly we often would want " instead of Ø, e.g. if we write ØØ it's to make "" and
map('i', "Ø", function()
    local r, c = util.get_cursor()
    local char, c1 = util.get_char(r, c)
    if char == 'Ø' then
        vim.api.nvim_buf_set_text(0, r, c1, r, c, { '""' })
        -- try to see if it's ok, otherwise change to only looking for uppercase.
    elseif util.get_char(r, c + 1):match('[A-Åa-å]') then
        util.put_char('"')
    else
        local put = char:match('[A-Åa-å.,!?]') and '"' or 'Ø'
        util.put_char(put)
    end
end)
-- for å we might want [] or {}, but with <C-6> pressed they're mapped to å; and Å:,
-- however it's a remap from ] and } so we write that.
map('i', 'å]', '[]')
map('i', 'Å}', '{}')

local function toggle_danglish(silent)
    if not silent then
        if vim.bo.iminsert == 0 then
            print("Danglish")
        else
            print("Code")
        end
    end
    if util.get_mode() == 'n' then
        util.press("a<C-^><Esc>")
    else
        util.press("<C-^>")
    end
    -- other related things to get triggered when toggling language.
    -- Has to be scheduled since it checks iminsert which isn't updated immediately.
    vim.schedule(function()
        vim.api.nvim_exec_autocmds("User", { pattern = "ToggleDansk" })
    end)
end

local function toggle_danish_imaps(silent)
    local is_mapped = vim.g.danish_imaps
    if not is_mapped then
        vim.g.danish_imaps = true
        vim.keymap.set('i', 'ae', 'æ', { desc = "ae -> æ" })
        vim.keymap.set('i', 'oe', 'ø', { desc = "oe -> ø" })
        vim.keymap.set('i', 'aa', 'å', { desc = "aa -> å" })
        vim.keymap.set('i', 'Ae', 'Æ', { desc = "Ae -> Æ" })
        vim.keymap.set('i', 'Oe', 'Ø', { desc = "Oe -> Ø" })
        vim.keymap.set('i', 'Aa', 'Å', { desc = "Aa -> Å" })
        vim.keymap.set('i', 'AE', 'Æ', { desc = "AE -> Æ" })
        vim.keymap.set('i', 'OE', 'Ø', { desc = "OE -> Ø" })
        vim.keymap.set('i', 'AA', 'Å', { desc = "AA -> Å" })
        if not silent then print("ae -> æ, osv.   aktiveret") end
    else
        vim.g.danish_imaps = false
        vim.keymap.del('i', 'ae')
        vim.keymap.del('i', 'oe')
        vim.keymap.del('i', 'aa')
        vim.keymap.del('i', 'Ae')
        vim.keymap.del('i', 'Oe')
        vim.keymap.del('i', 'Aa')
        vim.keymap.del('i', 'AE')
        vim.keymap.del('i', 'OE')
        vim.keymap.del('i', 'AA')
        if not silent then print("ae -> æ, osv. DEaktiveret") end
    end
end

map('i', '<C-6>', toggle_danglish, { desc = "Toggle dansk" })
-- Ins key was a bit useless just doing what i does so let's make it a language switch insertion:
map('n', '<Ins>', 'i<C-6>', { remap = true, desc = "Insert + Toggle dansk" })
-- switch to/from Danish æøå and to insert mode, which is convenient.
-- remap in order to utilise the remapped <C-6> which updates the cmp dictionary
map("n", "yod", "i<C-6>", { remap = true, desc = "Danish (<C-^>)" })

map('i', "<C-S-6>", toggle_danish_imaps, { desc = "Toggle Danish imaps" })
map('n', "<leader>D", function()
    toggle_danish_imaps(true)
    toggle_danglish(true)
    if vim.g.danish_imaps then
        print("Danglish + ae -> æ, osv.   aktiveret")
    else
        print("Code     + ae -> æ, osv. DEaktiveret")
    end
end, { desc = "Toggle Danish imaps + Danglish" })

-- And in case danglish keyboard is active:
map('i', "<A-æ>", ";")
map('i', "<A-S-æ>", ":")
map('i', "<A-ø>", "'")
map('i', "<A-S-ø>", '"')

-- Remap for dealing with word wrap
map('i', '<down>', [[v:count == 0 ? '<C-\><C-O>gj' : '<down>']], { expr = true, silent = true })
map('i', '<up>', [[v:count == 0 ? '<C-\><C-O>gk' : '<up>']], { expr = true, silent = true })
map('n', '<down>', [[v:count == 0 ? 'gj' : '<down>']], { expr = true, silent = true })
map('n', '<up>', [[v:count == 0 ? 'gk' : '<up>']], { expr = true, silent = true })

-- Typos. I don't use command window much but I often press q: or q; when I mean :q
map('n', 'q:', ':q')
map('n', 'q;', ':q')
map('n', '<C-;>', 'q:')

-- Typos.
vim.api.nvim_create_user_command("Q", "q", {})
vim.api.nvim_create_user_command("X", "x", {})
vim.api.nvim_create_user_command("WQ", "wq", {})
vim.api.nvim_create_user_command("Wq", "wq", {})
vim.api.nvim_create_user_command("Lw", "w", {})
-- abbrev instead of command since command has to start with uppercase
vim.cmd [[cnoreabbrev qq q]]

map('i', "<C-'>", "''<left>")
map('i', "<C-S-'>", '""<left>')
map('i', "<C-9>", '()<left>')
map('i', "<C-S-9>", '()<left>')
map('i', "<C-0>", '()<left>')
map('i', "<C-S-0>", '()<left>')
-- map('i', "<C-[>", '[]<left>') -- not possible since it's literally ESC
map('i', "<C-]>", '[]<left>')
map('i', "<C-S-[>", '{}<left>')
map('i', "<C-S-]>", '{}<left>')

-- hack map of shift+space
local bracketJumpCode = "\x1F"
local paired = { '""', "''", "``", "()", '[]', "{}", "<>", "$$" }
local paireddouble = { '[[]]', "\\{\\}" }     -- lua and tex
local singles = { "'", '"', '`', '(', ')', '[', ']', '{', '}', '<', '>', '$' }
local doubles = { "[[", "]]", "\\{", "\\}", } -- lua and tex
local triples = { '"""', "'''", "```" }
local function bracketJump(line, c)
    if vim.tbl_contains(triples, line:sub(c - 2, c)) then
        return "<left><left><left>"
    elseif vim.tbl_contains(triples, line:sub(c + 1, c + 3)) then
        return "<right><right><right>"
    elseif vim.tbl_contains(paireddouble, line:sub(c - 3, c)) then
        return "<left><left>"
    elseif vim.tbl_contains(paireddouble, line:sub(c + 1, c + 4)) then
        return "<right><right>"
        -- left priority over right for pairs so we go back first for e.g. Matrix{}|[]
    elseif vim.tbl_contains(paired, line:sub(c - 1, c)) then
        return "<left>"
    elseif vim.tbl_contains(paired, line:sub(c + 1, c + 2)) then
        return "<right>"
    elseif vim.tbl_contains(doubles, line:sub(c + 1, c + 2)) then
        return "<right><right>"
    elseif vim.tbl_contains(doubles, line:sub(c - 1, c)) then
        return "<left><left>"
        -- right priority over left for singles, since they are usually half of a
        -- filled out pair and we want to prioritize progressing in that case.
    elseif vim.tbl_contains(singles, line:sub(c + 1, c + 1)) then
        return "<right>"
    elseif vim.tbl_contains(singles, line:sub(c, c)) then
        return "<left>"
    else -- maybe accidentally still held shift
        return " "
    end
end
map({ 'i', 'n' }, bracketJumpCode, function()
    local line = vim.api.nvim_get_current_line()
    local r, c = unpack(vim.api.nvim_win_get_cursor(0))
    return bracketJump(line, c)
end, { expr = true, desc = "Move inside empty pair/triples or outside non-empty" })
map('c', bracketJumpCode, function()
    local line = vim.fn.getcmdline()
    local c = vim.fn.getcmdpos()
    return bracketJump(line, c - 1)
end, { expr = true, desc = "Move inside empty pair/triples or outside non-empty" })

-- useful with cursor | in {|} to get
-- {
--     |
-- }
map('i', '<S-CR>', "<CR><Esc>O", { desc = "Indented newline" })

-- small hack to remove excess whitespace possible since iw also captures
-- whitespace under cursor.
map("n", "di ", "ciw <Esc>", { desc = "Delete excess whitespace" })

-- like =p but for substitution
map("n", '=ss', 'ss=`]', { remap = true, silent = true, desc = "Substitute+reindent" })
map("n", "<leader>Sr", "<cmd>source $XDG_CONFIG_HOME/nvim/after/plugin/luasnip.lua<CR>", { desc = "Reload snippets" })

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
map('i', '<C-S-s>', spell_correct_closest, { expr = true, desc = "Spell correct closest" })
map('n', '<C-S-s>', [=[:set spell<CR>[s1z=e]=], { desc = "Spell correct closest" })

-- not the most elegant but it works.
-- LeftMouse to move cursor to pressed location.
-- Then set @/ to the current word (\< and \> are to search strictly).
-- Then enable hlsearch. This is all a way to search without going to the next
-- match (if we just pressed * for instance)
map('n', '<RightMouse>',
    [[<LeftMouse>:let @/='\<'.expand('<cword>').'\>'|set hlsearch<CR>]],
    { silent = true, desc = "Search pressed word" }
)

-- fallback search replace if both treesitter and LSP are not attached.
map('n', '<leader>rn', [[:%s/<C-r><C-w>/]], { desc = "Search/replace cword" })
-- use a selection that isn't a perfect cword, or just to use the simple search/replace when LSP is attached etc.
map('x', '<leader>rn', [["ry:%s/<C-r>r/]], { desc = "Search/replace" })

-- <sa = my keybind for enable setting autoformat
-- gww = autoformat line
-- 0   = goto column 0, so we scroll all the way back to the right
-- gi  = go to last insert location and enter insert mode. Works even with the change to the line.
map('i', '<C-S-A>', "<Esc><sagww0gi", { remap = true, desc = "Enable autoformat and apply it." })

-- map('n', "<leader>qo", "<Cmd>copen<CR>", { desc = "Open" })
-- map('n', "<leader>qq", "<Cmd>cclose<CR>", { desc = "Close" }) -- q for quit and is fast
-- Trying out "quicker.nvim" alt to stock quickfix
map('n', "<leader>qo", function() require("quicker").open() end, { desc = "Open quicker", })
map('n', "<leader>qq", function() require("quicker").close() end, { desc = "Close quicker", })
map('n', "<leader>Qo", function() require("quicker").open({ loclist = true }) end, { desc = "Open quicker loclist", })
map('n', "<leader>Qq", function() require("quicker").close({ loclist = true }) end, { desc = "Close quicker loclist", })
map('n', "<leader>q1", "<Cmd>cc 1<CR>", { desc = "Entry 1" })
map('n', "<leader>q2", "<Cmd>cc 2<CR>", { desc = "Entry 2" })
map('n', "<leader>q3", "<Cmd>cc 3<CR>", { desc = "Entry 3" })
-- we don't map :cnext etc here since we have ]q etc

map('n', "<leader>Qo", "<Cmd>lopen<CR>", { desc = "Open" })
map('n', "<leader>QQ", "<Cmd>lclose<CR>", { desc = "Close" }) -- Q for quit and is fast
map('n', "<leader>Q1", "<Cmd>ll 1<CR>", { desc = "Entry 1" })
map('n', "<leader>Q2", "<Cmd>ll 2<CR>", { desc = "Entry 2" })
map('n', "<leader>Q3", "<Cmd>ll 3<CR>", { desc = "Entry 3" })
-- we don't map :lnext etc here since we have ]l etc

---Delete buffer. Repeat for unnamed empty buffers.
---@param opts table with bool key force (passed to vim.api.nvim_buf_delete)
---@param lastbufnr integer? for recursion
local function bufdel(opts, lastbufnr)
    opts = opts or {}
    local bufnr = vim.api.nvim_get_current_buf()
    -- make sure to not retry if a previous call failed
    if bufnr == lastbufnr then return end
    vim.api.nvim_buf_delete(0, opts)
    -- repeat if next buffer is empty (stop annoying behaviour of vimtex)
    if not util.is_named() and util.is_empty() then
        vim.schedule(function() bufdel(opts, bufnr) end)
    end
end

map('n', '<leader>bd', bufdel, { desc = "Delete" })
map('n', '<leader>bD', function() bufdel { force = true } end, { desc = "Delete!" })
map('n', "<leader>bn", "<Cmd>enew<CR>", { desc = "New" })
map('n', "<leader>bc", "<Cmd>tabclose<CR>", { desc = "tabclose" })

-- poor fix for cmp replacing capitlisation of buffer words.
-- This long form is used over e.g. b~ea since the shorter form doesn't work for 1 or 2 char long words.
map('i', "<C-`>", "<Esc>viwo<Esc>~gvo<Esc>a", { desc = "Capitalise last word" })
map('n', "<C-`>", "viwo<Esc>~gvo<Esc>", { desc = "Capitalise last word" })

-- tired for accidentally jumping really far when pressing shift+down
-- Getting mapped by multicursor instead
map('v', "<S-down>", "<down>")
map('v', "<S-up>", "<up>")
-- map('n', "<S-down>", "v<down>")
-- map('n', "<S-up>", "v<up>")
-- shift+up and down jumping way to far for anything that would make sense in insert mode.
-- Changed to start/end of line but could do other things too.
map('i', "<S-up>", "<C-o>^")
map('i', "<S-down>", "<C-o>$")

map('i', "<C-l>", "<right>", { desc = "Right" })

map('c', '<A-left>', "<s-left>", { desc = "move back one word" })
map('c', '<A-right>', "<s-right>", { desc = "move forward one word" })
map('c', '<A-BS>', "<C-w>", { desc = "Delete back one word" })
-- Doesn't work well since we go one WORD to the right but only delete one word back.
map('c', '<A-delete>', "<S-right><C-w>", { desc = "Delete next word" })
-- this ignored if kitty handles it.
map('c', '<D-BS>', "<C-u>", { desc = "Delete to beginning of line" })

map('n', '<D-v>', 'p<C-=>', { desc = "Paste after, auto-indent, place cursor after", remap = true })
map('n', '<S-D-v>', 'P<C-=>', { desc = "Paste before, auto-indent, place cursor after", remap = true })
map('i', '<D-v>', function()
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
end, { desc = "Paste, auto-indent, place cursor after" })
map('c', '<D-v>', '<C-r>+', { desc = "Paste, place cursor after" })

-- treesitter mappings.
map('n', '<leader>th', "<Cmd>TSBufToggle highlight<CR>", { desc = "Toggle local highlight" })
map('n', '<leader>tH', "<Cmd>TSToggle highlight<CR>", { desc = "Toggle global highlight" })
map('n', '<leader>ti', vim.show_pos, { desc = "Inspect" }) -- Same as :Inspect
map('n', '<leader>tt', "<Cmd>InspectTree<CR>", { desc = "Inspect tree" })
map('n', '<leader>tI', "<Cmd>Capture TSInstallInfo<CR>", { desc = "Install info" })
map('n', '<leader>tn', function()
    local node = vim.treesitter.get_node()
    local text = vim.treesitter.get_node_text(node, 0)
    local type = node:type()
    print(text, "type=", type)
end, { desc = "node" })
map('n', '<leader>tN', function()
    local node = vim.treesitter.get_node():parent()
    local text = vim.treesitter.get_node_text(node, 0)
    local type = node:type()
    print(text, "type=", type)
end, { desc = "parent" })

-- window layout.
local grp = vim.api.nvim_create_augroup("WindowLayout", { clear = true })
vim.api.nvim_create_autocmd("Filetype", {
    pattern = "*",
    group = grp,
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
        vim.keymap.set('n', '<LocalLeader>1',
            rectangle { maximize = this },
            { buffer = true, desc = "Whole screen layout" })
        vim.keymap.set('n', '<LocalLeader>2',
            rectangle { ["right-half"] = other, ["left-half"] = this },
            { buffer = true, desc = "Half screen layout" })
        vim.keymap.set('n', '<LocalLeader>3',
            rectangle { ["last-third"] = other, ["first-two-thirds"] = this },
            { buffer = true, desc = "Two-thirds screen layout" })
    end
})
