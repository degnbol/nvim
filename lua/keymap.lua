#!/usr/bin/env lua
local map = vim.keymap.set

-- shift should have no effect on scroll
counts = {"", "2-", "3-", "4-"}
directions = {"Up", "Down", "Left", "Right"}
for _, i in ipairs(counts) do
    for _, d in ipairs(directions) do
        k = i .. "ScrollWheel" .. d
        map({'n','v','o','i'}, '<S-' .. k .. '>', '<' .. k .. '>', {remap=true})
    end
end

-- Keymaps for better default experience
map({ 'n', 'v' }, '<Space>', '<Nop>', { silent = true })

-- use U instead of <C-r> for redo.
-- TODO move line undo to Vu 
-- map('n', "U", "<C-r>")

-- For remote with problematic clipboard we replace the cutlass x that works 
-- locally with a mapping where we call y (which is mapped to smartyank) then 
-- d.
map('n', "xx", "yydd")
map('x', "x",  "ygvd")
function CutOperator(type, ...)
    if type == "line" then
        vim.cmd.normal "`[V`]ygvd"
    else
        vim.cmd.normal "`[v`]ygvd"
    end
end
map('n', "x", function ()
    -- TODO: use function rather than string after pull request is merged:
    -- https://github.com/neovim/neovim/pull/20187
    vim.opt.operatorfunc = "v:lua.CutOperator"
    return "g@"
end, {expr=true})

-- in the terminal map escape to changing from terminal mode (insert mode) to 
-- normal terminal mode <C-\><C-n> then change window to the left assuming that 
-- term is on the right
map('t', "<Esc>", [[<C-\><C-n><C-w>h]])


-- Danglish support. For when Danglish keyboard is selected.
-- Generally you should instead stay in code keyboard and use iminsert=2
-- This can also be done with langmap but since these are 
-- never used in normal mode then it doesn't hurt, and I tried and it didn't 
-- seem to map to [ with remapping even though langremap was on so I don't know
map({'n', 'x'}, "æ", ";", { remap=true })
map({'n', 'x'}, "Æ", ":", { remap=true })
map({'n', 'x'}, "ø", "'", { remap=true })
map({'n', 'x'}, "Ø", '"', { remap=true })
map({'n', 'x'}, "å", "[", { remap=true })
map({'n', 'x'}, "Å", "{", { remap=true })
-- In Danglish I moved : and ; to the }] button
-- But this messes with things when I'm not in Danglish.

-- Mapping for function keys available on mechanical keyboard
-- TODO: too far away for being useful for something like regular completion, 
-- tab is more convenient. Find another use, e.g. snippet 
-- completion and navigation? Then you can free up ctrl+j/k/l although we're 
-- not always on mech keyboard, and I don't mind using <s-c-K> for the rare 
-- case of needing special chars.
map('i', "<F13>", "<C-j>", {remap=true})
map('i', "<F14>", "<C-k>", {remap=true})
map('i', "<F15>", "<C-l>", {remap=true})

-- When language is set to danish we can get ;:'" with alt the same way as we 
-- do in danglish keyboard input, however that is with the left alt which is 
-- deeper in the OS and the mapping below is for when esc+key (^[) is detected which 
-- is what is sent to the terminal, e.g. with kitty's setting 'macos_option_as_alt right'.
-- Note that they are also available with the ]} and \| keys.
map('i', "<A-;>", ";")
map('i', "<A-S-;>", ":")
map('i', "<A-'>", "'")
map('i', "<A-S-'>", '"')
-- And in case danglish keyboard is active:
map('i', "<A-æ>", ";")
map('i', "<A-S-æ>", ":")
map('i', "<A-ø>", "'")
map('i', "<A-S-ø>", '"')

-- Remap for dealing with word wrap
map('i', '<down>', [[v:count == 0 ? '<C-\><C-O>gj' : '<down>']], { expr=true, silent = true })
map('i', '<up>',   [[v:count == 0 ? '<C-\><C-O>gk' : '<up>']],   { expr=true, silent = true })
map('n', '<down>', [[v:count == 0 ? 'gj' : '<down>']],           { expr=true, silent = true })
map('n', '<up>',   [[v:count == 0 ? 'gk' : '<up>']],             { expr=true, silent = true })

-- Typos. I don't use command window much but I often press q: or q; when I mean :q
map('n', 'q:', ':q')
map('n', 'q;', ':q')
map('n', '<C-;>', 'q:')

-- Typos.
vim.api.nvim_create_user_command("Q", "q", {})
vim.api.nvim_create_user_command("X", "x", {})
vim.api.nvim_create_user_command("WQ", "wq", {})
vim.api.nvim_create_user_command("Wq", "wq", {})
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

local triples = { '"""', "'''", "```" }
local pairs   = { '""', "''", "``", "()", '[]', "{}", "<>", "$$" }
local singles = { "'", '"', '`', '(', ')', '[', ']', '{', '}', '<', '>', '$' }
-- hack map of shift+space
map({'i', 'n'}, "\x1F", function ()
    local r, c = unpack(vim.api.nvim_win_get_cursor(0))
    local line = vim.api.nvim_get_current_line()
    if vim.tbl_contains(triples, line:sub(c-2,c)) then
        return "<left><left><left>"
    elseif vim.tbl_contains(triples, line:sub(c+1,c+3)) then
        return "<right><right><right>"
    elseif c == #line then
        return "<left>"
    -- left priority over right for pairs so we go back first for e.g. Matrix{}|[]
    elseif vim.tbl_contains(pairs, line:sub(c-1,c)) then
        return "<left>"
    elseif vim.tbl_contains(pairs, line:sub(c+1,c+2)) then
        return "<right>"
    -- right priority over left for singles, since they are usually half of a 
    -- filled out pair and we want to prioritize progressing in that case.
    elseif vim.tbl_contains(singles, line:sub(c+1,c+1)) then
        return "<right>"
    elseif vim.tbl_contains(singles, line:sub(c,c)) then
        return "<left>"
    else
        return "<right>"
    end
end, {expr=true, desc="Move inside empty pair/triples or outside non-empty"})

-- useful with cursor | in {|} to get
-- {
--     |
-- }
vim.keymap.set('i', '<S-CR>', "<CR><Esc>O", { desc="Indented newline" })

-- switch to/from Danish æøå and to insert mode, which is convenient.
-- remap in order to utilise the remapped <C-^> which updates the cmp dictionary
map("n", "yod", "i<C-^>", { remap=true, desc="Danish (<C-^>)" })

-- small hack to remove excess whitespace possible since iw also captures 
-- whitespace under cursor.
map("n", "di ", "ciw <Esc>", { desc="Delete excess whitespace" })

-- like =p but for substitution
map("n", '=ss', 'ss=`]', {remap=true, silent=true, desc="Substitute+reindent"})
map("n", "<leader>Sr", "<cmd>source $XDG_CONFIG_HOME/nvim/after/plugin/luasnip.lua<CR>", {desc="Reload snippets"})

-- use the following two commands to enable spelling
-- setlocal spell
-- set spelllang=en_us
-- ctrl+s to fix last misspelled word
-- credit: https://castel.dev/post/lecture-notes-1/
-- with git: https://github.com/gillescastel/latex-snippets
-- "i_ctrl-g u" = something about undo history, probably to keep the edit part of the insert edit?
-- [s = goto last misspelled, 1z= = replace with first spell suggestion. 
map('i', '<C-s>', [[<c-g>u<Esc>[s1z=`]i<c-g>u]], { desc="Spell correct closest" })
map('n', '<C-s>', [=[[s1z=``]=], { desc="Spell correct closest" })

-- not the most elegant but it works.
-- LeftMouse to move cursor to pressed location.
-- Then set @/ to the current word (\< and \> are to search strictly).
-- Then enable hlsearch. This is all a way to search without going to the next 
-- match (if we just pressed * for instance)
map('n', '<RightMouse>',
    [[<LeftMouse>:let @/='\<'.expand('<cword>').'\>'|set hlsearch<CR>]],
    { silent=true, desc="Search pressed word" }
)

-- fallback search replace if both treesitter and LSP are not attached.
map('n', '<leader>rn', [[:%s/<C-r><C-w>/]], { desc="Search/replace cword" })
-- use a selection that isn't a perfect cword, or just to use the simple search/replace when LSP is attached etc.
map('x', '<leader>rn', [["ry:%s/<C-r>r/]], { desc="Search/replace" })

map('i', '<C-6>',
    [[<C-^><C-\><C-o>:doautocmd User ToggleDansk<CR>]],
    { silent=true, desc="Toggle dansk" }
)
-- Ins key was a bit useless just doing what i does so let's make it a language switch insertion:
map('n', '<Ins>', 'i<C-6>', { silent=true, desc="Insert + Toggle dansk" })
-- <sa = my keybind for enable setting autoformat
-- gww = autoformat line
-- 0   = goto column 0, so we scroll all the way back to the right
-- gi  = go to last insert location and enter insert mode. Works even with the change to the line.
map('i', '<C-S-A>', "<Esc><sagww0gi", { remap=true, desc="Enable autoformat and apply it." })

map('n', "<leader>qo", "<Cmd>copen<CR>", {desc="Open"})
map('n', "<leader>qq", "<Cmd>cclose<CR>", {desc="Close"}) -- q for quit and is fast
map('n', "<leader>q1", "<Cmd>cc 1<CR>", {desc="Entry 1"})
map('n', "<leader>q2", "<Cmd>cc 2<CR>", {desc="Entry 2"})
map('n', "<leader>q3", "<Cmd>cc 3<CR>", {desc="Entry 3"})
-- we don't map :cnext etc here since we have ]q etc

map('n', "<leader>Qo", "<Cmd>lopen<CR>", {desc="Open"})
map('n', "<leader>QQ", "<Cmd>lclose<CR>", {desc="Close"}) -- Q for quit and is fast
map('n', "<leader>Q1", "<Cmd>ll 1<CR>", {desc="Entry 1"})
map('n', "<leader>Q2", "<Cmd>ll 2<CR>", {desc="Entry 2"})
map('n', "<leader>Q3", "<Cmd>ll 3<CR>", {desc="Entry 3"})
-- we don't map :lnext etc here since we have ]l etc

map('n', '<leader>bd', "<Cmd>bd<CR>", { desc="Delete" })
map('n', '<leader>bD', "<Cmd>bd!<CR>", { desc="Delete!" })
map('n', "<leader>bn", "<Cmd>enew<CR>", { desc="New" })
map('n', "<leader>bc", "<Cmd>tabclose<CR>", { desc="tabclose" })

-- poor fix for cmp replacing capitlisation of buffer words.
-- default C-c is cancel, so exit to normal mode. Not very useful.
map('i', "<C-c>", "<Esc>b~gi", { desc="Capitalise last word" })

