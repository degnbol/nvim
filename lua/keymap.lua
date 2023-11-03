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
map('n', "U", "<C-r>")

-- For remote with problematic clipboard we replace the cutlass x that works 
-- locally with a mapping where we call y (which is mapped to smartyank) then 
-- d.
map('n', "xx", "yydd", {remap=true})
map('x', "x",  "ygvd", {remap=true})
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
map('n', 'j', "v:count == 0 ? 'gj' : 'j'", { expr = true, silent = true })
map('n', 'k', "v:count == 0 ? 'gk' : 'k'", { expr = true, silent = true })
map('i', '<down>', [[v:count == 0 ? '<C-\><C-O>gj' : '<down>']], { expr = true, silent = true })
map('i', '<up>',   [[v:count == 0 ? '<C-\><C-O>gk' : '<up>']],   { expr = true, silent = true })

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

-- for when you are at "...|" and want to exit the quotes and your fingers 
-- naturally find the quote key
map('i', "<C-'>", "<right>")
map('i', "<C-0>", "<right>")

-- switch to/from Danish æøå and to insert mode, which is convenient.
-- remap in order to utilise the remapped <C-^> which updates the cmp dictionary
map("n", "yod", "i<C-^>", { remap=true, desc="Danish (<C-^>)" })

-- small hack to remove excess whitespace possible since iw also captures 
-- whitespace under cursor.
map("n", "di ", "ciw <Esc>", { desc="Delete excess whitespace" })

-- like =p but for substitution
map("n", '=ss', 'ss=`]', {remap=true, silent=true, desc="Substitute+reindent"})
map("n", "<leader>Sr", "<cmd>source $XDG_CONFIG_HOME/nvim/after/plugin/luasnip.lua<CR>", {desc="Reload snippets"})

-- <leader>b can be for buffer related things, (bufdel, but :bd is already 
-- fast) or bibliography things e.g. with papis but this is not added yet and 
-- could also just fit under latex (<leader>l) since I don't really use it in 
-- markdown.

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

map('n', '<leader>bd', ":bd<CR>", { desc="Delete" })
map('n', '<leader>bD', ":bd!<CR>", { desc="Delete!" })

