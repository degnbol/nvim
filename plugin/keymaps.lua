local util = require "utils/init"
local map = require "utils/keymap"
local ts = require "utils/treesitter"

require "keymaps/options"
require "keymaps/danglish"
require "keymaps/surround"
require "keymaps/blockim"
require "keymaps/comments"

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
vim.api.nvim_create_user_command("X", "x", {})
vim.api.nvim_create_user_command("WQ", "wq", {})
vim.api.nvim_create_user_command("Wq", "wq", {})
vim.api.nvim_create_user_command("Lw", "w", {})
-- abbrev instead of command since command has to start with uppercase
vim.cmd [[cnoreabbrev qq q]]

-- small hack to remove excess whitespace possible since iw also captures
-- whitespace under cursor.
map.n("di ", "ciw <Esc>", "Delete excess whitespace")

-- like =p but for substitution
map.n('=ss', 'ss=`]', "Substitute+reindent", { remap = true, silent = true })
map.n("<leader>Sr", "<cmd>source $XDG_CONFIG_HOME/nvim/after/plugin/luasnip.lua<CR>", "Reload snippets")

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

-- Disabled since other things might need to override ESC.
-- map.n('<Esc>', ":noh<CR><Esc>", "Disable search highlights", { silent = true, remap = false })

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

map.n("<leader>Qo", "<Cmd>lopen<CR>", "Open")
map.n("<leader>QQ", "<Cmd>lclose<CR>", "Close") -- Q for quit and is fast
map.n("<leader>Q1", "<Cmd>ll 1<CR>", "Entry 1")
map.n("<leader>Q2", "<Cmd>ll 2<CR>", "Entry 2")
map.n("<leader>Q3", "<Cmd>ll 3<CR>", "Entry 3")
-- we don't map :lnext etc here since we have ]l etc

map.n('<leader>bd', function () require"mini.bufremove".delete(vim.v.count) end, "bdel without win close")
map.n('<leader>bD', function() require"mini.bufremove".delete(vim.v.count, true) end, "bdel! without win close")
map.n('<leader>bw', function () require"mini.bufremove".wipeout(vim.v.count) end, "bwipe without win close")
map.n('<leader>bW', function() require"mini.bufremove".wipeout(vim.v.count, true) end, "bwipe! without win close")
map.n('<leader>bu', function () require"mini.bufremove".unshow(vim.v.count) end, "unshow without win close")
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

map.c('<A-left>', "<s-left>", "move back one word")
map.c('<A-right>', "<s-right>", "move forward one word")
map.c('<A-BS>', "<C-w>", "Delete back one word")
-- Doesn't work well since we go one WORD to the right but only delete one word back.
map.c('<A-delete>', "<S-right><C-w>", "Delete next word")
-- this ignored if kitty handles it.
map.c('<D-BS>', "<C-u>", "Delete to beginning of line")

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
map.n('<leader>th', "<Cmd>TSBufToggle highlight<CR>", "Toggle local highlight")
map.n('<leader>tH', "<Cmd>TSToggle highlight<CR>", "Toggle global highlight")
map.n('<leader>ti', vim.show_pos, "Inspect") -- Same as :Inspect
map.n('<leader>tt', "<Cmd>InspectTree<CR>", "Inspect tree")
map.n('<leader>tI', "<Cmd>Capture TSInstallInfo<CR>", "Install info")
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
map.n("<leader>li", "<Cmd>LspInfo<CR>", "Info")
-- can't use backspace since it is hardcoded by mini.clue for up one level
map.n("<leader>l0", "<Cmd>LspStop<CR>", "Stop")
map.n("<leader>l1", "<Cmd>LspStart<CR>", "Start")
map.n("<leader>l!", "<Cmd>LspRestart<CR>", "Restart")
map.n("<leader>lL", "<Cmd>LspLog<CR>", "Log")
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

-- Goto references excluding the line it's called from.
map.n('grr', function()
    vim.lsp.buf.references(nil, map.filter_lsp_items(function(item)
        return not map.qf_item_is_self(item)
    end))
end, "Filtered references")
map.desc('n', 'gra', "Code actions")
map.desc('n', 'gri', "Implementations")
map.desc('n', 'grn', "Rename")
map.desc('n', 'grt', "Type definitions")

-- TODO: the LSP nmap K should also show the float aligned to the function name start and not cursor.
map.i('<C-s>', function()
    -- TODO: modify float to remove empty lines at top and bottom.
    -- TODO: update signature help when pressing comma or deleting a comma or moving cursor.
    -- Also decide if repeated <C-s> should cycle the signatures, as is default.
    local call_expression = ts.get_parent('call_expression')
    if call_expression == nil then return end
    local r, c = unpack(vim.api.nvim_win_get_cursor(0))
    local start_row, start_col, end_row, end_col = vim.treesitter.get_node_range(call_expression)
    vim.lsp.buf.signature_help { title = nil, offset_x = start_col - c, close_events = { "WinScrolled", "ModeChanged" } }
end, "Signature help")

vim.api.nvim_create_autocmd("LspAttach", {
    group = vim.api.nvim_create_augroup("my_lsp_attach", { clear = true }),
    callback = function()
        -- TODO: complete this after work
        -- map.n("K", function()
        --     local r, c = unpack(vim.api.nvim_win_get_cursor(0))
        --     local cword_start, _ = util.cword_cols()
        --     vim.lsp.buf.hover { offset_x = cword_start - c, }
        -- end, "LSP hover", { buffer = true })
    end
})

-- custom gx function that supports more website links.
map.n("gx", function()
    -- Go to github for plugins easily.
    if vim.bo.filetype == "lua" then
        -- First check if we are editing a file read by lazy.nvim
        local rtp = vim.opt.runtimepath:get()[1]
        local filepath = vim.api.nvim_buf_get_name(0)
        if filepath:match(rtp .. '/lua/plugins/') then
            -- get first string on the line, assumes we don't list multiple plugins on one line.
            local line = vim.api.nvim_get_current_line()
            local repo = line:match([["([%w%p]+/[%w%p]+)"]])
            repo = repo or line:match([['([%w%p]+/[%w%p]+)']])
            if repo then
                if repo:match("http") then
                    return vim.ui.open(repo)
                else
                    return vim.ui.open("https://github.com/" .. repo)
                end
            end
        end
    end

    -- for latex packages...
    if vim.bo.filetype == "tex" then
        local line = vim.api.nvim_get_current_line()
        local pac = line:match("\\usepackage.*{([%w_-]+)}")
        if pac ~= nil then
            local ctan = "https://ctan.org/pkg/" .. pac .. "?lang=en"
            return vim.ui.open(ctan)
        end
    end

    -- Using various-textobjs url finder, we will detect url further ahead,
    -- e.g. useful when being lazy and the url is right there on the line.
    -- visually select URL
    require("various-textobjs").url()
    -- plugin only switches to visual mode when textobj found
    local foundURL = vim.fn.mode():find("v")
    -- retrieve URL with the z-register as intermediary
    vim.cmd.normal { '"zy', bang = true }
    local url = vim.fn.getreg("z")
    vim.ui.open(url)
end, "Smart URL opener")
