#!/usr/bin/env lua
local cmd = vim.cmd
local api = vim.api
local fn = vim.fn

vim.keymap.set("n", "tsm", "<plug>(vimtex-env-toggle-math)")
vim.keymap.set("n", "dsm", "<plug>(vimtex-env-delete-math)")
vim.keymap.set("n", "csm", "<plug>(vimtex-env-change-math)")
vim.keymap.set("n", "xad", "yaddad", {remap=true})
vim.keymap.set("n", "xid", "yiddid", {remap=true})
-- item with i instead of m and math with m
vim.keymap.set({"o", "x"}, "ai", "<Plug>(vimtex-am)")
vim.keymap.set({"o", "x"}, "ii", "<Plug>(vimtex-im)")
vim.keymap.set({"o", "x"}, "am", "<Plug>(vimtex-a$)")
vim.keymap.set({"o", "x"}, "im", "<Plug>(vimtex-i$)")
-- shorthand to $ just using 4 ($ without shift)
vim.keymap.set({"o", "x"}, "a4", "<Plug>(vimtex-a$)")
vim.keymap.set({"o", "x"}, "i4", "<Plug>(vimtex-i$)")

-- async overleaf sync on save and load
local overleaf = api.nvim_create_augroup("overleaf", {clear=true})
local remote = fn.system("git remote get-url origin")
if string.sub(remote, 1, 25) == "https://git.overleaf.com/" then
    print("Adding overleaf aucmd")
    api.nvim_create_autocmd({"BufRead", "BufWritePost"}, {
        group = overleaf,
        buffer = 0, -- only for current buffer. Mutually exclusive with pattern arg.
        callback = function ()
            fn.jobstart({fn.expand("$XDG_CONFIG_HOME/nvim/latex/overleaf/gitsync.sh"), vim.api.nvim_buf_get_name(0)})
        end
    })
end


vim.keymap.set("n", "<plug>alignTable", "gaie<C-f>v/\\multicol<CR>*&", {silent=true, remap=true})

local function in_env(name)
    local is_inside = vim.fn['vimtex#env#is_inside'](name)
    return (is_inside[1] > 0 and is_inside[2] > 0)
end

-- return tex table column 0-index by counting (unescaped) ampersands (&) before cursor "column" c.
local function getCurrentColumn(c)
    local line = vim.api.nvim_get_current_line()
    return #line:sub(1, c):gsub('\\&', ''):gsub('[^&]', '')
end
function texCol()
    return getCurrentColumn(vim.api.nvim_win_get_cursor(0)[2])
end

-- inspired by https://gist.github.com/tpope/287147
vim.keymap.set("i", "&", function()
    local r, c = unpack(vim.api.nvim_win_get_cursor(0))
    
    -- -1 since nvim_win_get_cursor is (1,0)-indexed and nvim_buf_set_text is 0-indexed.
    if not in_env("table") or vim.api.nvim_buf_get_text(0, r-1, c-1, r-1, c, {})[1] == "\\" then
        -- just write a regular & and exit
        vim.api.nvim_buf_set_text(0, r-1, c, r-1, c, {'&'})
        vim.api.nvim_win_set_cursor(0, {r, c+1})
        return
    end

    vim.api.nvim_buf_set_text(0, r-1, c, r-1, c, {'& '})
    local column = getCurrentColumn(c + 2)
    
    cmd [[silent! exe "normal \<Plug>alignTable"]]
    
    local line = vim.fn.getline(r)
    -- replace escaped ampersands with something else so they don't get counted
    _, c = line:gsub('\\&', '  '):find(('[^&]*&'):rep(column))
    vim.api.nvim_win_set_cursor(0, {r, c+1})
    -- Hacks to clear message area.
    -- Needed in combination with silent! above to not see any prints.
    print " "
    cmd "echo ' '"
end, {remap=false, silent=true, buffer=true})





