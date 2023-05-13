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
    local line = vim.api.nvim_get_current_line():gsub('\\&', '  ')
    local nAmp = #line:sub(1,c):gsub('[^&]', '')
    if not line:match('multicol') then return nAmp else
        local beforeCell = line:sub(1,c):match('.*&') -- matches the last & before cursor
        if beforeCell then
            for n in beforeCell:gmatch[[\multicolumn{(%d+)}]] do
                nAmp=nAmp + n-1
            end
        end
        return nAmp
    end
end
-- return tex table column 0-index assuming the cursor (with "column" c) is at the preamble line.
local function getCurrentColumnPre(c)
    local line = vim.api.nvim_get_current_line()
    -- find preamble range using balanced bracket patterns.
    -- Empty () to get location (instead of an extracted substring).
    -- \begin{tabularx}{\textwidth}{@{}lcX@{}}
    local pre_begin, pre_end = line:match('\\begin{tabularx}%b{}()%b{}()')
    pre_begin=pre_begin+1
    pre_end=pre_end-2
    return math.max(0, #line:sub(pre_begin, c+1):gsub('[^%a]', '')-1)
end
-- jump to or from the relevant column in the preamble for a tabularx table.
function gotoPreTable()
    local r_tabularx, c_tabularx = unpack(vim.fn['vimtex#env#is_inside']("tabularx"))
    if r_tabularx == 0 then return end
    
    local r, c = unpack(vim.api.nvim_win_get_cursor(0))
    if r == r_tabularx then
        local col = getCurrentColumnPre(c)
        -- goto first line with the most cells found.
        -- scan 9 lines for now. 0-indexing means r excludes preamble.
        local scan = vim.api.nvim_buf_get_lines(0, r, r+10, false)
        local iGoto, maxAmp = r, 0
        for i, v in ipairs(scan) do
            scan[i] = v:gsub('\\&', '  ')
            _, nAmp = scan[i]:gsub('&', '')
            if nAmp > maxAmp then
                iGoto, maxAmp = i, nAmp
            end
        end
        local lGoto = scan[iGoto]
        local _, cGoto = lGoto:find(("[^&]*&"):rep(col))
        if cGoto == nil then
            -- there was no match so the line is lacking &s
            -- we go to EOL although before optional \\
            if lGoto:sub(-2) == [[\\]] then cGoto = #lGoto-4
            else cGoto = #lGoto end
        end
        -- go to first non-whitespace as long as not & (indicating the next cell)
        local _, firstNonW = lGoto:sub(cGoto+1):find("^%s*[^%s&]")
        cGoto = cGoto + (firstNonW or 1)-1
        -- Mark current location for jumplist to work (e.g. <c-o> to go back in table)
        -- lua version doesn't work
        cmd "normal m`"
        vim.api.nvim_win_set_cursor(0, {r+iGoto, cGoto})
    else
        local col = getCurrentColumn(c)
        -- match col number of alphanumerics
        -- -1 due to zero-indexing.
        local linePre = vim.api.nvim_buf_get_lines(0, r_tabularx-1, r_tabularx, true)[1]
        -- +1 and -1 so the 0-th column will lead to the first alphanumeric 
        -- char, rather than the @{} that might precede
        local _, cGoto = linePre:find("\\begin{tabularx}%b{}{" .. ("[^%a]*%a"):rep(col+1))
        -- Mark current location for jumplist to work (e.g. <c-o> to go back in table)
        -- lua version doesn't work
        cmd "normal m`"
        vim.api.nvim_win_set_cursor(0, {r_tabularx, cGoto-1})
    end
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
    local line = vim.api.nvim_get_current_line()
    local column = #line:sub(1, c + 2):gsub('\\&', ''):gsub('[^&]', '')
    
    cmd [[silent! exe "normal \<Plug>alignTable"]]
    
    -- replace escaped ampersands with something else so they don't get counted
    _, c = line:gsub('\\&', '  '):find(('[^&]*&'):rep(column))
    vim.api.nvim_win_set_cursor(0, {r, c+1})
    -- Hacks to clear message area.
    -- Needed in combination with silent! above to not see any prints.
    print " "
    cmd "echo ' '"
end, {remap=false, silent=true, buffer=true})





