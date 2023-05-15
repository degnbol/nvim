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
-- get tabular-like env start r and c
local function get_tabular()
    local r, c = unpack(vim.fn['vimtex#env#is_inside']("tabularx"))
    if r > 0 or c > 0 then return r, c end
    r, c = unpack(vim.fn['vimtex#env#is_inside']("tabulary"))
    if r > 0 or c > 0 then return r, c end
    r, c = unpack(vim.fn['vimtex#env#is_inside']("tabular"))
    if r > 0 or c > 0 then return r, c end
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

-- get r, c, contents of preamble for a containing tabular/tabularx/tabulary env.
local function getPreTabular()
    local r_tabular, _ = get_tabular()
    if r_tabular == 0 then return end
    local line = vim.fn.getline(r_tabular)
    -- find preamble range using balanced bracket patterns.
    -- Empty () to get location (instead of an extracted substring).
    -- example to match:
    -- \begin{tabular}[c]{@{}lcX@{}}
    -- \begin{tabularx}{\textwidth}{@{}lcX@{}}
    -- \begin{tabulary}{1.0\textwidth}{L|L|L|L|L|L}
    local pre_begin, pre_end = line:match('\\begin{tabular.}%b{}()%b{}()')
    if pre_begin == nil then -- tabular
        pre_begin, pre_end = line:match('\\begin{tabular}%b[]()%b{}()')
    end -- cannot test pre_begin from first if-statement in second if it is an elseif
    if pre_begin == nil then
        pre_begin, pre_end = line:match('\\begin{tabular}()%b{}()')
    end
    pre_begin=pre_begin+1
    pre_end=pre_end-2
    return r_tabular, pre_begin, line:sub(pre_begin, pre_end)
end

-- return tex table column 0-index assuming the cursor (with "column" c) is at the preamble line.
local function getCurrentColumnPre(c, c_pre, pre)
    return math.max(0, #pre:sub(1, c+1-c_pre):gsub('[^%a]', ''))
end
-- jump to or from the relevant column in the preamble for a tabular/tabularx/tabulary table.
local function gotoPreTable()
    local r_tabular, c_pre, pre = getPreTabular()
    if r_tabular == nil then return end
    
    local r, c = unpack(vim.api.nvim_win_get_cursor(0))
    if r == r_tabular then
        local col = getCurrentColumnPre(c, c_pre, pre)
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
            if lGoto:sub(-2) == [[\\]] then cGoto = #lGoto-2
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
        -- char, rather than the @{} that might precede
        local _, cGoto = pre:find(("[^%a]*%a"):rep(col+1))
        if cGoto == nil then
            -- the preamble is lacking entries so we go to its end.
            cGoto = c_pre + #pre
        else
            cGoto=cGoto+c_pre-2
        end
        -- Mark current location for jumplist to work (e.g. <c-o> to go back in table)
        -- lua version doesn't work
        cmd "normal m`"
        vim.api.nvim_win_set_cursor(0, {r_tabular, cGoto})
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
    
    -- the alignTable call moves the cursor and modifies lines.
    -- The cursor is now at the top of inside the table env.
    line = vim.fn.getline(r)
    -- replace escaped ampersands with something else so they don't get counted
    _, c = line:gsub('\\&', '  '):find(('[^&]*&'):rep(column))
    vim.api.nvim_win_set_cursor(0, {r, c+1})
    -- Hacks to clear message area.
    -- Needed in combination with silent! above to not see any prints.
    print " "
    cmd "echo ' '"
end, {remap=false, silent=true, buffer=true})





