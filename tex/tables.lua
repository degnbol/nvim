#!/usr/bin/env lua
local vtu = require "vimtex_util"
local is_inside = vtu.is_inside
local in_env = vtu.in_env
local cmd = vim.cmd


vim.keymap.set("n", "<plug>AlignTable", "gaie<C-f>v/\\multicol<CR>*&", {silent=true, remap=true})

-- get tabular-like env start r and c
local function get_tabular()
    local r, c = unpack(is_inside("tabularx"))
    if r > 0 or c > 0 then return r, c end
    r, c = unpack(is_inside("tabulary"))
    if r > 0 or c > 0 then return r, c end
    r, c = unpack(is_inside("tabular"))
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
vim.keymap.set("n", "<plug>TexJumpPre", function()
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
end)

-- TODO: let multicolumn align to the rest but not the other way. What to do if not one row per line?
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

-- TODO: deal with \& and maybe add more tabs after multicolumn (which we will have to remove later in the inverse function.)
-- also, is the \\ missing after toprule etc? if so, edit the snippet template.
vim.keymap.set("n", "<plug>YankTable", function()
    local delim = '\t'
    cmd [[normal "tyie]]
    local content = vim.fn.getreg("t")
    content = content:
    	gsub('[\n\r]', ''):
    	gsub('^%s*', ''):
    	gsub('%s*$', ''):
    	gsub([[%s*\\%s*]], '\n'):
    	gsub('[ \t]*&[ \t]*', delim) -- can't replace %s since it also removes \n
    -- set system clipboard
    vim.fn.setreg('+', content)
end)

vim.keymap.set("n", "<plug>PasteTable", function()
    local delim = '\t'
    -- TODO inverse of above
end)



