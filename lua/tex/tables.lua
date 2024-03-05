#!/usr/bin/env lua
local util = require "utils/init"
local vtu = require "utils/vimtex"
local is_inside = vtu.is_inside
local in_env = vtu.in_env
local cmd = vim.cmd

local M = {}

-- TODO:
-- functions or ways of moving/deleting/adding rows/columns.
-- Maybe each a function but ideally a more modular (vim) way, e.g. selection of row/col,
-- moving between cells, etc. But might not be possible without elaborate 
-- function due to things like header row, and not having one row per line.

-- TODO: can we get sideways scroll to lag less

-- TODO: maybe conceal makecell, and somehow conceal overflowing cell or 
-- consider conceal setting when deciding space a cell takes up (after conceal)

-- get tabular-like env start r and c
local function getTabular()
    local r, c = unpack(is_inside("tabularx"))
    if r > 0 or c > 0 then return r, c end
    r, c = unpack(is_inside("tabulary"))
    if r > 0 or c > 0 then return r, c end
    r, c = unpack(is_inside("tabular"))
    if r > 0 or c > 0 then return r, c end
end

-- return tex table column 0-index by counting (unescaped) ampersands (&) before cursor "column" c.
-- Deprecate, since it doesn't read multi. Use getCurrentCell instead with 
-- parseTable
local function getCurrentColumn()
    local _, c = unpack(vim.api.nvim_win_get_cursor(0))
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

-- get {r, c, string} of preamble for a containing tabular/tabularx/tabulary env.
-- where string is e.g. "@{}p{1.1cm}rrlp{1.1cm}rrrrlX@{}"
-- TODO: handle multiline
local function getPreTabular()
    local r_tabular, _ = getTabular()
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

--- get {r, c, string, string[], int[]}, where the latter is an array of each element in the preamble, e.g.
--- { "@{}", "p{1.1cm}", "r", "|", "r", "||", " ", "p{1.1cm}", "r", "r", ">{l}", "r", "l", "X", "@{}" }
local function parsePre()
    local r, c, pre = getPreTabular()
    local parsed = {}
    local _pre = pre
    while #_pre > 0 do
        local part = _pre:match("^.%b{}")
        if part ~= nil then
            table.insert(parsed, part)
            _pre = _pre:sub(1+#part)
        else
            -- for e.g. ||
            local nonalpha = _pre:match("^[^%a]+")
            if nonalpha ~= nil then
                table.insert(parsed, nonalpha)
                _pre = _pre:sub(1+#nonalpha)
            else
                table.insert(parsed, _pre:sub(1,1))
                _pre = _pre:sub(2)
            end
        end
    end
    -- attach modifiers <{...} and >{...} to what they modify
    local parsedMod = {}
    local i=0
    while i<#parsed do
        i=i+1
        if parsed[i]:match("^[><]") then
            table.insert(parsedMod, parsed[i]..parsed[i+1])
            i=i+1
        else
            table.insert(parsedMod, parsed[i])
        end
    end
    -- create column mapping to parsed parts
    local columns = {}
    for i, part in ipairs(parsedMod) do
        if part:match("^[><%a]") then
            table.insert(columns, i)
        end
    end
    return r, c, pre, parsedMod, columns
end

-- return tex table column 0-index assuming the cursor (with "column" c) is at the preamble line.
local function getCurrentColumnPre(c, c_pre, pre)
    local until_cursor = pre:sub(1, c+2-c_pre)
    local nAlpha = #until_cursor:gsub("%b{}", ""):gsub('[^%a]', '')
    -- 0-indexing
    return math.max(0, nAlpha-1)
end

-- jump to or from the relevant column in the preamble for a tabular/tabularx/tabulary table.
vim.keymap.set("n", "<plug>TableJumpPre", function()
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
        if firstNonW == nil then -- the cell is empty
            if cGoto == 0 then -- the very first cell
                -- we want to go to indentation.
                -- Take indentation from previous line.
                -- ERRORs if we are jumping to first line in table.
                -- Maybe rethink the scan 10 lines logic to jump to first line with >1 cell (jump over hline and without multicol)
                firstNonW = scan[iGoto-1]:find("%S")
            else
                firstNonW = 1
            end
        end
        cGoto = cGoto + firstNonW-1
        -- Mark current location for jumplist to work (e.g. <c-o> to go back in table)
        -- lua version doesn't work
        cmd "normal m`"
        vim.api.nvim_win_set_cursor(0, {r+iGoto, cGoto})
    else
        local col = getCurrentColumn()
        -- go across preamble string for "col" alphanumerics ignoring {...}
        local bracketDepth = 0
        local cGoto = 0
        while cGoto < #pre do
            cGoto=cGoto+1
            local char = pre:sub(cGoto, cGoto)
            if char == '{' then bracketDepth=bracketDepth+1
            elseif char == '}' then bracketDepth=bracketDepth-1
            elseif bracketDepth == 0 and char:match("%a") then
                col=col-1
                if col < 0 then break end
            end
        end
        cGoto=cGoto+c_pre-2
        -- Mark current location for jumplist to work (e.g. <c-o> to go back in table)
        -- lua version doesn't work
        cmd "normal m`"
        vim.api.nvim_win_set_cursor(0, {r_tabular, cGoto})
    end
end, {desc="Goto preamble"})



function parseTable()
    local envname, r1, _, r2, _ = unpack(vtu.get_env())
    if not envname:match("^tabular") then return end
    -- -1 since we are using 1-indexed r1, r2 in 0-indexed end-exclusive function
    -- to get the lines withinh the end, so excluding the begin/end lines
    r2=r2-1
    local lines = vim.api.nvim_buf_get_lines(0, r1, r2, true)
    -- final newline necessary to not lose final cell
    local text = table.concat(lines, '\n') .. '\n'
    local indentation = text:match('^[%s\t]*')
    local texts = {}
    local rows = {} -- row in table
    local cols = {} -- column in table
    local presufs = {} -- column prefix (if whitespace), suffix (if ends in newline) or empty string.
    local row = 0
    local col = 0
    local lasti = 0
    local presuf = indentation
    -- rs and cs marks the first char of a cell, e.g. right after & or screen col 0.
    local rs = {} -- 1-indexed
    local cs = {} -- 0-indexed
    local r = r1+1 -- +1 to move past preamble line. TODO: this will break with multiline
    local c = 0 -- keeps track of the beginning of the section in screen columns
    local maxrow = 0 -- max among cells, i.e. excluding \bottomrule

    -- find the position of each & and \\, the latter assumed to be at end of line.
    -- This assumption avoids \\ inside cell with forced linebreak.
    for i, m in text:gmatch('()([&\n])') do
        local section = text:sub(lasti+1, i-1)
        if m == '&' and text:sub(i,i) ~= '\\' then
            -- the minus works like *, except match as shortly as possible, which let's us trim whitespace.
            local text = section:match('\n?%s*([^\n]-)%s*$')
            table.insert(texts, text)
            table.insert(rows, row)
            table.insert(cols, col)
            table.insert(presufs, presuf)
            table.insert(rs, r)
            table.insert(cs, c)
            local colspan = text:match("^\\multicolumn{(%d+)}") or 1
            col=col + colspan
            presuf = ""
            c=c+i-lasti
        elseif m == '\n' then
            local cell = section:match('\n?%s*([^\n]-)%s*\\\\%s*$')
            -- doesn't have \\, if
            -- unfinished cell,
            -- or a line rule (\toprule etc)
            if cell == nil then
                table.insert(texts, section)
                table.insert(rows, row)
                table.insert(cols, -1) -- mark to ignore
                table.insert(presufs, '\n')
                table.insert(rs, r)
                table.insert(cs, c)
            else
                table.insert(texts, cell)
                table.insert(rows, row)
                table.insert(cols, col)
                table.insert(presufs, ' \\\\\n')
                table.insert(rs, r)
                table.insert(cs, c)
                maxrow=row
                row=row+1
                col=0
            end
            r=r+1
            c = 0
            presuf = indentation
        end
        lasti = i
    end
    local maxcol = math.max(unpack(cols))
    return {
        texts=texts,
        rows=rows,
        cols=cols,
        presufs=presufs,
        maxrow=maxrow,
        maxcol=maxcol,
        r1=r1,
        r2=r2,
        rs=rs,
        cs=cs,
    }
end

---tab: table {r1=int, r2=int, maxcol=int, texts=string[], cols=int[], presufs=string[] }
---maxwidth: maximum allowed column width
---usemax: in case a cell exceeds the max, should it be ignored (default), or 
---should the column be set to maxwidth (usemax=true).
local function writeTable(tab, opts)
    opts = opts or {}
    maxwidth = opts.maxwidth or 20
    usemax = opts.usemax or false
    -- column width preallocate zeros
    local widths = {}
    for i = 0, tab.maxcol do widths[i] = 0 end
    for i, col in ipairs(tab.cols) do
        if col >= 0 then
            local cell = tab.texts[i]
            -- ignore multicolumn
            if not cell:match("^\\multicolumn") then
                -- get column width, restricted by a max allowed,
                if #cell < maxwidth then
                    widths[col] = math.max(widths[col], #cell)
                elseif usemax then
                    widths[col] = maxwidth
                end
            end
        end
    end
    -- build new text
    local text = {}
    local idealw = 0
    local currentw = 0
    for i, cell in ipairs(tab.texts) do
        local col = tab.cols[i]
        local presuf = tab.presufs[i]
        -- flagged for ignoring, which is final cell without & or \\ or e.g. 
        -- \toprule
        if col == -1 then
            table.insert(text, cell .. presuf)
            idealw, currentw = 0, 0
        else
            idealw=idealw + widths[col]
            local multicol = cell:match("^\\multicolumn{(%d+)}")
            if multicol ~= nil then
                for i = 1, multicol-1 do
                    idealw=idealw + 3 + widths[col+i]
                end
            end
            if presuf == "" or presuf:match('\n$') then
                if cell ~= "" then cell = " " .. cell end
            else -- is indentation
                cell = presuf..cell
                idealw=idealw+#presuf
            end
            currentw=currentw + #cell
            if currentw < idealw then
                cell = cell .. string.rep(' ', idealw - currentw)
                currentw = idealw
            end
            if presuf:match("\n$") then
                cell = cell .. presuf
                idealw, currentw = 0, 0
            else
                if cell ~= "" or currentw+1 < idealw then
                    cell=cell.." "
                    currentw=currentw+1
                end
                cell=cell.."&"
                currentw=currentw+1
                idealw=idealw+3
            end
            table.insert(text, cell)
        end
    end
    local lines = vim.split(table.concat(text):gsub('\n$', ''), '\n', {plain=true})
    vim.api.nvim_buf_set_lines(0, tab.r1, tab.r2, true, lines)
end


--- Work for either a visual selection or for the current env.
---maxwidth: maximum allowed column width
---usemax: in case a cell exceeds the max, should it be ignored (default), or 
---should the column be set to maxwidth (usemax=true).
local function alignTable(opts)
    if not in_env("table") then return end
    writeTable(parseTable(), opts)
end

---return row, col. 0-indexed
local function getCurrentCell(tab)
    -- both this and tab.rs, tab.cs are (1,0)-indexed
    local rCur, cCur = unpack(vim.api.nvim_win_get_cursor(0))
    for i = 1, #tab.texts do
        local r = tab.rs[i]
        local c = tab.cs[i]
        -- check that we have gone past a cell
        if r > rCur or r == rCur and c > cCur then
            return tab.rows[i-1], tab.cols[i-1]
        end
    end
    -- last entry
    return tab.rows[#tab.rows], tab.cols[#tab.cols]
end

-- Move cursor to cell at given 0-indexed table row, col.
-- Moves to first non-whitespace, if any
local function gotoCell(tab, row, col)
    row = math.max(row, 0)
    col = math.max(col, 0)
    row = math.min(row, tab.maxrow)
    col = math.min(col, tab.maxcol)
    for i = 1, #tab.texts do
        if row == tab.rows[i] and col == tab.cols[i] then
            local r = tab.rs[i]
            local c = tab.cs[i]
            if tab.texts[i] ~= "" then
                local presuf = tab.presufs[i]
                if presuf == "" or presuf:match("\n$") then
                    c=c+1
                else -- indentation
                    c=c+#presuf
                end
            end
            vim.api.nvim_win_set_cursor(0, {r, c})
            return true
        end
    end
end

-- NOTE: what about when conceal is on an \textbf makes a cell look shorter than it is. Mostly relevant for header row.
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
    
    cmd [[silent! exe "normal \<Plug>TableAlign"]]
    
    -- the alignTable call moves the cursor and modifies lines.
    -- The cursor is now at the top of inside the table env.
    line = vim.fn.getline(r)
    -- replace escaped ampersands with something else so they don't get counted
    _, c = line:gsub('\\&', '  '):find(('[^&]*&'):rep(column))
    vim.api.nvim_win_set_cursor(0, {r, c+1})
    -- Hacks to clear message area.
    -- Needed in combination with silent! above to not see any prints.
    -- Another solution I saw somewhere temporarily redefines some print functions to not do anything.
    print " "
    cmd "echo ' '"
end, {remap=false, silent=true, buffer=true})

local function deleteColumn(opts)
    if not in_env("table") then return end
    local tab = parseTable()

    -- 0-indexed
    local _, colDel = getCurrentCell(tab)
    local texts, cols, presufs = {}, {}, {}
    for i, text in ipairs(tab.texts) do
        local col = tab.cols[i]
        local presuf = tab.presufs[i]
        if col < colDel then
            table.insert(texts, text)
            table.insert(cols,  col)
            table.insert(presufs, presuf)
        elseif col > colDel then
            table.insert(texts, text)
            table.insert(cols,  col-1)
            -- keep indent
            if col == 1 then
                table.insert(presufs, tab.presufs[i-1])
            else
                table.insert(presufs, presuf)
            end
        end
    end

    tab.texts = texts
    tab.cols = cols
    tab.presufs = presufs
    tab.maxcol=tab.maxcol-1
    writeTable(tab, opts)

    -- remove entry from preamble
    local r, c, pre, parts, premap = parsePre()
    -- +1 for index conv
    table.remove(parts, premap[colDel+1])
    -- from 1 to 0 indexing
    vim.api.nvim_buf_set_text(0, r-1, c-1, r-1, c+#pre-1, {table.concat(parts)})
end

-- Delete everything between & and &, as opposed to changeInCell that leaves 
-- whitespace.
-- register: name of register to store the contents of the deleted cell. 
-- Default=unnamed register.
function M.deleteInCell(register)
    register = register or ""
    if not in_env("table") then return end
    local tab = parseTable()
    local row, col = getCurrentCell(tab)
    for i = 1, #tab.texts do
        if row == tab.rows[i] and col == tab.cols[i] then
            local r = tab.rs[i]
            local c = tab.cs[i]
            -- don't delete indent
            local presuf = tab.presufs[i]
            if presuf:match("^[\t%s]+$") then c=c+#presuf end
            -- find nearest &, \\ or newline
            local line = vim.api.nvim_buf_get_lines(0, r-1, r, true)[1]
            local c2 = c + #line:sub(c+1):gsub("\\\\$", ""):gsub("\\&", "  "):gsub("&.*", "")
            -- r-1 for index conv
            vim.api.nvim_buf_set_text(0, r-1, c, r-1, c2, {})
            vim.api.nvim_win_set_cursor(0, {r, c})
            vim.fn.setreg(register, tab.texts[i])
            return true
        end
    end
end

-- select contents without any padding spaces
local function selectInCell()
    util.end_visual()
    local tab = parseTable()
    local row, col = getCurrentCell(tab)
    for i = 1, #tab.texts do
        if row == tab.rows[i] and col == tab.cols[i] then
            local text = tab.texts[i]
            local r = tab.rs[i]
            local c = tab.cs[i]
            -- don't include indent or leading single space
            local presuf = tab.presufs[i]
            if presuf:match("^[\t%s]+$") then c=c+#presuf
            elseif text ~= "" then c=c+1 end
            vim.api.nvim_win_set_cursor(0, {r, c})
            vim.cmd.normal 'v'
            -- assuming single line
            vim.api.nvim_win_set_cursor(0, {r, c + #text - 1})
            return true
        end
    end
end

-- change text in cell except for single surrounding spaces.
-- register: name of register to store the contents of the deleted cell 
-- (without space padding). Default=unnamed register.
function M.changeInCell(register)
    local tab = parseTable()
    local row, col = getCurrentCell(tab)
    for i = 1, #tab.texts do
        if row == tab.rows[i] and col == tab.cols[i] then
            local r = tab.rs[i]
            local c = tab.cs[i]
            local replacement = "  "
            -- don't delete indent and no leading space in that case
            local presuf = tab.presufs[i]
            if presuf:match("^[\t%s]+$") then
                c=c+#presuf
                replacement = " "
            end
            -- find nearest &, \\ or newline
            local line = vim.api.nvim_buf_get_lines(0, r-1, r, true)[1]
            local c2 = c + #line:sub(c+1):gsub("\\\\$", ""):gsub("\\&", "  "):gsub("&.*", "")
            -- r-1 for index conv
            vim.api.nvim_buf_set_text(0, r-1, c, r-1, c2, {replacement})
            vim.api.nvim_win_set_cursor(0, {r, c+#replacement-1})
            vim.fn.setreg(register, tab.texts[i])
            vim.cmd.startinsert()
            return true
        end
    end
end
-- TODO: vix to select contents without any space, cix to change except a 
-- single leading and final space. yix to yank same selection as vix.

local function swapColumn(left, opts)
    left = left or false
    if not in_env("table") then return end
    local tab = parseTable()

    -- first is 0-indexed index of column to swap to the right
    local rowCur, colCur = getCurrentCell(tab)
    local first = colCur
    if left then
        first=first-1
        if first < 0 then return end
    elseif first+1 > tab.maxcol then return end

    local texts2, cols2, presufs2 = {}, {}, {}
    local i = 0
    while i <= #tab.texts do
        i=i+1
        local text = tab.texts[i]
        local col = tab.cols[i]
        local presuf = tab.presufs[i]
        if col == first then
            table.insert(texts2, tab.texts[i+1])
            table.insert(cols2,  col)
            table.insert(presufs2, presuf)
            table.insert(texts2, text)
            table.insert(cols2,  tab.cols[i+1])
            table.insert(presufs2, tab.presufs[i+1])
            i=i+1
        else
            table.insert(texts2, text)
            table.insert(cols2,  col)
            table.insert(presufs2, presuf)
        end
    end
    tab.texts = texts2
    tab.cols = cols2
    tab.presufs = presufs2
    
    writeTable(tab, opts)

    -- also mod preamble
    local r, c, pre, parts, premap = parsePre()
    local leftPart = parts[premap[first+1]] -- +1 for index conv
    parts[premap[first+1]] = parts[premap[first+2]]
    parts[premap[first+2]] = leftPart
    -- from 1 to 0 indexing
    vim.api.nvim_buf_set_text(0, r-1, c-1, r-1, c+#pre-1, {table.concat(parts)})

    -- move cursor to stay on the same column. Reparse, since there may have 
    -- been changes.
    if left then colCur=colCur-1
    else colCur=colCur+1 end
    gotoCell(parseTable(), rowCur, colCur)
end

-- simply add empty cells and move cursor to preamble location to insert
function M.addCol(tab, index, opts)
    local texts = {}
    local cols = {}
    local presufs = {}

    if index <= tab.maxcol then

        for i = 1, #tab.texts do
            local col = tab.cols[i]
            local presuf = tab.presufs[i]
            if col == index then
                -- insert the new column here and shift the current one by 1
                table.insert(texts, "")
                table.insert(cols, index)
                -- steal the prefix in case it is indent
                if presuf:match('\n$') then
                    table.insert(presufs, "")
                else
                    table.insert(presufs, presuf)
                    presuf = ""
                end
                col=col+1
            elseif col > index then
                -- all later columns shifted 1
                col=col+1
            end
            table.insert(texts, tab.texts[i])
            table.insert(cols, col)
            table.insert(presufs, presuf)
        end

    else
        -- also have to be able to add new column at very end
        for i = 1, #tab.texts do
            local col = tab.cols[i]
            local presuf = tab.presufs[i]
            table.insert(texts, tab.texts[i])
            table.insert(cols, col)
            if presuf:match('\\\\\n$') then
                table.insert(presufs, "")
                table.insert(texts, "")
                table.insert(cols, index)
                table.insert(presufs, presuf)
            else
                table.insert(presufs, presuf)
            end
        end
    end

    --tab={r1=int, r2=int, maxcol=int, texts=string[], cols=int[], presufs=string[] }
    -- Remaining parsed entries do not need to be updated.
    tab.texts = texts
    tab.cols = cols
    tab.presufs = presufs
    tab.maxcol=tab.maxcol+1
    writeTable(tab, opts)

    -- jump cursor to preamble to describe new column
    local r, c, pre, parts, premap = parsePre()
    if index > tab.maxcol-1 then -- -1 since the value has just been updated
        for i = 1, premap[#premap] do
            c=c+#parts[i]
        end
    else
        for i = 1, premap[index+1]-1 do
            c=c+#parts[i]
        end
    end
    vim.api.nvim_win_set_cursor(0, {r, c-1})
    vim.cmd.startinsert()
end

-- TODO: deal with \& and maybe add more tabs after multicolumn (which we will have to remove later in the inverse function.)
-- also, is the \\ missing after toprule etc? if so, edit the snippet template.
vim.keymap.set("n", "<plug>TableYank", function()
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
end, {desc="Yank table"})

vim.keymap.set("n", "<plug>TablePaste", function()
    local delim = '\t'
    -- TODO: inverse of above
end, {desc="Paste table"})


vim.keymap.set(
    "n", "<plug>TableAlign", alignTable,
    { silent=true, desc="Align columns in tex table", }
)
vim.keymap.set(
    "n", "<plug>TableDelCol", deleteColumn,
    { silent=true, desc="Delete current column (also from preamble)", }
)
vim.keymap.set(
    "n", "<plug>TableSwapRight", swapColumn,
    { silent=true, desc="Swap table column right (also mod preamble)", }
)
vim.keymap.set(
    "n", "<plug>TableSwapLeft", function () swapColumn(true) end,
    { silent=true, desc="Swap table column left (also mod preamble)", }
)
vim.keymap.set("n", "<plug>TableGoLeft", function ()
    local count = vim.v.count
    if count == 0 then count = 1 end
    local tab = parseTable()
    local row, col = getCurrentCell(tab)
    gotoCell(tab, row, col-count)
end)
vim.keymap.set("n", "<plug>TableGoRight", function ()
    local count = vim.v.count
    if count == 0 then count = 1 end
    local tab = parseTable()
    local row, col = getCurrentCell(tab)
    gotoCell(tab, row, col+count)
end)
vim.keymap.set("n", "<plug>TableGoUp", function ()
    local count = vim.v.count
    if count == 0 then count = 1 end
    local tab = parseTable()
    local row, col = getCurrentCell(tab)
    gotoCell(tab, row-count, col)
end)
vim.keymap.set("n", "<plug>TableGoDown", function ()
    local count = vim.v.count
    if count == 0 then count = 1 end
    local tab = parseTable()
    local row, col = getCurrentCell(tab)
    gotoCell(tab, row+count, col)
end)

vim.keymap.set('n', '<Plug>TableAddColLeft', function ()
    local tab = parseTable()
    local _, col = getCurrentCell(tab)
    M.addCol(tab, col)
end, { desc="Add new empty column to the left" })

vim.keymap.set('n', '<Plug>TableAddColRight', function ()
    local tab = parseTable()
    local _, col = getCurrentCell(tab)
    M.addCol(tab, col+1)
end, { desc="Add new empty column to the right" })

vim.keymap.set('v', '<Plug>TableSelInCell', selectInCell, { silent=true, desc="Select in cell" })

return M

