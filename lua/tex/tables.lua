#!/usr/bin/env lua
local util = require "utils/init"
local vtu = require "utils/vimtex"
local is_inside = vtu.is_inside
local in_env = vtu.in_env
local cmd = vim.cmd

--- Work for either a visual selection or for the current env.
---maxwidth: maximum allowed column width
---usemax: in case a cell exceeds the max, should it be ignored (default), or 
---should the column be set to maxwidth (usemax=true).
function alignTable(opts)
    opts = opts or {}
    maxwidth = opts.maxwidth or 20
    usemax = opts.usemax or false

    if not in_env("table") then return end

    local mode = vim.api.nvim_get_mode().mode
    local r1, r2
    if mode == "n" then
        -- get lines of env using vimtex
        -- should never return nil since we test in_env above
        _, r1, _, r2, _ = unpack(vtu.get_env())
    else
        util.end_visual()
        r1, _, r2, _ = util.get_visual_range()
    end
    -- -1 since we are using 1-indexed r1, r2 in 0-indexed end-exclusive function
    -- to get the lines withinh the end, so excluding the begin/end lines
    local lines = vim.api.nvim_buf_get_lines(0, r1, r2-1, true)
    -- final newline necessary to not lose final cell
    local text = table.concat(lines, '\n') .. '\n'
    local indentation = text:match('^[s%\t]*')
    local cells = {}
    local cellc = {} -- column in table
    local cellcl = {} -- column in line
    local col = 0
    local coll = 0
    local lasti = 0
    -- find the position of each & and \\, the latter assumed to be at end of line.
    -- This assumption avoids \\ inside cell with forced linebreak.
    for i, m in text:gmatch('()([&\n])') do
        local section = text:sub(lasti, i-1)
        if m == '&' and text:sub(i,i) ~= '\\' then
            -- the minus works like *, except match as shortly as possible, which let's us trim whitespace.
            local cell = section:match('\n?%s*([^\n]-)%s*$')
            table.insert(cells, cell)
            table.insert(cellc, col)
            table.insert(cellcl, coll)
            local colspan = cell:match("^\\multicolumn{(%d+)}") or 1
            col=col + colspan
            coll=coll + colspan
            lasti = i+1
        elseif m == '\n' then
            local cell = section:match('\n?%s*([^\n]-)%s*\\\\%s*$')
            -- might not have the \\, e.g. if cell is a \toprule line
            if cell == nil then
                table.insert(cells, section .. '\n')
                table.insert(cellc, -1) -- mark to ignore
                table.insert(cellcl, coll)
            else
                table.insert(cells, cell)
                table.insert(cellc, col)
                table.insert(cellcl, coll)
                col=0
            end
            coll=0
            lasti = i+1
        end
    end
    local maxcol = math.max(unpack(cellc))
    -- column width preallocate zeros
    local widths = {}
    for i = 0, maxcol do widths[i] = 0 end
    for i, col in ipairs(cellc) do
        if col >= 0 then
            local cell = cells[i]
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
    for i, cell in ipairs(cells) do
        local col = cellc[i]
        local coll = cellcl[i]
        if col == -1 then
            table.insert(text, cell)
            idealw, currentw = 0, 0
        else
            idealw=idealw + widths[col]
            local multicol = cell:match("^\\multicolumn{(%d+)}")
            if multicol ~= nil then
                for i = 1, multicol-1 do
                    idealw=idealw + 3 + widths[col+i]
                end
            end
            if coll == 0 then
                cell = indentation .. cell
                idealw=idealw+#indentation
            elseif cell ~= "" then
                cell = " " .. cell
            end
            currentw=currentw + #cell
            if currentw < idealw then
                cell = cell .. string.rep(' ', idealw - currentw)
                currentw = idealw
            end
            if i == #cellc or cellc[i+1] == 0 or col == maxcol then
                cell = cell .. " \\\\\n"
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
    vim.api.nvim_buf_set_lines(0, r1, r2-1, true, lines)
end

vim.keymap.set(
    "n", "<plug>TableAlign", alignTable,
    { silent=true, desc="Align columns in tex table", }
)
-- align &s using plugin (that also aligns \\ since it understands it is for 
-- latex), ignoring lines that contain "multicol".
-- Fails for complicated cells, such as \makecell[l]{a \\ b},
-- hence I wrote the alignTable function above.
vim.keymap.set("n", "<plug>TableAlign2", "gaie<C-f>v/\\multicol<CR>*&", {
    silent=true,
    remap=true,
    desc="Align table",
})


-- TODO:
-- functions or ways of moving/deleting/adding rows/columns.
-- Maybe each a function but ideally a more modular (vim) way, e.g. selection of row/col,
-- moving between cells, etc. But might not be possible without elaborate 
-- function due to things like header row, and not having one row per line.

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

-- get {r, c, string} of preamble for a containing tabular/tabularx/tabulary env.
-- where string is e.g. "@{}p{1.1cm}rrlp{1.1cm}rrrrlX@{}"
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
    local until_cursor = pre:sub(1, c+2-c_pre)
    local nAlpha = #until_cursor:gsub("%b{}", ""):gsub('[^%a]', '')
    -- 0-indexing
    return math.max(0, nAlpha-1)
end

-- TODO: when jumping to 0-th that is empty we go to 0-th column instead of cell indent.
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
        local col = getCurrentColumn(c)
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

-- TODO: let multicolumn align to the rest but not the other way. What to do if not one row per line?
-- TODO: have a config parameter for max cell length, so anything longer will not force other cells to be silly long.
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

-- TODO: <leader>-a and <leader>-A could swap columns like they normally swap args.


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


