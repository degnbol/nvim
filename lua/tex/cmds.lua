
local util = require "utils/init"
local vtu = require "utils/vimtex"

local before = util.before
local get_visual_range = util.get_visual_range
local in_math = vtu.in_math

-- bold, italics, and similar.
-- TODO: color text, where color text should probs be <C-something>r for red, and so on.
-- <c-c> is already for cancel and <c-u> for delete to start of line. These are uncommon so maybe just uppercase U and C?

local M = {}

local cmd_mappings = {
    bold = {
        imap = "<D-b>",
        text = "textbf",
        -- unlike mathbf that only bolds latin, symbf with unicode-math works with greek.
        math = "symbf"
    },
    italic = {
        imap = "<D-i>",
        text = "textit",
        -- useful to reduce spacing slightly in multi letter variable names
        math = "mathit"
    },
    emphasis = {
        imap = "<D-e>",
        text = "emph",
    },
    upright = {
        imap = "<C-S-u>",
        text = "normalfont",
        -- Unlike mathrm that only un-italizes latin, symrm with unicode-math works with greek.
        math = "symrm"
    },
    teletype = {
        imap = "<C-S-t>",
        text = "texttt",
    },
    underline = {
        imap = "<C-S-->",
        text = "underline",
    },
    color = {
        imap = "<C-S-c>",
        text = "textcolor",
        -- There may be a need to redefine the textcolor cmd in math, see
        -- https://tex.stackexchange.com/questions/21598/how-to-color-math-symbols
    }
}

-- adjust a c(olumn) to before a cmd if it is inside it.
-- Implication is we don't split a cmd.
-- indexing is (1,0)
local before_cmd = function (r, c)
    local line = vim.fn.getline(r)
    -- +1 to include the char right after cursor bar
    local cmdStart = line:sub(1, c+1):match('\\%a+$')
    if cmdStart ~= nil then
        return c-#cmdStart+1 -- +1 to undo the +1
    end
    return c
end

-- max number of lines to look before cursor for an un-closed cmd
local CONTEXT = 20

local _insert_or_del_cmd = function(name)
    local cmd = vtu.inside_cmd(name)
    if cmd ~= nil then
        return vim.fn['vimtex#cmd#delete'](cmd[1], cmd[2])
    end

    name = '\\' .. name

    local r, c = unpack(vim.api.nvim_win_get_cursor(0))

    -- if the cmd is just opened, it will not be detected by vimtex
    local iDel, jDel = nil, nil
    -- look at current line and then a few lines before
    for rDel = r, math.max(1, r-CONTEXT), -1 do
        local line = vim.fn.getline(rDel)
        if rDel == r then line = line:sub(1, c + #name + 1) end
        for i,j in line:gmatch('()' .. name .. '{?()') do
            local cmd = vim.fn['vimtex#cmd#get_at'](rDel, i)
            if vim.tbl_isempty(cmd) or cmd["name"] ~= name then
                iDel, jDel = i, j
                -- keep looping so the last match is the one that gets replaced in 
                -- the edge-case where there are multiple
            end
        end
        if iDel ~= nil then
            -- -1 for 0-index
            vim.api.nvim_buf_set_text(0, rDel-1, iDel-1, rDel-1, jDel-1, {})
            if r == rDel then vim.api.nvim_win_set_cursor(0, {r, c-(jDel-iDel)}) end
            return
        end
    end

    -- are we inside a (any) cmd?
    local cCurs = c
    -- we specically DON'T want to use the vimtex version
    local c = before_cmd(r, c)
    vim.api.nvim_buf_set_text(0, r-1, c, r-1, c, {name .. '{'})
    if c == cCurs then
        -- move cursor to after the insertion
        vim.api.nvim_win_set_cursor(0, {r, c+ #name + 1})
    end
end

-- TODO: handle square brackets, e.g. in desc env merging \item[\texttt{reps1}, \texttt{reps2}] without surrounding the \item.
local function cmd_contained(cmd, r1, c1, r2, c2)
    local args = cmd["args"]
    -- args is empty for simple commands without {}
    if vim.tbl_isempty(args) then return false end
    args = args[1] -- always one element?
    local ro, co = args["open"]["lnum"], args["open"]["cnum"]
    local rc, cc = args["close"]["lnum"], args["close"]["cnum"]
    return before(ro, co, r1, c1) and before(r2, c2, rc, cc)
end

-- surround a selection with a cmd that is joined if overlapping, e.g. textbf 
-- where it isn't meaningful to bold text twice.
local function _surround_visual(name)
    local r1, c1, r2, c2 = get_visual_range()
    -- (1,1)-indexed.
    c1,c2=c1+1,c2+1
    -- extend start of selection if it is in the middle of a (any) cmd
    -- EXCEPT if it is fully contained
    local cmd1 = vim.fn['vimtex#cmd#get_at'](r1, c1)
    if not vim.tbl_isempty(cmd1) and not cmd_contained(cmd1, r1, c1, r2, c2) then
        r1, c1 = cmd1["pos_start"]['lnum'], cmd1["pos_start"]['cnum']
    end
    local cmd2 = vim.fn['vimtex#cmd#get_at'](r2, c2)
    if not vim.tbl_isempty(cmd2) and not cmd_contained(cmd2, r1, c1, r2, c2) then
        r2, c2 = cmd2["pos_end"]['lnum'], cmd2["pos_end"]['cnum']
    end

    -- find cmds already present.
    -- We expand the "search" by 1 so that we consume (merge with) adjacent cmds as well.
    -- the max on c1 is needed. It would be easy enough to search wrapping lines as well, e.g. if c1==1 then r1-=1, c1=EOL end
    local cmds = vtu.inside_cmd_range(name, r1, math.max(1,c1-1), r2, c2+1)
    -- if any go beyond the selection we expand the selection to merge them
    for _, cmd in ipairs(cmds) do
        local _r1, _c1, _r2, _c2 = unpack(cmd)
        if before(_r1, _c1, r1, c1) then
            r1, c1 = _r1, _c1
        end
        if before(r2, c2, _r2, _c2) then
            r2, c2 = _r2, _c2
        end
    end
    -- easier to modify and replace the selection text rather than the buffer 
    -- when changes might affect indexes for other cmds.
    -- -1 for 0-index, end-exclusive so not on r2
    local lines = vim.api.nvim_buf_get_lines(0, r1-1, r2, true)
    local shift, r_last = 0, r1 -- keep track of the shift in the next indexes
    for _, cmd in ipairs(cmds) do
        local _r1, _c1, _r2, _c2 = unpack(cmd)
        if _r1 ~= r_last then shift = 0 else -- reset on new line
            _c1=_c1-shift
        end
        -- remove beginning
        lines[_r1-r1+1] = lines[_r1-r1+1]:sub(1, _c1-1) .. lines[_r1-r1+1]:sub(_c1 + #name + 2)
        if _r2 ~= _r1 then shift = 0 else
            shift=shift+#name+2
            _c2=_c2-shift
        end
        -- remove end
        lines[_r2-r1+1] = lines[_r2-r1+1]:sub(1, _c2-1) .. lines[_r2-r1+1]:sub(_c2+1)
        -- keep note for next iteration
        shift=shift+1
        r_last = _r2
    end
    -- then surround the whole selection with the cmd
    lines[1] = lines[1]:sub(1, c1-1) .. '\\' .. name .. '{' .. lines[1]:sub(c1)
    -- if we are adding to a single line then the end will be shifted opposite
    if #lines == 1 then c2=c2 - shift + #name + 2
    elseif r2 == r_last then -- we may be multiple lines but have already modified the last line
        c2=c2-shift
    end
    lines[#lines] = lines[#lines]:sub(1, c2) .. '}' .. lines[#lines]:sub(c2+1)
    vim.api.nvim_buf_set_lines(0, r1-1, r2, true, lines)
end

-- set keymaps

local function _get_cmd_name(cmd_grp)
    if in_math() then
        return cmd_grp.math or cmd_grp.text
    else
        return cmd_grp.text
    end
end

function M.map_keys()
    for desc, cmd_grp in pairs(cmd_mappings) do
        local cmd_name = _get_cmd_name(cmd_grp)
        if cmd_grp.imap ~= nil then
            vim.keymap.set("i", cmd_grp.imap, function ()
                _insert_or_del_cmd(cmd_name)
            end, {buffer=true, desc=desc})
        end
        local xmap = cmd_grp.xmap or cmd_grp.imap
        if xmap ~= nil then
            vim.keymap.set("x", xmap, function ()
                _surround_visual(cmd_name)
            end, {buffer=true, desc=desc})
        end
    end
end

return M
