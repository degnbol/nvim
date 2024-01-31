
function tabNew()
    count = vim.v.count == 0 and 5 or vim.v.count

    local lines = {
        "e|" .. ('-'):rep(count) .. "|",
        "B|" .. ('-'):rep(count) .. "|",
        "G|" .. ('-'):rep(count) .. "|",
        "D|" .. ('-'):rep(count) .. "|",
        "A|" .. ('-'):rep(count) .. "|",
        "E|" .. ('-'):rep(count) .. "|",
    }

    vim.api.nvim_put(lines, "l", true, false)
end

local function isTabLine(line)
    return line:match("^%a?|%S*|")
end
function inTab()
    return isTabLine(vim.api.nvim_get_current_line())
end
---Get index of guitar string.
---@return integer: in range [1;6]
local function getCurrentString()
    local curs, _ = unpack(vim.api.nvim_win_get_cursor(0))
    for r = curs-1, 1, -1 do
        local line = vim.api.nvim_buf_get_lines(0, r-1, r, true)[1]
        if not isTabLine(line) then
            return curs - r
        end
    end
    return curs
end

-- TODO: not in visual block mode
function tabDash()
    if inTab() then
        local r, c = unpack(vim.api.nvim_win_get_cursor(0))
        vim.api.nvim_win_set_cursor(0, {r - getCurrentString() + 1, c})
        vim.api.nvim_put({'-', '-', '-', '-', '-', '-'}, 'b', false, false)
        vim.api.nvim_win_set_cursor(0, {r, c+1})
    else
        vim.api.nvim_put({'-'}, 'c', true, true)
    end
end
function tabBS()
    if inTab() then
        local r, c = unpack(vim.api.nvim_win_get_cursor(0))
        local n = getCurrentString()
        local lines = vim.api.nvim_buf_get_lines(0, r-n, r-n+6+1, true)
        local replacement = {}
        for i, line in ipairs(lines) do
            replacement[i] = line:sub(1,c-1) .. line:sub(c+1)
        end
        vim.api.nvim_buf_set_lines(0, r-n, r-n+6+1, true, replacement)
        vim.api.nvim_win_set_cursor(0, {r, c-1})
    else
        local keys = vim.api.nvim_replace_termcodes('<BS>', true,false,true)
        vim.api.nvim_feedkeys(keys, 'n', false)
    end
end

function tabSelect()
    local n = getCurrentString() - 1
    if n > 0 then
        return n .. "k<C-v>5jo"
    else
        return       "<C-v>5jo"
    end
end
