local map = require "utils/keymap"

-- Write brackets and other doubled symbols quicker
map.i("<C-'>", "''<left>")
map.i("<C-S-'>", '""<left>')
map.i("<C-9>", '()<left>')
map.i("<C-S-9>", '()<left>')
map.i("<C-0>", '()<left>')
map.i("<C-S-0>", '()<left>')
-- imap("<C-[>", '[]<left>') -- not possible since it's literally ESC
map.i("<C-]>", '[]<left>')
map.i("<C-S-[>", '{}<left>')
map.i("<C-S-]>", '{}<left>')
-- Seem to activate signature help more smoothly than ()<S-space> which may reguire a trigger char such as ,
map.i('<C-b>', "()<left>", "()<left>")

-- hack map of shift+space
local bracketJumpCode = "\x1F"
local paired = { '""', "''", "``", "()", '[]', "{}", "<>", "$$" }
local paireddouble = { '[[]]', "\\{\\}" }     -- lua and tex
local singles = { "'", '"', '`', '(', ')', '[', ']', '{', '}', '<', '>', '$' }
local doubles = { "[[", "]]", "\\{", "\\}", } -- lua and tex
local triples = { '"""', "'''", "```" }
local function bracketJump(line, c)
    if vim.tbl_contains(triples, line:sub(c - 2, c)) then
        return "<left><left><left>"
    elseif vim.tbl_contains(triples, line:sub(c + 1, c + 3)) then
        return "<right><right><right>"
    elseif vim.tbl_contains(paireddouble, line:sub(c - 3, c)) then
        return "<left><left>"
    elseif vim.tbl_contains(paireddouble, line:sub(c + 1, c + 4)) then
        return "<right><right>"
        -- left priority over right for pairs so we go back first for e.g. Matrix{}|[]
    elseif vim.tbl_contains(paired, line:sub(c - 1, c)) then
        return "<left>"
    elseif vim.tbl_contains(paired, line:sub(c + 1, c + 2)) then
        return "<right>"
    elseif vim.tbl_contains(doubles, line:sub(c + 1, c + 2)) then
        return "<right><right>"
    elseif vim.tbl_contains(doubles, line:sub(c - 1, c)) then
        return "<left><left>"
        -- right priority over left for singles, since they are usually half of a
        -- filled out pair and we want to prioritize progressing in that case.
    elseif vim.tbl_contains(singles, line:sub(c + 1, c + 1)) then
        return "<right>"
    elseif vim.tbl_contains(singles, line:sub(c, c)) then
        return "<left>"
    else -- maybe accidentally still held shift
        return " "
    end
end
map({ 'i', 'n' }, bracketJumpCode, function()
    local line = vim.api.nvim_get_current_line()
    local r, c = unpack(vim.api.nvim_win_get_cursor(0))
    return bracketJump(line, c)
end, { expr = true, desc = "Move inside empty pair/triples or outside non-empty" })
map.c(bracketJumpCode, function()
    local line = vim.fn.getcmdline()
    local c = vim.fn.getcmdpos()
    return bracketJump(line, c - 1)
end, "Move inside empty pair/triples or outside non-empty", { expr = true })

-- useful with cursor | in {|} to get
-- {
--     |
-- }
map.i('<S-CR>', "<CR><Esc>O", "Indented newline")
