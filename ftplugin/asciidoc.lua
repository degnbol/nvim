

-- instead of to Constant, which is currently italic.
vim.api.nvim_set_hl(0, "AsciidoctorCode", {link="@string", default=true})

local function get_xref_id()
    local r, c = unpack(vim.api.nvim_win_get_cursor(0))
    local c=c+1 -- 1-index
    local line = vim.api.nvim_get_current_line()
    -- check for <<id,text>> pattern
    for startpos, xref_match, endpos in line:gmatch('()(%b<>)()') do
        if startpos <= c and c <= endpos then
            local tag = xref_match:match('^<<([a-z-_]+).*>>$')
            return tag
        end
    end
    -- check for xref:id[text] pattern
    for startpos, tag, endpos in line:gmatch('()xref:([a-z-_]+)%b[]()') do
        if startpos <= c and c <= endpos then
            return tag
        end
    end
end
-- Goto tag that xref points to under cursor
local function goto_xref_tag()
    local tag = get_xref_id()
    if tag then
        if not pcall(vim.cmd, "tag " .. tag) then
            -- assuming we have not changed default ctags call that stores tags by 
            -- natural naming rather than the ids, we can simply modify the search 
            -- here to work for both natural naming and id naming (:tag is 
            -- case-insensitive)
            vim.cmd("silent! tag " .. tag:gsub("_", " "))
        end
    end
end

vim.keymap.set('n', 'gd', goto_xref_tag, { desc="Goto tag definition" })

vim.keymap.set('i', '<S-CR>', " +<CR>", {
    desc=[[Hard line break, similar to \\ in tex.
    Can also be achieved with paragraph option [%hardbreaks] or document option :hardbreaks-option:]]
}
)

-- if :hardbreaks-option: is on we shouldn't auto break lines for obvious 
-- reasons, which also means we most likely want to wrap lines.
local lines = vim.api.nvim_buf_get_lines(0, 0, -1, false)
for _, line in ipairs(lines) do
    if line == ":hardbreaks-option:" then
        vim.opt_local.formatoptions:remove('a')
        vim.opt_local.wrap = true
        break
    end
end

-- Using blink.cmp instead.
-- TODO: write the blink_ equivalent file if LSP hasn't caught up.
-- require"completion.asciidoc.cmp_asciidoc".setup()
