local util = require "utils/init"

vim.opt_local.conceallevel = 1
-- adoc is for typing prose. But the autowrap works poorly for lists and most
-- syntax beyond basic prose.
vim.opt_local.wrap = true
-- NOTE: things that depend on vim-asciidoc being loaded (compiler, 'comments'
-- override, buffer keymaps, plugin-specific highlight groups) live in
-- after/ftplugin/asciidoc.lua so they run once the plugin's own ftplugin has.

local function get_xref_id()
    local _, c = util.get_cursor()
    local line = vim.api.nvim_get_current_line()
    -- check for <<id,text>> pattern
    local xref = util.match_covering(line, c + 1, '%b<>')
    if xref then return xref:match('^<<([a-z-_]+).*>>$') end
    -- check for xref:id[text] pattern
    xref = util.match_covering(line, c + 1, 'xref:[a-z-_]+%b[]')
    if xref then return xref:match('^xref:([a-z-_]+)') end
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

vim.keymap.set('n', 'gd', goto_xref_tag, { buf=0, desc="Goto tag definition" })

vim.keymap.set('i', '<S-CR>', " +<CR>", {
    buf=0,
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
-- require"completion.asciidoc.cmp_asciidoc".setup()
