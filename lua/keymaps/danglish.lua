local util = require "utils/init"
local map = require "utils/keymap"

-- Danglish support. For when Danglish keyboard is selected.
-- Generally you should instead stay in code keyboard and use iminsert=2
-- This can also be done with langmap but since these are
-- never used in normal mode then it doesn't hurt, and I tried and it didn't
-- seem to map to [ with remapping even though langremap was on so I don't know
map({ 'n', 'x' }, "æ", ";", { remap = true })
map({ 'n', 'x' }, "Æ", ":", { remap = true })
map({ 'n', 'x' }, "ø", "'", { remap = true })
map({ 'n', 'x' }, "Ø", '"', { remap = true })
map({ 'n', 'x' }, "å", "[", { remap = true })
map({ 'n', 'x' }, "Å", "{", { remap = true })
-- In Danglish I moved : and ; to the }] button
-- But this messes with things when I'm not in Danglish.

-- When language is set to danish we can get ;:'" with alt the same way as we
-- do in danglish keyboard input, however that is with the left alt which is
-- deeper in the OS and the mapping below is for when esc+key (^[) is detected which
-- is what is sent to the terminal, e.g. with kitty's setting 'macos_option_as_alt right'.
-- Note that they are also available with the ]} and \| keys.
map.i("<A-;>", ";")
map.i("<A-S-;>", ":")
map.i("<A-'>", "'")
map.i("<A-S-'>", '"')
-- when writing text with Danish we might try to write : but the key is mapped to Æ.
-- : is written at ends of word where we would never write capital Æ, so we can check if we are at end of word.
-- The only exception would be if the entire word is uppercase. Currently choosing to ignore that edge case.
map.i("Æ", function()
    local char = util.get_current_char()
    local put = char:match('[A-Åa-å0-9.,!?]') and ':' or 'Æ'
    util.put_char(put)
end)
-- similarly we often would want " instead of Ø, e.g. if we write ØØ it's to make "" and
map.i("Ø", function()
    local r, c = util.get_cursor()
    local char, c1 = util.get_char(r, c)
    if char == 'Ø' then
        vim.api.nvim_buf_set_text(0, r, c1, r, c, { '""' })
        -- try to see if it's ok, otherwise change to only looking for uppercase.
    elseif util.get_char(r, c + 1):match('[A-Åa-å]') then
        util.put_char('"')
    else
        local put = char:match('[A-Åa-å.,!?]') and '"' or 'Ø'
        util.put_char(put)
    end
end)
-- for å we might want [] or {}, but with <C-6> pressed they're mapped to å; and Å:,
-- however it's a remap from ] and } so we write that.
map.i('å]', '[]')
map.i('Å}', '{}')

local function toggle_danglish(silent)
    if not silent then
        if vim.bo.iminsert == 0 then
            print("Danglish")
        else
            print("Code")
        end
    end
    if util.get_mode() == 'n' then
        util.press("a<C-^><Esc>")
    else
        util.press("<C-^>")
    end
    -- other related things to get triggered when toggling language.
    -- Has to be scheduled since it checks iminsert which isn't updated immediately.
    vim.schedule(function()
        vim.api.nvim_exec_autocmds("User", { pattern = "ToggleDansk" })
    end)
end

local function toggle_danish_imaps(silent)
    local is_mapped = vim.g.danish_imaps
    if not is_mapped then
        vim.g.danish_imaps = true
        map.i('ae', 'æ', "ae -> æ")
        map.i('oe', 'ø', "oe -> ø")
        map.i('aa', 'å', "aa -> å")
        map.i('Ae', 'Æ', "Ae -> Æ")
        map.i('Oe', 'Ø', "Oe -> Ø")
        map.i('Aa', 'Å', "Aa -> Å")
        map.i('AE', 'Æ', "AE -> Æ")
        map.i('OE', 'Ø', "OE -> Ø")
        map.i('AA', 'Å', "AA -> Å")
        if not silent then print("ae -> æ, osv.   aktiveret") end
    else
        vim.g.danish_imaps = false
        vim.keymap.del('i', 'ae')
        vim.keymap.del('i', 'oe')
        vim.keymap.del('i', 'aa')
        vim.keymap.del('i', 'Ae')
        vim.keymap.del('i', 'Oe')
        vim.keymap.del('i', 'Aa')
        vim.keymap.del('i', 'AE')
        vim.keymap.del('i', 'OE')
        vim.keymap.del('i', 'AA')
        if not silent then print("ae -> æ, osv. DEaktiveret") end
    end
end

map.i('<C-6>', toggle_danglish, "Toggle dansk")
-- Ins key was a bit useless just doing what i does so let's make it a language switch insertion:
map.n('<Ins>', 'i<C-6>', "Insert + Toggle dansk", { remap = true })
-- switch to/from Danish æøå and to insert mode, which is convenient.
-- remap in order to utilise the remapped <C-6> which updates the cmp dictionary
map.n("yod", "i<C-6>", "Danish (<C-^>)", { remap = true })

map.i("<C-S-6>", toggle_danish_imaps, "Toggle Danish imaps")
map.n("<leader>D", function()
    toggle_danish_imaps(true)
    toggle_danglish(true)
    if vim.g.danish_imaps then
        print("Danglish + ae -> æ, osv.   aktiveret")
    else
        print("Code     + ae -> æ, osv. DEaktiveret")
    end
end, "Toggle Danish imaps + Danglish")

-- And in case danglish keyboard is active:
map.i("<A-æ>", ";")
map.i("<A-S-æ>", ":")
map.i("<A-ø>", "'")
map.i("<A-S-ø>", '"')
