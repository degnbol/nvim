#!/usr/bin/env lua
local function hl(name, val)
    vim.api.nvim_set_hl(0, name, val)
end
local function fg(name, color)
    hl(name, {fg=color})
end
local function bg(name, color)
    hl(name, {bg=color})
end
local function fgbg(name, fgcol, bgcol)
    hl(name, {fg=fgcol, bg=bgcol})
end
local function link(name, linkto)
    hl(name, {link=linkto})
end
local function rev(name)
    hl(name, {reverse=true})
end
function gethl(name)
    return vim.api.nvim_get_hl(0, {name=name, link=false})
end
-- update subset of settings for a highlight group instead of replacing them all
function modhl(name, val)
    hl(name, vim.tbl_extend("force", gethl(name), val))
end

---- some sane defaults that a colorscheme may overwrite

-- Italic highlight group doesn't actually make terminal text italic by default.
hl("Italic", {italic=true})

-- remove bg on concealchars
bg("Conceal", nil)

-- set completion menu bg to main bg and make scrollbar minimal
bg("Pmenu", nil)
bg("PmenuSbar", nil)
rev("PmenuSel")

fg("LineNr", "grey")
fg("CursorLineNr", nil)

-- @variable seemed to start being linked to @identifier after update?
hl("@variable", {gui=nil})


local function afterColorscheme()
    -- fix issue where colorscheme change removes all telescope highlight groups
    vim.cmd 'silent Lazy reload telescope.nvim'

    ---- git signs link to Gitsigns equivalent which gets defined
    link("DiffAddNr", "GitsignsAddLn")
    link("DiffChangeNr", "GitsignsChangeLn")
    link("DiffTopDeleteNr", "GitsignsDeleteVirtLn") -- edge-case where first line(s) of file is deleted
    link("DiffChangeDeleteNr", "GitsignsChangedeleteLn")
    -- special decides the color for the underline
    hl("DiffDeleteNr", {underline=true, special=gethl("DiffDelete")["fg"], fg=nil})

    local spellBad = gethl("spellBad")
    -- the default already has red undercurl so we don't wanna mess with that but 
    -- other themes sets the fg instead
    if spellBad["special"] then
        spellBad = spellBad["special"]
    else
        spellBad = spellBad["fg"]
    end
    hl("SpellBad", {undercurl=true, special=spellBad})

    -- replace whatever search was doing completely but only modify IncSearch's 
    -- reverse setting so it remains a different color from Search, probably 
    -- green.
    rev("Search")
    modhl("IncSearch", {reverse=true})

    ---- bufferline
    -- bold instead of italic+bold selected
    hl("BufferLineBufferSelected", {bold=true, italic=false})
    hl("BufferLineNumbersSelected", {bold=true, italic=false})

    ---- Completion. Some links to IncSearch which is too distracting
    hl("CmpItemAbbrMatch", {bold=true})
    hl("CmpItemAbbrMatchFuzzy", {bold=true})

    -- never italic comments but italize builtin stuff
    -- :h group-name
    -- Italics are not default since terms aren't assumed to support it,
    -- however bold is sometimes set for highlight groups for some languages.
    -- italic for builtin reserved words makes sense since they are few and 
    -- don't need to be too legible.
    -- highlight groups starting with @ are from treesitter and others from regex 
    -- syntax groups. Treesitter highlights take priority over regex.
    -- Which one is edited matters, since editing one with color 
    -- will combine the color with the italic, while editing the other may replace 
    -- the color with italic. Before making changes look at test files in 
    -- testfiles/
    modhl("Comment", {italic=false})
    modhl("Operator", {bold=true})
    modhl("Include", {italic=true})
    modhl("Repeat", {italic=true})
    modhl("Label", {italic=true})
    modhl("Conditional", {italic=true})
    modhl("Exception", {italic=true})
    modhl("@include", {italic=true})
    modhl("Keyword", {italic=true, bold=false})
    modhl("@keyword", {italic=true, bold=false})
    modhl("@keyword.function", {italic=true})
    modhl("@keyword.return", {italic=true})
    modhl("@keyword.operator", {italic=true})
    modhl("@parameter", {italic=false})
    modhl("@function", {italic=false})
    modhl("@boolean", {italic=true})

    -- NonText shouldn't be exactly like comments
    if gethl("Comment")['fg'] == gethl("NonText")['fg'] then
        modhl("NonText", {fg="gray"})
    end
    -- it seems pretty good if they're bold even if the color is similar.
    modhl("NonText", {bold=true})
end

local defaultDark = 'fluoromachine'
local defaultLight = 'kanagawa-lotus'

local grp = vim.api.nvim_create_augroup("afterColorscheme", {clear=true})
vim.api.nvim_create_autocmd("Colorscheme", {
    pattern = "*", group = grp,
    callback = afterColorscheme,
})
vim.api.nvim_create_autocmd("VimEnter", {
    pattern = "*", group = grp,
    callback = function ()

        vim.defer_fn(function ()
            -- call twice for the bufferline backgrounds to be set for some reason.
            if vim.o.background == "dark" then
                vim.cmd('colorscheme ' .. defaultDark)
                vim.cmd('colorscheme ' .. defaultDark)
            else
                vim.cmd('colorscheme ' .. defaultLight)
                vim.cmd('colorscheme ' .. defaultLight)
            end
        end, 0) -- Not sure why it needs a 0 ms delay.
    end
})

-- commands :Light and :Dark
vim.api.nvim_create_user_command("Dark", function ()
    vim.fn.system("~/dotfiles/dark.sh")
    vim.cmd('colorscheme ' .. defaultDark)
end, {})
vim.api.nvim_create_user_command("Light", function ()
    vim.fn.system("~/dotfiles/light.sh")
    vim.cmd('colorscheme ' .. defaultLight)
end, {})

