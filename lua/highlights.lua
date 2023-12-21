#!/usr/bin/env lua
local hl = require "utils/highlights"

---- some sane defaults that a colorscheme may overwrite

-- mods to default nvim colorscheme
-- TODO: maybe move these to a colorscheme that is an extension of default.
hl.set("@conditional",      {fg="#9a808f", italic=true})
hl.set("@repeat",           {fg="#9a808f", italic=true})
hl.set("@include",          {fg="#9a808f", italic=true})
hl.set("Exception",         {fg="#9a808f", italic=true})
hl.set("Keyword",           {fg="#9a808f", italic=true})
hl.set("Operator",          {fg="#9a808f", bold=true})
hl.set("@keyword.function", {fg="#9a808f"})
hl.set("Type",              {fg="#9a808f", bold=true})
hl.set("Preproc",           {fg="#8888ff"})
-- hl.set("Type", {fg="#88ef88", bold=true})
hl.link("Identifier", "@variable")
hl.link("Special", "PreProc")
hl.link("Number", "String")
hl.link("Constant", "String")
hl.link("@constant.builtin", "String")

-- Italic highlight group doesn't actually make terminal text italic by default.
hl.set("Italic", {italic=true})

-- remove bg on concealchars
hl.bg("Conceal", nil)

-- set completion menu bg to main bg and make scrollbar minimal
hl.bg("Pmenu", nil)
hl.bg("PmenuSbar", nil)
hl.rev("PmenuSel")

hl.fg("LineNr", "grey")
hl.fg("CursorLineNr", nil)

-- @variable seemed to start being linked to @identifier after update?
hl.set("@variable", {gui=nil})


local function afterColorscheme()
    -- fix issue where colorscheme change removes all telescope highlight groups
    -- vim.cmd 'silent Lazy reload telescope.nvim'

    ---- git signs link to Gitsigns equivalent which gets defined
    hl.link("DiffAddNr", "GitsignsAddLn")
    hl.link("DiffChangeNr", "GitsignsChangeLn")
    hl.link("DiffTopDeleteNr", "GitsignsDeleteVirtLn") -- edge-case where first line(s) of file is deleted
    hl.link("DiffChangeDeleteNr", "GitsignsChangedeleteLn")
    -- special decides the color for the underline
    hl.set("DiffDeleteNr", {underline=true, special=hl.get("DiffDelete")["fg"], fg=nil})

    local spellBad = hl.get("spellBad")
    -- the default already has red undercurl so we don't wanna mess with that but 
    -- other themes sets the fg instead
    if spellBad["special"] then
        spellBad = spellBad["special"]
    else
        spellBad = spellBad["fg"]
    end
    hl.set("SpellBad", {undercurl=true, special=spellBad})

    -- replace whatever search was doing completely but only modify IncSearch's 
    -- reverse setting so it remains a different color from Search, probably 
    -- green.
    hl.rev("Search")
    hl.mod("IncSearch", {reverse=true})

    -- highlight bg defeats the purpose of folding for me.
    -- It is still plenty clear that text is folded.
    hl.bg("Folded", nil)

    ---- bufferline
    -- bold instead of italic+bold selected
    hl.set("BufferLineBufferSelected", {bold=true, italic=false})
    hl.set("BufferLineNumbersSelected", {bold=true, italic=false})

    ---- Completion. Some links to IncSearch which is too distracting
    hl.set("CmpItemAbbrMatch", {bold=true})
    hl.set("CmpItemAbbrMatchFuzzy", {bold=true})

    -- there is enough things indicated by color, matched parenthesis etc. with underline is great
    hl.set("MatchParen", {underline=true})

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
    hl.mod("Comment", {italic=false})
    hl.mod("Operator", {bold=true})
    hl.mod("Include", {italic=true})
    hl.mod("Repeat", {italic=true})
    hl.mod("Label", {italic=true})
    hl.mod("Type", {italic=false})
    hl.mod("Conditional", {italic=true})
    -- for e.g. julia there is the `a = bool ? b : c` notation. It's weird to 
    -- have ? and : italic since that is meant for words, but it does help 
    -- distinguish them from : used in e.g. ranges.
    -- hl.mod("@conditional.ternary", {italic=false})
    hl.mod("Identifier", {italic=false})
    hl.mod("Number", {italic=false})
    hl.mod("Exception", {italic=true})
    hl.mod("@include", {italic=true})
    hl.link("@keyword", "Keyword")
    hl.mod("Keyword", {italic=true, bold=false})
    hl.mod("@keyword.function", {italic=true})
    hl.mod("@keyword.return", {italic=true})
    hl.mod("@keyword.operator", {italic=true})
    hl.mod("@parameter", {italic=false})
    hl.mod("@function", {italic=false, bold=true})
    hl.mod("@function.call", {bold=false})
    hl.link("@function.macro", "@function.call")
    hl.mod("Conditional", {italic=true, bold=false})
    hl.mod("@conditional.ternary", {italic=false})
    hl.mod("@repeat", {italic=true, bold=false})
    hl.link("@boolean", "Boolean")
    hl.mod("Boolean", {italic=true, bold=false})
    hl.mod("@variable.builtin", {italic=true})
    hl.mod("@function.builtin", {italic=true})
    hl.mod("@constant.builtin", {italic=true})
    hl.mod("@type.builtin",     {italic=true})

    -- delim. Currently using the default for parentheses.
    hl.link("Delimiter", "RainbowDelimiterViolet")
    hl.link("@punctuation.delimiter", "Delimiter")
    
    -- NonText shouldn't be exactly like comments
    if hl.get("Comment")['fg'] == hl.get("NonText")['fg'] then
        hl.mod("NonText", {fg="gray"})
    end
    -- it seems pretty good if they're bold even if the color is similar.
    hl.mod("NonText", {bold=true})

    -- semantic tokens
    hl.mod("@lsp.mod.strong", {bold=true})
    hl.mod("@lsp.mod.emph", {italic=true})
    hl.link("@lsp.type.operator", "@operator")
    hl.link("@lsp.type.keyword", "@keyword")

    -- variable is the default.
    -- I set it to copy normal instead of using clear since if I clear it will 
    -- be overrideen by other colors, but I want it to appear. Example in julia:
    -- "$variable" will color variable as string with clear and as normal with this approach.
    hl.fg("@variable", hl.get("Normal")["fg"])
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

        vim.schedule(function ()
            -- pretend we called colorscheme in order to trigger all autcmds 
            -- that fire after setting a new colorscheme.
            vim.api.nvim_exec_autocmds("Colorscheme", {})
            -- vim.cmd "hi clear"
            -- if vim.o.background == "dark" then
            --     vim.cmd('colorscheme ' .. defaultDark)
            -- else
            --     vim.cmd('colorscheme ' .. defaultLight)
            -- end
        end)
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

