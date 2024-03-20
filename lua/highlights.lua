#!/usr/bin/env lua
local hl = require "utils/highlights"

local function afterColorscheme()
    -- fix issue where colorscheme change removes all telescope highlight groups
    -- vim.cmd 'silent Lazy reload telescope.nvim'

    ---- git signs link to Gitsigns equivalent which gets defined
    hl.def("DiffAddNr", "GitsignsAddLn")
    hl.def("DiffChangeNr", "GitsignsChangeLn")
    hl.def("DiffTopDeleteNr", "GitsignsDeleteVirtLn") -- edge-case where first line(s) of file is deleted
    hl.def("DiffChangeDeleteNr", "GitsignsChangedeleteLn")
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
    -- cursor looks to much like the current search word under cursor
    hl.set("CurSearch", {fg="gray", bg=hl.get("IncSearch")["fg"], standout=true, bold=true})
    -- hl.mod("CurSearch", {standout=true})

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
    hl.mod("Conditional", {italic=true, bold=false})
    -- for e.g. julia there is the `a = bool ? b : c` notation. It's weird to 
    -- have ? and : italic since that is meant for words, but it does help 
    -- distinguish them from : used in e.g. ranges.
    -- hl.mod("@conditional.ternary", {italic=false})
    hl.mod("Identifier", {italic=false})
    hl.mod("Number", {italic=false})
    -- vim in lua is @lsp.typemod.variable.global.lua linked to Constant.
    -- If you find other constants that you don't want to make italic then mod the semantic @lsp global instead.
    hl.mod("Constant", {italic=true})
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
    hl.mod("@conditional", {italic=true, bold=false})
    hl.mod("@conditional.ternary", {italic=false})
    hl.mod("@repeat", {italic=true, bold=false})
    hl.link("@boolean", "Boolean")
    hl.mod("Boolean", {italic=true, bold=false})
    hl.mod("@variable.builtin", {italic=true})
    hl.mod("@function.builtin", {italic=true, bold=false})
    hl.mod("@constant.builtin", {italic=true})
    hl.mod("@type.builtin",     {italic=true})

    -- delim. Currently using the default for parentheses.
    hl.link("Delimiter", "RainbowDelimiterViolet")
    hl.link("@punctuation.bracket", "RainbowDelimiterViolet")
    hl.link("@punctuation.delimiter", "Delimiter")
    hl.link("@punctuation.special", "Special")
    
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
    -- instead of to @function since we only want function definitions to be 
    -- bold and @type.function is often not, e.g. in this very file.
    hl.link("@lsp.type.function", "@function.call")
    hl.link("@lsp.type.method", "@function.call")
    hl.link("@lsp.type.string", "@string")

    --- Fix lua. Colon before function call is captured by @constant making it italic.
    if vim.bo.filetype == "lua" then
        hl.set("@constant", {})
        -- remap things that were mapped to @constant to preserve their functioning highlights
        hl.link("@lsp.typemod.variable.global.lua", "Constant")
        hl.link("@lsp.typemod.variable.defaultLibrary.lua", "Constant")
    end

    -- variable is the default capture for most things in code so we want it to 
    -- be neutral, although not if the language is markdown etc.
    local showVar = {typst=true, markdown=true, asciidoc=true, latex=true}
    if showVar[vim.bo.filetype] then
        hl.fg("@variable", hl.get("@variable.builtin")["fg"])
    else
        -- I set it to copy normal instead of using clear since if I clear it will 
        -- be overrideen by other colors, but I want it to appear. Example in julia:
        -- "$variable" will color variable as string with clear and as normal with this approach.
        hl.fg("@variable", hl.get("Normal")["fg"])
    end

    -- extmarks
    -- By default colors. Underline variants makes more sense.
    hl.set("DiagnosticUnderlineInfo", {})
    -- subtle. underline and underdashed are stronger but the warn is often 
    -- wrong, e.g. missing reference to things LSP doesn't understand is 
    -- imported.
    hl.set("DiagnosticUnderlineWarn", {underdotted=true})
    hl.set("DiagnosticUnderlineError", {undercurl=true})
end

-- local defaultDark = 'fluoromachine'
local defaultDark = 'delta'
-- local defaultDark = 'neutral'
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
            -- vim.api.nvim_exec_autocmds("Colorscheme", {})
            vim.cmd "hi clear"
            if vim.o.background == "dark" then
                vim.cmd('colorscheme ' .. defaultDark)
            else
                vim.cmd('colorscheme ' .. defaultLight)
            end
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

