#!/usr/bin/env lua
local hi = require "utils/highlights"

local function afterColorscheme()
    -- GitSigns
    local linenr = hi.get("LineNr")['fg']
    local delete = hi.get("DiffDelete")["bg"]
    local stageddelete = hi.get("GitSignsDelete")['fg']
    -- for the sake of consistency with Changedelete
    -- special decides the color for the underline
    hi.mod("GitSignsDelete", {underline=true, special=delete})
    -- bar and underline better than tilde
    hi.mod("GitSignsChangedelete", {underline=true, special=delete})
    hi.set("GitSignsAddNr", {fg=linenr, bg=hi.get("DiffAdd")['fg']})
    hi.set("GitSignsStagedAddNr", {fg=linenr, bg=hi.get("GitSignsStagedAdd")['fg']})
    hi.set("GitSignsChangeNr", {fg=linenr, bg=hi.get("DiffChange")['fg']})
    hi.set("GitSignsStagedChangeNr", {fg=linenr, bg=hi.get("GitSignsChangeAdd")['fg']})
    hi.set("GitSignsChangedeleteNr", {fg=linenr, bg=hi.get("DiffChange")['fg'], underline=true, special=delete})
    hi.set("GitSignsStagedChangedeleteNr", {fg=linenr, bg=hi.get("GitSignsChangedelete")['fg'], underline=true, special=stageddelete})
    -- edge-case where first line(s) of file is deleted
    hi.set("GitSignsTopDeleteNr", {fg=linenr, bg=delete})
    hi.set("GitSignsStagedTopDeleteNr", {fg=linenr, bg=stageddelete})
    hi.set("GitSignsDeleteNr", {fg=linenr, underline=true, special=delete})
    hi.set("GitSignsStagedDeleteNr", {fg=linenr, underline=true, special=stageddelete})

    local spellBad = hi.get("spellBad")
    -- the default already has red undercurl so we don't wanna mess with that but 
    -- other themes sets the fg instead
    if spellBad["special"] then
        spellBad = spellBad["special"]
    else
        spellBad = spellBad["fg"]
    end
    hi.set("SpellBad", {undercurl=true, special=spellBad})

    -- replace whatever search was doing completely but only modify IncSearch's 
    -- reverse setting so it remains a different color from Search, probably 
    -- green.
    hi.rev("Search")
    hi.mod("IncSearch", {reverse=true})
    -- cursor looks to much like the current search word under cursor
    hi.set("CurSearch", {fg="gray", bg=hi.get("IncSearch")["fg"], standout=true, bold=true})
    -- hl.mod("CurSearch", {standout=true})

    -- Default is comment fg and similar strong bg for the whole line. It grabs too much attention.
    -- NonText is bold gray. It feels like a good balance of attention grabbing.
    -- Treesitter hl with NonText "… 35 …" could be too little attention and get overlooked.
    hi.link("Folded", "NonText")

    ---- bufferline
    -- bold instead of italic+bold selected
    hi.set("BufferLineBufferSelected", {bold=true, italic=false})
    hi.set("BufferLineNumbersSelected", {bold=true, italic=false})

    ---- Completion. Some links to IncSearch which is too distracting
    hi.set("CmpItemAbbrMatch", {bold=true})
    hi.set("CmpItemAbbrMatchFuzzy", {bold=true})

    -- there is enough things indicated by color, matched parenthesis etc. with underline is great
    hi.set("MatchParen", {underline=true})

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
    hi.mod("Comment", {italic=false})
    hi.mod("Operator", {bold=true})
    hi.mod("Include", {italic=true})
    hi.mod("Repeat", {italic=true})
    hi.mod("Label", {italic=true})
    hi.mod("Type", {italic=false})
    hi.mod("Conditional", {italic=true, bold=false})
    -- for e.g. julia there is the `a = bool ? b : c` notation. It's weird to 
    -- have ? and : italic since that is meant for words, but it does help 
    -- distinguish them from : used in e.g. ranges.
    -- hl.mod("@conditional.ternary", {italic=false})
    hi.mod("Identifier", {italic=false})
    hi.mod("Number", {italic=false})
    -- vim in lua is @lsp.typemod.variable.global.lua linked to Constant.
    -- If you find other constants that you don't want to make italic then mod the semantic @lsp global instead.
    hi.set("Constant", {italic=true})
    hi.mod("Exception", {italic=true})
    hi.mod("@include", {italic=true})
    hi.link("@keyword", "Keyword")
    hi.mod("Keyword", {italic=true, bold=false})
    hi.mod("@keyword.function", {italic=true})
    hi.mod("@keyword.return", {italic=true})
    hi.mod("@keyword.operator", {italic=true})
    hi.mod("@parameter", {italic=false})
    hi.mod("@function", {italic=false, bold=true})
    hi.mod("@function.call", {bold=false})
    hi.link("@function.macro", "@function.call")
    hi.mod("@conditional", {italic=true, bold=false})
    hi.mod("@conditional.ternary", {italic=false})
    hi.mod("@repeat", {italic=true, bold=false})
    hi.link("@boolean", "Boolean")
    hi.mod("Boolean", {italic=true, bold=false})
    hi.mod("@variable.builtin", {italic=true})
    hi.mod("@function.builtin", {italic=true, bold=false})
    hi.mod("@constant.builtin", {italic=true})
    hi.mod("@type.builtin",     {italic=true})
    hi.mod("@markup.link.url",  {italic=false}) -- underscore is enough distinction
    -- By default links to Keyword which we highlights with both italic and background glow in current theme.
    -- Making documentation text stand out more than normal comments is fine, but not that much, and def not italic, since doc is there to be read as long text.
    hi.fg("@string.documentation", hi.get("Keyword")['fg'])
    -- no color
    hi.link("@markup.italic", "Italic")
    hi.link("@markup.strong", "Bold")

    -- delim. Currently using the default for parentheses.
    hi.link("Delimiter", "RainbowDelimiterViolet")
    hi.link("@punctuation.bracket", "RainbowDelimiterViolet")
    hi.link("@punctuation.delimiter", "Delimiter")

    -- NonText shouldn't be exactly like comments
    if hi.get("Comment")['fg'] == hi.get("NonText")['fg'] then
        hi.mod("NonText", {fg="gray"})
    end
    -- it seems pretty good if they're bold even if the color is similar.
    hi.mod("NonText", {bold=true})

    -- semantic tokens
    hi.mod("@lsp.mod.strong", {bold=true})
    hi.mod("@lsp.mod.emph", {italic=true})
    hi.link("@lsp.type.operator", "@operator")
    hi.link("@lsp.type.keyword", "@keyword")
    -- instead of to @function since we only want function definitions to be 
    -- bold and @type.function is often not, e.g. in this very file.
    hi.link("@lsp.type.function", "@function.call")
    hi.link("@lsp.type.method", "@function.call")
    hi.link("@lsp.type.string", "@string")

    --- Colon before function call is captured by @constant making it italic.
    hi.clear("@constant.lua")
    -- remap things that were mapped to @constant to preserve their functioning highlights
    hi.link("@lsp.typemod.variable.global.lua", "Constant")
    hi.link("@lsp.typemod.variable.defaultLibrary.lua", "Constant")
    -- wrong highlight by treesitter
    hi.set("@type.sql", {})
    -- default links to constant which we make italic so no thanks.
    -- @string.special -> Special -> same fg as function call. This works well since @string.special is for e.g. `` cmd in julia.
    hi.link("@string.special", "Special")
    -- symbols in julia are differentiated clearly enough by a preceding colon colored according to delimiters.
    -- We clear its default link here to Constant, since we don't want it italic.
    hi.clear("@string.special.symbol.julia")
    -- By default it's linked to Keyword which gives it builtin color and italic.
    -- Highlight like regular function is clear enough that it's macro, since it ends in !.
    -- We link to @function.call, not @function since we use the latter for function definition, which are in bold.
    hi.link("@lsp.type.macro.rust", "@function.call")

    -- variable is the default capture for most things in code so we want it to 
    -- be neutral
    -- I set it to copy normal instead of using clear since if I clear it will 
    -- be overrideen by other colors, but I want it to appear. Example in julia:
    -- "$variable" will color variable as string with clear and as normal with this approach.
    hi.fg("@variable", hi.get("Normal")["fg"])
    -- Except for some languages where variable shouldn't be neutral:
    for _, showVar in ipairs{"typst", "markdown", "asciidoc", "latex"} do
        hi.fg("@variable." .. showVar, hi.get("@variable.builtin")["fg"])
    end
    -- wrong double annotation as variable for functions
    hi.clear("@variable.wgsl")

    -- extmarks
    -- By default colors. Underline variants makes more sense.
    hi.set("DiagnosticUnderlineHint", {underdotted=true})
    hi.set("DiagnosticUnderlineInfo", {underdotted=true})
    -- subtle. underline and underdashed are stronger but the warn is often 
    -- wrong, e.g. missing reference to things LSP doesn't understand is 
    -- imported.
    hi.set("DiagnosticUnderlineWarn", {underdotted=true})
    hi.set("DiagnosticUnderlineError", {undercurl=true})

    -- Make a hl group we can link to that hides text
    hi.hide("Background")
end

-- local defaultDark = 'fluoromachine'
-- local defaultDark = 'delta'
local defaultDark = 'carbonfox'
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

