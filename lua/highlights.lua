#!/usr/bin/env lua
local hi = require "utils/highlights"

local function afterColorscheme()
    -- GitSigns
    local linenr = hi.get("LineNr")['fg']
    local delete = hi.get("DiffDelete")['bg']
    local stageddelete = hi.get("GitSignsDelete")['fg']
    -- for the sake of consistency with Changedelete
    -- special decides the color for the underline
    hi.mod("GitSignsDelete", { underline = true, special = delete })
    -- bar and underline better than tilde
    hi.mod("GitSignsChangedelete", { underline = true, special = delete })
    hi.set("GitSignsAddNr", { fg = linenr, bg = hi.get("DiffAdd")['fg'] })
    hi.set("GitSignsStagedAddNr", { fg = linenr, bg = hi.get("GitSignsStagedAdd")['fg'] })
    hi.set("GitSignsChangeNr", { fg = linenr, bg = hi.get("DiffChange")['fg'] })
    hi.set("GitSignsStagedChangeNr", { fg = linenr, bg = hi.get("GitSignsChangeAdd")['fg'] })
    hi.set("GitSignsChangedeleteNr", { fg = linenr, bg = hi.get("DiffChange")['fg'], underline = true, special = delete })
    hi.set("GitSignsStagedChangedeleteNr",
        { fg = linenr, bg = hi.get("GitSignsChangedelete")['fg'], underline = true, special = stageddelete })
    -- edge-case where first line(s) of file is deleted
    hi.set("GitSignsTopDeleteNr", { fg = linenr, bg = delete })
    hi.set("GitSignsStagedTopDeleteNr", { fg = linenr, bg = stageddelete })
    hi.set("GitSignsDeleteNr", { fg = linenr, underline = true, special = delete })
    hi.set("GitSignsStagedDeleteNr", { fg = linenr, underline = true, special = stageddelete })

    local spellBad = hi.get("spellBad")
    -- the default already has red undercurl so we don't wanna mess with that but
    -- other themes sets the fg instead
    if spellBad["special"] then
        spellBad = spellBad["special"]
    else
        spellBad = spellBad["fg"]
    end
    hi.set("SpellBad", { undercurl = true, special = spellBad })

    hi.set("IncSearch", { standout = true, bold = true })
    hi.set("Search", { fg = "gray", bg = hi.get("IncSearch")["fg"], standout = true })
    hi.link("CurSearch", "Search")

    -- Default is comment fg and similar strong bg for the whole line. It grabs too much attention.
    -- NonText is bold gray. It feels like a good balance of attention grabbing.
    -- Treesitter hl with NonText "… 35 …" could be too little attention and get overlooked.
    hi.link("Folded", "NonText")

    ---- bufferline
    -- bold instead of italic+bold selected.
    -- Grey instead of white from using the fg from unselected tab colour ("visible")
    hi.set("BufferLineBufferSelected",
        {
            bold = true,
            italic = false,
            bg = hi.get("BufferLineBufferSelected")['bg'],
            fg = hi.getfg(
                "BufferLineBufferVisible")
        })
    hi.set("BufferLineNumbersSelected",
        {
            bold = false,
            italic = false,
            bg = hi.get("BufferLineNumbersSelected")['bg'],
            fg = hi.getfg(
                "BufferLineBufferVisible")
        })

    ---- Completion. Some links to IncSearch which is too distracting
    hi.set("CmpItemAbbrMatch", { bold = true })
    hi.set("CmpItemAbbrMatchFuzzy", { bold = true })

    -- there is enough things indicated by color, matched parenthesis etc. with underline is great
    hi.set("MatchParen", { underline = true })

    -- Should be default
    hi.set("@text.underline", { underline = true })
    hi.link("@markup.underline", "@text.underline")
    hi.set("@text.strong", { bold = true })

    -- temp, fix using lush plugin
    hi.fg("function.call", "#73a3b7") -- link to function fg colour
    hi.fg("@module", "#5e5050")       -- dimmed down version of @import / Include / PreProc. Use darkening with lush in dark mode and lighten in light mode.
    hi.fg("@constructor", hi.getfg("TSRainbowBlue"))
    hi.link("@constructor.lua", "@constructor")

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
    hi.mod("Comment", { italic = false })
    hi.set("Operator", { bold = true, fg = hi.get("Keyword")['fg'] })
    hi.mod("Include", { italic = true })
    hi.mod("Repeat", { italic = true })
    hi.mod("Label", { italic = true })
    hi.mod("Type", { italic = false })
    hi.mod("Conditional", { italic = true, bold = false })
    -- for e.g. julia there is the `a = bool ? b : c` notation. It's weird to
    -- have ? and : italic since that is meant for words, but it does help
    -- distinguish them from : used in e.g. ranges.
    -- hl.mod("@conditional.ternary", {italic=false})
    hi.mod("Identifier", { italic = false })
    hi.mod("Number", { italic = false })
    -- If you find constants that you don't want to make italic then mod the semantic @lsp global instead.
    hi.set("Constant", { italic = true })
    hi.mod("Exception", { italic = true })
    hi.mod("@include", { italic = true })
    hi.link("@keyword", "Keyword")
    hi.mod("Keyword", { italic = true, bold = false })
    hi.mod("@keyword.function", { italic = true })
    hi.mod("@keyword.return", { italic = true })
    -- all same as  keyword except bold since operators are bold.
    hi.set("@keyword.operator", { italic = true, bold = true, fg = hi.get("Keyword")['fg'] })
    hi.mod("@parameter", { italic = false })
    hi.mod("@function", { italic = false, bold = true })
    hi.fg("@function.call", hi.getfg("Function"))
    hi.link("@function.method", "function.call")
    hi.link("@function.method.call", "function.call")
    hi.link("@function.macro", "@function.call")
    hi.mod("@conditional", { italic = true, bold = false })
    hi.mod("@conditional.ternary", { italic = false })
    hi.mod("@repeat", { italic = true, bold = false })
    hi.link("@boolean", "Boolean")
    hi.mod("Boolean", { italic = true, bold = false })
    -- italic is enough distinction, and fits with the pattern. No need for a different colour.
    hi.set("@variable.builtin", { italic = true, fg = hi.get("@variable")['fg'] })
    hi.mod("@variable.parameter.builtin", { italic = true })
    hi.set("@function.builtin", { italic = true, fg = hi.get("@function.call")['fg'] })
    hi.set("@attribute.builtin", { italic = true, fg = hi.get("@attribute")['fg'] })
    hi.mod("@constant.builtin", { italic = true, fg = hi.get("@constant")['fg'] })
    hi.set("@type.builtin", { italic = true, fg = hi.get("@type")['fg'] })
    hi.mod("@module.builtin", { italic = true, fg = hi.get("@module")['fg'] })
    -- underline tags since they are kinda like links
    hi.mod("Tag", { underline = true })
    hi.mod("@tag", { underline = true })
    hi.mod("@tag.builtin", { underline = true, italic = true })
    hi.mod("@tag.attribute", { underline = true })
    hi.mod("@tag.delimiter", { underline = true })
    hi.mod("@lsp.type.ref", { underline = true })
    hi.mod("@lsp.type.link", { underline = true })
    hi.set("@markup.raw", { underline = false, fg = hi.get("@markup.raw")['fg'] }) -- trying to just remove italic
    hi.set("@markup.link", { underline = true, fg = hi.get("@markup.link.label")['fg'] })
    hi.mod("@markup.link.url", { italic = false })                                 -- underscore is enough distinction
    hi.mod("@markup.link.url", { italic = false })                                 -- underscore is enough distinction
    hi.mod("@string.special.url", { italic = false })                              -- underscore is enough distinction
    -- By default links to Keyword which we highlights with both italic and background glow in current theme.
    -- Making documentation text stand out more than normal comments is fine, but not that much, and def not italic, since doc is there to be read as long text.
    hi.fg("@string.documentation", hi.get("Keyword")['fg'])
    -- no color
    hi.link("@markup.italic", "Italic")
    hi.link("@markup.strong", "Bold")

    -- delim. Currently using the default for parentheses.
    hi.link("Delimiter", "TSRainbowViolet") -- or RainbowDelimitersViolet
    hi.link("@punctuation.bracket", "Delimiter")
    hi.link("@punctuation.delimiter", "Delimiter")
    -- Was overwriting the rainbow ext marks:
    hi.clear("@lsp.type.punct.typst")
    -- Not sure what "pol" is but it was lined to @variable which is neutral color globally but not for typst.
    hi.link("@lsp.type.pol.typst", "@variable.typst")

    -- NonText shouldn't be exactly like comments
    if hi.get("Comment")['fg'] == hi.get("NonText")['fg'] then
        hi.mod("NonText", { fg = "gray" })
    end
    -- it seems pretty good if they're bold even if the color is similar.
    hi.mod("NonText", { bold = true })

    -- semantic tokens
    hi.mod("@lsp.mod.strong", { bold = true })
    hi.mod("@lsp.mod.emph", { italic = true })
    hi.link("@lsp.type.operator", "@operator")
    hi.link("@lsp.type.keyword", "@keyword")
    -- instead of to @function since we only want function definitions to be
    -- bold and @type.function is often not, e.g. in this very file.
    hi.link("@lsp.type.function", "@function.call")
    hi.link("@lsp.type.method", "@function.call")
    hi.link("@lsp.type.string", "@string")
    hi.link("@function.typst", "@function.call")

    --- Colon before function call is captured by @constant making it italic.
    hi.clear("@constant.lua")
    hi.link("@constant.numeric", "Number") -- typst. Remove italic.
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

    -- fixes typst highlighting operator in == heading
    hi.link("@text.title", "Title")
    hi.set("Title", { bold = true, underdouble = true, fg = hi.get("Normal")['fg'] })

    -- variable is the default capture for most things in code so we want it to
    -- be neutral. I set it to copy normal instead of using clear since if I clear it will
    -- be overrideen by other colors, but I want it to appear. Example in julia:
    -- "$variable" will color variable as string with clear and as normal with this approach.
    -- local variable = hi.get("@variable")["fg"]
    local variable = hi.get("@parameter")["fg"]
    hi.fg("@variable", hi.get("Normal")["fg"])
    -- Except for some languages where variable shouldn't be neutral:
    for _, showVar in ipairs { "typst", "markdown", "asciidoc", "latex" } do
        hi.fg("@variable." .. showVar, variable)
    end
    -- wrong double annotation as variable for functions
    hi.clear("@variable.wgsl")

    -- we want members and properties to be neutral colour as well
    hi.link("@variable.member", "@variable")
    hi.link("@property", "@variable")

    -- extmarks
    -- By default colors. Underline variants makes more sense.
    hi.set("DiagnosticUnderlineHint", { underdotted = true })
    hi.set("DiagnosticUnderlineInfo", { underdotted = true })
    -- subtle. underline and underdashed are stronger but the warn is often
    -- wrong, e.g. missing reference to things LSP doesn't understand is
    -- imported.
    hi.set("DiagnosticUnderlineWarn", { underdotted = true })
    hi.set("DiagnosticUnderlineError", { undercurl = true })
    -- instead of mildly red text, do red underline.
    hi.set("Error", { undercurl = true, special = "red" })

    -- Make a hl group we can link to that hides text
    hi.hide("Background")

    -- remove title link, which makes it bold
    hi.link("FidgetTitle", "Normal")

    -- default was comment but it's stuff that can be more important, not less
    hi.link("helpCommand", "@string.special")
    hi.link("helpExample", "@string.documentation")

    -- vim in lua is @lsp.typemod.variable.global.lua linked to Constant. We link to italic module instead to get italic and faded colour.
    hi.link("@lsp.typemod.variable.global", "@module.builtin")
    hi.link("@lsp.typemod.variable.global.lua", "@module.builtin")

    -- By default links to Error and gets red undercurl.
    hi.set("ConflictMarkerBegin", { bg = "#0f3625" })
    hi.link("ConflictMarkerOurs", "ConflictMarkerBegin")
    hi.set("ConflictMarkerCommonAncestors", { bg = "#1b4849" })
    hi.link("ConflictMarkerCommonAncestorsHunk", "ConflictMarkerCommonAncestors")
    hi.link("ConflictMarkerSeparator", "Normal") -- Remove red undercurl
    hi.set("ConflictMarkerEnd", { bg = "#36324e" })
    hi.link("ConflictMarkerTheirs", "ConflictMarkerEnd")
end

-- local defaultDark = 'fluoromachine'
-- local defaultDark = 'delta'
local defaultDark = 'terafox'
-- local defaultDark = 'neutral'
-- local defaultLight = 'kanagawa-lotus'
local defaultLight = 'dawnfox'

local grp = vim.api.nvim_create_augroup("afterColorscheme", { clear = true })

vim.api.nvim_create_autocmd("Colorscheme", {
    pattern = "*",
    group = grp,
    callback = afterColorscheme,
})

vim.api.nvim_create_autocmd("VimEnter", {
    pattern = "*",
    group = grp,
    callback = function()
        vim.schedule(function()
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
vim.api.nvim_create_user_command("Dark", function()
    vim.fn.system("~/dotfiles/dark.sh")
    vim.cmd('colorscheme ' .. defaultDark)
end, {})
vim.api.nvim_create_user_command("Light", function()
    vim.fn.system("~/dotfiles/light.sh")
    vim.cmd('colorscheme ' .. defaultLight)
end, {})
