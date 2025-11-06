local hi = require "utils/highlights"

local function _basic_styling()
    -- Assume we will be working in terminal that supports underline, bold, italic, etc.
    hi.set("@text.underline", { underline = true })
    hi.link("@markup.underline", "@text.underline")
    hi.link("@text.strong", "Bold")
    -- Remove color and use the relevant styles.
    hi.link("@markup.italic", "Italic")
    hi.link("@markup.strong", "Bold")
    hi.mod("@lsp.mod.strong", { bold = true })
    hi.mod("@lsp.mod.emph", { italic = true })
end

local function _GitSigns()
    local linenr = hi.fg("LineNr")
    local delete = hi.bg("DiffDelete")
    local stageddelete = hi.fg("GitSignsDelete")
    -- For the sake of consistency with Changedelete.
    -- Special decides the color for the underline.
    hi.mod("GitSignsDelete", { underline = true, special = delete })
    -- Bar and underline better than tilde.
    hi.mod("GitSignsChangedelete", { underline = true, special = delete })
    hi.set("GitSignsAddNr", { fg = linenr, bg = hi.fg("DiffAdd") })
    hi.set("GitSignsStagedAddNr", { fg = linenr, bg = hi.fg("GitSignsStagedAdd") })
    hi.set("GitSignsChangeNr", { fg = linenr, bg = hi.fg("DiffChange") })
    hi.set("GitSignsStagedChangeNr", { fg = linenr, bg = hi.fg("GitSignsChangeAdd") })
    hi.set("GitSignsChangedeleteNr", { fg = linenr, bg = hi.fg("DiffChange"), underline = true, special = delete })
    hi.set("GitSignsStagedChangedeleteNr",
        { fg = linenr, bg = hi.fg("GitSignsChangedelete"), underline = true, special = stageddelete })
    -- Edge-case where first line(s) of file is deleted.
    hi.set("GitSignsTopDeleteNr", { fg = linenr, bg = delete })
    hi.set("GitSignsStagedTopDeleteNr", { fg = linenr, bg = stageddelete })
    hi.set("GitSignsDeleteNr", { fg = linenr, underline = true, special = delete })
    hi.set("GitSignsStagedDeleteNr", { fg = linenr, underline = true, special = stageddelete })
end

local function _diagnostics()
    -- By default indicated by colors.
    -- Underline etc. makes more sense when supported.
    hi.set("DiagnosticUnderlineHint", { underdotted = true, special = hi.fg("DiagnosticHint") })
    hi.set("DiagnosticUnderlineInfo", { underdotted = true, special = hi.fg("DiagnosticInfo") })
    -- subtle. underline and underdashed are stronger but the warn is often
    -- wrong, e.g. missing reference to things LSP doesn't understand is
    -- imported.
    hi.set("DiagnosticUnderlineWarn", { underdotted = true, special = hi.fg("DiagnosticWarn") })
    hi.set("DiagnosticUnderlineError", { undercurl = true, special = hi.fg("DiagnosticError") })
    -- instead of mildly red text, do red underline.
    hi.set("Error", { undercurl = true, special = hi.fg("ErrorMsg") })
end

local function _spellBad()
    -- Some colorschemes sets the fg for spellBad. We make sure the colour
    -- chosen is used, but always as undercurl.
    local spellBad = hi.get("spellBad")
    spellBad = spellBad["special"] or spellBad["fg"]
    hi.set("SpellBad", { undercurl = true, special = spellBad })
end

---Function to apply after setting a colorscheme so that we can use any
---colorscheme without customising it, and have certain rules always apply.
local function afterColorscheme()
    -- Aim to make this function organised from things that should be default,
    -- to the most opinionated at the bottom.
    _basic_styling()
    _diagnostics()
    _spellBad()

    -- TEMP: the colour for statement should look more like function. Trying something closer, e.g. blue for function, purple for statement.
    local statement_fg = "#8577af"
    hi.mod("Statement", {fg=statement_fg})
    hi.mod("Keyword", {fg=statement_fg})
    -- Maybe make this slightly different in future if we think it's meaningful
    hi.mod("Conditional", {fg=statement_fg})

    hi.mod("NonText", { bold = true })
    -- NonText shouldn't be exactly like comments
    if hi.fg("Comment") == hi.fg("NonText") then
        hi.mod("NonText", { fg = "gray" })
    end
    -- MoreMsg is shown for text in multiline messages, e.g. :Inspect.
    -- By default it is a bold blue color like a function def which is confusing.
    -- NonText link makes more sense, it's not text that can be edited.
    hi.link("MoreMsg", "NonText")
    -- Default is comment fg and similar strong bg for the whole line. It grabs too much attention.
    -- NonText is bold gray. It feels like a good balance of attention grabbing.
    -- Treesitter hl with NonText "… 35 …" could be too little attention and get overlooked.
    hi.link("Folded", "NonText")

    -- Background colour for cmdline.
    hi.setbg("MsgArea", hi.bg("StatusLine"))

    hi.mod("DiffDelete", {fg=hi.bg("Normal")})
    _GitSigns()

    -- Matching parenthesis is indicated by colour sometimes.
    -- Reverse is common but takes a lot of focus.
    -- Using underline is more subtle and allows the original colors to stay.
    hi.set("MatchParen", { underline = true })
    hi.set("IncSearch", { standout = true, bold = true })
    hi.set("Search", { fg = "gray", bg = hi.fg("IncSearch"), standout = true })
    hi.link("CurSearch", "Search")

    -- TEMP, fix using lush plugin
    -- link to function fg colour
    hi.setfg("function.call", "#73a3b7")
    -- dimmed down version of @import / Include / PreProc. Use darkening with lush in dark mode and lighten in light mode.
    hi.setfg("@module", "#5e5050")

    -- Definitions are bold, while the subsequent usage of these classes, types, functions etc are not.
    -- There's also @lsp.mod.definition, which is used when defining e.g. arguments to a function.
    hi.mod("@lsp.typemod.class.definition", { bold = true })
    hi.mod("@lsp.typemod.method.definition", { bold = true })
    hi.set("@lsp.typemod.function.declaration", { bold = true, fg = hi.fg("Function") })

    ---- cmp completion.
    -- Some colorschemes links to IncSearch which is too distracting
    hi.set("CmpItemAbbrMatch", { bold = true })
    hi.set("CmpItemAbbrMatchFuzzy", { bold = true })

    -- Curly brackets in lua.
    hi.setfg("@constructor.lua", hi.fg("Macro"))

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
    hi.set("Operator", { bold = true, fg = hi.fg("Function") })
    hi.set("Include", { italic = true, bold=true, fg=statement_fg })
    hi.mod("Repeat", { italic = true })
    -- Labels are meant to be read, they should definitely not be italic.
    hi.mod("Label", { italic = false, bold = true })
    hi.mod("Type", { italic = false })
    hi.mod("Conditional", { italic = true, bold = false })
    -- for e.g. julia there is the `a = bool ? b : c` notation. It's weird to
    -- have ? and : italic since that is meant for words, but it does help
    -- distinguish them from : used in e.g. ranges.
    -- We can also replace italics with bold, similar to operator, but it's still different because of the colour.
    hi.mod("@keyword.conditional.ternary", {italic=false, bold=true})
    hi.mod("Identifier", { italic = false })
    -- TEMP from blending string colour with the red-ish for more builtin things.
    local value_fg = "#7eb1b1"
    hi.mod("Number", { italic = false, fg= value_fg})
    -- If you find constants that you don't want to make italic then mod the semantic @lsp global instead.
    hi.set("Constant", { italic = false, fg=value_fg })
    hi.set("@constant", { italic = false, fg=value_fg })
    -- Nothing is constant in python and it's just based on if chars are uppercase.
    hi.set("@constant.python", {})
    -- In python format string "{x:.3f}" the .3f is captured as @none.python, which is ignored by default showing .3f as string.
    hi.link("@none.python", "Macro")
    hi.mod("Exception", { italic = true })
    hi.mod("@include", { italic = true })
    hi.link("@keyword", "Keyword")
    hi.mod("Keyword", { italic = true, bold = false })
    hi.mod("@keyword.function", { italic = true, fg=statement_fg })
    hi.set("@keyword.return", { italic = true, bold=true, fg=statement_fg })
    -- all same as keyword except bold since operators are bold.
    hi.set("@keyword.operator", { italic = true, bold = true, fg = hi.fg("@operator") })
    -- TEMP hardcoded colour.
    hi.mod("@parameter", { italic = false, fg = "#ac9ba1" })
    -- attribute is by default linked to constant.
    -- Attributes are not constant.
    -- TODO: what's the difference between attribute, parameter, and property?
    hi.link("@attribute", "@parameter")
    hi.link("@variable.attribute", "@attribute")
    hi.mod("@function", { italic = false, bold = true })
    hi.setfg("@function.call", hi.fg("Function"))
    hi.link("@function.method", "function.call")
    hi.link("@function.method.call", "function.call")
    hi.link("@function.macro", "@function.call")
    hi.mod("@conditional", { italic = true, bold = false })
    hi.mod("@conditional.ternary", { italic = false })
    hi.mod("@repeat", { italic = true, bold = false })
    hi.link("@boolean", "Boolean")
    hi.mod("Boolean", { italic = true, bold = false })
    -- italic is enough distinction, and fits with the pattern. No need for a different colour.
    hi.set("@variable.builtin", { italic = true, })
    hi.mod("@variable.attribute.builtin", { italic = true })
    hi.set("@function.builtin", { italic = true, fg = hi.fg("@function.call") })
    -- hi.link("@attribute", "PreProc")
    hi.set("@attribute.builtin", { italic = true, fg = hi.fg("@attribute") })
    hi.set("@constant.builtin", { bold=false, italic = true, fg = value_fg }) -- The value_fg is enough to indicate constant.
    hi.set("@type.builtin", { italic = true, fg = hi.fg("@type") })
    hi.set("@identifier.builtin", { italic = true, fg = hi.fg("Identifier") })
    -- @lsp understands types better than TS. TS annotates def type(...) in class as @type.builtin.
    hi.set("@type.builtin.python", { italic = false })
    hi.mod("@module.builtin", { italic = true, fg = hi.fg("@module") })
    -- underline tags since they are kinda like links in some filetypes like asciidoc I think, but they are more basic elements in xml and html.
    hi.mod("Tag", { underline = true })
    hi.link("@tag", "Tag")
    hi.mod("@tag.html", {underline=false})
    hi.mod("@tag.xml", {underline=false})
    hi.mod("@tag.builtin", { underline = true, italic = true })
    -- E.g. the "version" in xml `<release version="...`, i.e. it's not both a tag and an attribute, it's an attribute of a tag.
    hi.set("@tag.attribute", { underline = false, italic = false, fg=hi.fg("@attribute") })
    hi.link("@tag.delimiter", "@comment")
    hi.mod("@lsp.type.ref", { underline = true })
    hi.mod("@lsp.type.link", { underline = true })
    -- E.g. in typst with \<
    hi.link("@lsp.type.escape", "@string.escape")
    hi.set("@markup.raw", { underline = false, fg = hi.fg("@markup.raw") }) -- trying to just remove italic
    hi.set("@markup.link", { underline = true, fg = hi.fg("@markup.link.label") })
    hi.mod("@markup.link.url", { italic = false })                          -- underscore is enough distinction
    hi.mod("@markup.link.url", { italic = false })                          -- underscore is enough distinction
    hi.mod("@string.special.url", { italic = false })                       -- underscore is enough distinction

    -- I like having @string.documentation different colour from regular string to make it clear it has a different special role and is recognised as such.
    -- By default it was linked to keyword which is implying builtin, e.g. italic.
    -- TEMP: between string, comment, function, and statement/keyword colours
    hi.set("@string.documentation", {fg="#738fa6"})

    -- delim.
    hi.setfg("Delimiter", hi.fg("Keyword"))
    hi.link("@punctuation.bracket", "Delimiter")
    hi.link("@punctuation.delimiter", "Delimiter")
    -- In python this is curly braces in f"{...}" which are like delimiters, 
    -- except they are neutral. They don't form a dict like regular {} and the 
    -- different highlighting within them vs the string around makes them 
    -- visually redundant. For this reason we want them to look like delimiter but dimmed.
    -- TEMP: this is hardcoded here.
    hi.set("@punctuation.special", {fg="#504768"})
    -- Was overwriting the rainbow ext marks:
    hi.clear("@lsp.type.punct.typst")
    -- Not sure what "pol" is but it was lined to @variable which is neutral color globally but not for typst.
    hi.link("@lsp.type.pol.typst", "@variable.typst")

    -- semantic tokens
    hi.link("@lsp.type.operator", "@operator")
    hi.link("@lsp.type.keyword", "@keyword")
    -- E.g. `Self` in `def fn() -> Self` in python. Might have to remove if it also covers non-builtins.
    -- hi.link("@lsp.type.typeparameter", "@type.builtin")
    -- instead of to @function since we only want function definitions to be
    -- bold and @type.function is often not, e.g. in this very file.
    hi.link("@lsp.type.function", "@function.call")
    -- Python has it on e.g. `@dataclass(...)` before `class`, which is not Function coloured if called like `@dataclass`.
    -- However it will also have this hi on `return *fn(...)` so we should keep it for that.
    -- hi.clear("@lsp.type.function.python")
    hi.link("@lsp.type.method", "@function.call")
    hi.link("@lsp.type.string", "@string")
    -- E.g. in importing macro in rust. By default no colour.
    hi.link("@lsp.type.procMacro", "@function.macro")
    hi.link("@function.typst", "@function.call")

    --- Colon before function call is captured by @constant making it italic.
    hi.clear("@constant.lua")
    hi.link("@constant.numeric", "Number") -- typst. Remove italic.
    -- remap things that were mapped to @constant to preserve their functioning highlights
    hi.link("@lsp.typemod.variable.global.lua", "Constant")
    hi.link("@lsp.typemod.variable.defaultLibrary.lua", "Constant")
    -- wrong highlight by treesitter
    hi.clear("@type.sql")
    -- default links to constant which we make italic so no thanks.
    -- @string.special -> Special -> same fg as function call. This works well since @string.special is for e.g. `` cmd in julia.
    -- @punctuation.special can also be used, it has a similar colour in terafox.
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
    hi.set("Title", { bold = true, underdouble = true, fg = hi.fg("Normal") })

    -- variable is the default capture for most things in code so we want it to
    -- be neutral. I set it to copy normal instead of using clear since if I clear it will
    -- be overrideen by other colors, but I want it to appear. Example in julia:
    -- "$variable" will color variable as string with clear and as normal with this approach.
    -- However, there are situations where the opposite is true, e.g. in
    -- fortran some types will also be captured as variables, and we would ideally
    -- keep the type colouring.
    -- local variable = hi.getfg("@variable")
    local variable = hi.fg("@attribute")
    hi.setfg("@variable", hi.fg("Normal"))
    -- Except for some languages where variable shouldn't be neutral:
    for _, showVar in ipairs { "typst", "markdown", "asciidoc", "latex" } do
        hi.setfg("@variable." .. showVar, variable)
    end
    -- wrong double annotation as variable for functions
    hi.clear("@variable.wgsl")

    -- we want members and properties to be neutral colour as well
    hi.link("@variable.member", "@variable")
    hi.link("@property", "@variable")


    -- remove title link, which makes it bold
    hi.link("FidgetTitle", "Normal")

    -- default was comment but it's stuff that can be more important, not less
    hi.link("helpCommand", "@string.special")
    hi.link("helpExample", "@string.documentation")

    -- vim in lua is @lsp.typemod.variable.global.lua linked to Constant. We link to italic module instead to get italic and faded colour.
    hi.link("@lsp.typemod.variable.global", "@module.builtin")
    hi.link("@lsp.typemod.variable.global.lua", "@module.builtin")
    hi.link("@lsp.type.selfParameter", "@module.builtin") -- The word self in a class.
    -- wrongfully overrides delimiter ':' in lua, e.g.:
    -- vim.opt_local.formatoptions:remove('o')
    hi.clear("@lsp.type.variable.lua")
    -- In python some variables are highlighted as types because they start with an uppercase letter.
    -- It seems we can just clear that highlight since there is also the @lsp.type.class.python which highlights the types.
    hi.clear("@type.python")

    -- By default links to Error and gets red undercurl.
    hi.set("ConflictMarkerBegin", { bg = "#0f3625" })
    hi.link("ConflictMarkerOurs", "ConflictMarkerBegin")
    hi.set("ConflictMarkerCommonAncestors", { bg = "#1b4849" })
    hi.link("ConflictMarkerCommonAncestorsHunk", "ConflictMarkerCommonAncestors")
    hi.link("ConflictMarkerSeparator", "Normal") -- Remove red undercurl
    hi.set("ConflictMarkerEnd", { bg = "#36324e" })
    hi.link("ConflictMarkerTheirs", "ConflictMarkerEnd")

    -- shell
    -- vim doesn't seem to be combining underline and color from two separate
    -- groups? Maybe that's only for TS which we don't have working for shell.
    -- We do this hack solution of a new hl group which has the underline and the colour.
    hi.set("FunctionPath", { underline = true, fg = hi.fg("Function") })
    -- Same hack:
    hi.set("zshShortDerefPath", { underline = true, fg = hi.fg("zshShortDeref") })
    hi.set("Wildcard", hi.get("@operator"))
    hi.mod("Wildcard", { underline = true })

    -- Normally linked to title, which is fair enough but I have double
    -- underline for titles, which makes the dashboard (nvim opened without any
    -- file etc) look weird.
    hi.link("DashboardHeader", "Bold")
    hi.link("@lsp.type.decorator", "PreProc")

    -- Only relevant to json without treesitter
    -- hi.link("jsonKeyword", "@variable")
    -- hi.link("jsonQuote", "String")
    -- for e.g. '"' around string in dockerfile.
    -- TODO: used a dimmed string color?
    hi.link("Quote", "String")
    -- -- Why are commas called noise?
    -- hi.link("jsonNoise", "Delimiter")
    -- Dim things that would be concealed if conceallevel>0
    hi.link("@conceal", "comment")
    -- The quotation marks are technically indicating strings
    hi.link("@conceal.json", "String")

    -- underline active parameter in signature help rather than colour it in some pale unhelpful colour.
    -- blink.cmp should be defaulting to this hl with its BlinkCmpSignatureHelpActiveParameter.
    hi.link("LspSignatureActiveParameter", "Underlined")

    -- local StatusLineNC = hi.get("StatusLineNC")
    -- hi.set("StatusLine", {
    --     -- underline = true,
    --     bold = true,
    --     fg = StatusLineNC["fg"],
    --     bg = StatusLineNC["bg"],
    -- })
    -- For reformatting statusline as a simple border with laststatus=0.
    -- See lua/plugins/UI.lua
    hi.set('Statusline', { bg = nil, fg = hi.bg("Statusline") })
    hi.set('StatuslineNC', { bg = nil, fg = hi.bg("StatuslineNC") })

    hi.mod("TabLine", { fg = "grey" })
    hi.mod("TabLineFill", { fg = "grey" })

    -- Indicated for the whole regex string. By default shown as a colour like `Type` has, which is just confusing.
    -- Should be same as String or at least similar.
    hi.link("@string.regexp", "@string")
    -- @string.escape.regex linked to this by default.
    -- Could be linked to Character in the future when this is different from string.
    hi.link("@string.escape", "Macro")
    -- E.g. "\." for literal dot in regexp.
    hi.link("@string.regexp.regex", "@string.escape")
    -- Compromise between Operator and other things?
    -- hi.set("@operator.regex", {fg=hi.fg("Macro"), bold=true})

    -- Make a prompt hl group that matches shell, REPL, etc. prompts that various prompts in nvim can link to.
    hi.setfg("Prompt", "#FF6AC1")
    hi.link("MiniPickPrompt", "Prompt")
    hi.link("SnacksPickerPrompt", "Prompt")
    hi.link("FzfLuaFzfPrompt", "Prompt")
    -- hi.link("FzfLuaLivePrompt", "Prompt")
    hi.link("TelescopeSelectionCaret", "Prompt")
end

-- local defaultDark = 'fluoromachine'
-- local defaultDark = 'delta'
local defaultDark = 'terafox'
-- local defaultDark = 'neutral'
-- local defaultLight = 'kanagawa-lotus'
local defaultLight = 'dawnfox'

-- Don't clear so we can add to the same group for plugin and filetype specific hls.
local grp = vim.api.nvim_create_augroup("afterColorscheme", { clear = false })

vim.api.nvim_create_autocmd("ColorScheme", {
    group = grp,
    callback = afterColorscheme,
})

vim.api.nvim_create_autocmd("VimEnter", {
    group = grp,
    callback = function()
        vim.schedule(function()
            -- pretend we called colorscheme in order to trigger all autcmds
            -- that fire after setting a new colorscheme.
            -- vim.api.nvim_exec_autocmds("Colorscheme", {})
            -- vim.cmd "hi clear"
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
