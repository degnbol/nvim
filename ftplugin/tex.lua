-- NOTE: only currently attached to the first tex file opened (since require runs things once).
require "tex.overleaf"
local tbl = require "tex.tables"
require "tex.cmds".map_keys()
require "tex.textcolor".color_textcolor()
local hi = require "utils/highlights"
local util = require "utils/init" -- string.contains and schedule_notify
local latexmk = require "tex.latexmk"
local map = require "utils/keymap"

-- also defined in lua/plugins/mini.lua as shift+` but we can just use ` since
-- I don't think it has use in latex, maybe except in some verbatim code block or something?
require 'mini.surround'.config.custom_surroundings['`'] = {
    input = { "``().-()''" },
    output = { left = '``', right = "''" },
}

-- requires kana/vim-textobj-user, see lua/plugins/init.lua
vim.fn["textobj#user#plugin"]("tex", {
    ['latex-ticks'] = {
        pattern = { '``', "''" },
        ['select-a'] = '<buffer> a`',
        ['select-i'] = '<buffer> i`',
    },
})

-- keymaps assume vimtex is used.

map.desc('n', '<LocalLeader>a', "Context menu")
map.n('<LocalLeader>A', "<plug>TableAlign", "Align table", {buffer=true})
map.n("<LocalLeader>y", "<plug>TableYank", "Yank as TSV", {buffer=true})
map.nx("<LocalLeader>p", "<plug>TablePaste", "Paste TSV", {buffer=true})
map.n('da|', "<Plug>TableDelCol", "Delete a table column", {buffer=true})
map.n("dix", tbl.deleteInCell, "Delete in cell", { silent = true, buffer=true })
map.n("cix", tbl.changeInCell, "Change in cell", { silent = true, buffer=true })
map.n("xix", function() tbl.deleteInCell("+") end, "Cut in cell", { silent = true , buffer=true})
map.o("ix", "<Plug>TableSelInCell", "Select in cell", { silent = true, })
map.n('<LocalLeader>[', "<Plug>TableSwapLeft", "Swap table column left", {buffer=true})
map.n('<LocalLeader>]', "<Plug>TableSwapRight", "Swap table column right", {buffer=true})
map.n('<LocalLeader>{', "<Plug>TableAddColLeft", "Add new empty column to the left", {buffer=true})
map.n('<LocalLeader>}', "<Plug>TableAddColRight", "Add new empty column to the right", {buffer=true})
map.n('<LocalLeader><left>', "<Plug>TableGoLeft", "Goto left cell", {buffer=true})
map.n('<LocalLeader><right>', "<Plug>TableGoRight", "Goto right cell", {buffer=true})
map.n('<LocalLeader><up>', "<Plug>TableGoUp", "Goto up cell", {buffer=true})
map.n('<LocalLeader><down>', "<Plug>TableGoDown", "Goto down cell", {buffer=true})
-- can't use backspace since it is hardcoded for up one level in mini.clue
map.n('<LocalLeader><del>', "<Plug>(vimtex-clean)", "Clean (rm aux)", {buffer=true})
map.n('<LocalLeader><S-del>', "<Plug>(vimtex-clean-all)", "Clean all (rm aux+out)", {buffer=true})
map.desc('n', '<LocalLeader>e', "Errors")
-- hacky. VimtexErrors puts errors found by Vimtex in quickfix (should be
-- running, use <leader>Lb) then cclose closes quickfix, and then Telescope
-- opens the quickfix in a nicer view.
map.n('<space>E', "<Cmd>VimtexErrors<CR>|:cclose|<Cmd>Telescope quickfix<CR>", "Errors", {buffer=true})
-- Avoid accidentally deleting aux files with default keymap that is very similar to the compile keymaps.
-- no-operation, instead of del since deleting throws error and this means we don't start a vim change motion etc.
map.n("<LocalLeader>c", "<nop>", nil, { buffer = true })
map.n("<LocalLeader>C", "<nop>", nil, { buffer = true })

map.n('<Leader>cg', function()
    local auxs = vim.fs.find("aux", { upward = true, limit = 5 })
    if #auxs == 0 then
        vim.api.nvim_err_writeln("Couldn't makeglossaries. No aux/ dir found (searched upward)")
        return
    end
    local makeglossaries = function(on_exit)
        vim.system({ "makeglossaries", "main" }, { text = true, cwd = auxs[1] }, on_exit)
    end
    makeglossaries(function(obj)
        if obj.code == 0 then
            print("makeglossaries complete")
            return
        end
        -- Retry once. Fixes a specific error.
        makeglossaries(function(obj)
            if obj.code == 0 then
                print("makeglossaries complete")
                return
            end
            local stderr = obj.stderr:gsub("\n$", "")
            vim.schedule(function() -- notify when we are ready
                vim.notify(stderr)  -- vim.notify instead of print to see multiple lines
            end)
        end)
    end)
end, "Compile glossary", {buffer=false}) -- No need to make buffer local

---Run `biber --cache` to get biber cache dir, so far found in /var/folders/...
---Then call `rm -rf` on it.
---Then call the supplied on_exit function if given.
---@param on_exit function
local function biber_clear_cache(on_exit)
    vim.system({ "biber", "--cache" }, { text = true }, function(obj)
        if obj.code ~= 0 then
            util.schedule_notify(obj)
        else
            local cachepath = obj.stdout:gsub('\n$', '')
            if cachepath:match("^/var/folders/") then
                vim.system({ "rm", "-rf", cachepath }, { text = true }, function(obj)
                    if obj.code ~= 0 then
                        util.schedule_notify(obj)
                    else
                        print("biber cache cleared")
                        if on_exit ~= nil then
                            on_exit()
                        end
                    end
                end)
            else
                print("Unexpected cache path: " .. cachepath)
            end
        end
    end)
end
map.cmd("BiberClearCache", biber_clear_cache)

map.n('<leader>cb', function(main)
    if main == nil then main = "main" end
    local biber = function(on_exit)
        vim.system({ "biber", main }, { text = true }, on_exit)
    end
    biber(function(obj)
        if obj.code == 0 then
            print("biber complete")
            return
        end
        -- Might have failed due to lack of pdflatex/lualatex/etc. run
        if obj.stdout:contains("ERROR - Cannot find '" .. main .. ".bcf'!") then
            print("biber failed: no main.bcf")
            return
        end
        -- Otherwise clear cache then retry once.
        -- Clearing cache helped with cryptic error with message:
        -- Unicode::UCD: failed to find unicore/version in /var/folders/...
        biber_clear_cache(function()
            biber(function(obj)
                if obj.code == 0 then
                    print("biber complete")
                    return
                end
                util.schedule_notify(obj)
            end)
        end)
    end)
end, "Compile bibliography", {buffer=false}) -- No need for buffer local

map.n('gK', function()
    local line = vim.api.nvim_get_current_line()
    local pac = line:match("\\usepackage.*{([%w_-]+)}")
    if pac == nil then
        return print("No package name found on line.")
    end
    local obj = vim.system({ 'texdoc', pac }):wait()
    if obj.code ~= 0 then
        print("texdoc " .. pac .. " failed.")
    end
end, "texdoc help", {buffer=true})
-- gx can open ctan main site, see plugins/init.lua
map.n('gX', function()
    local line = vim.api.nvim_get_current_line()
    local pac = line:match("\\usepackage.*{([%w_-]+)}")
    if pac ~= nil then
        local ctan = "https://ctan.org/pkg/" .. pac .. "?lang=en"
        -- open manual pdf directly.
        -- We read it from ctan site since the url may vary (e.g. font doc at https://au.mirrors.cicku.me/ctan/fonts/baskervillef/doc/baskervillef-doc.pdf)
        print("Opening manual pdf(s)...")
        return vim.system({ 'curl', ctan }, {
            text = true,
            stdout = function(err, data)
                if err ~= nil then
                    print("Error opening ctan manual pdf.")
                else
                    vim.system({ 'grep', '-o', [["http[^"]*\.pdf"]] }, { text = true, stdin = data }, function(obj)
                        if obj.code ~= 0 then
                            util.schedule_notify(obj)
                        else
                            for _, quoted_url in ipairs(vim.split(obj.stdout, '\n')) do
                                util.open(quoted_url:match('^"(.*)"$'))
                            end
                            return
                        end
                    end)
                end
            end
        })
    end
end, "Open CTAN manual(s) for package", {buffer=true})

map.n('<LocalLeader>s', "<Plug>(vimtex-status)", "Status", {buffer=true})
map.n('<LocalLeader>S', "<Plug>(vimtex-status-all)", "Status all", {buffer=true})
map.desc('n', '<LocalLeader>i', "Info")
map.desc('n', '<LocalLeader>I', "Info full")
map.desc('n', '<LocalLeader>q', "Log of VimTeX actions")
-- j for jump. Not using p for preamble since I want to use it for pasting tables.
map.n("<LocalLeader>j", "<Plug>TableJumpPre", "goto/from preamble (table)", {buffer=true})
map.n('<LocalLeader>m', "<Plug>(vimtex-toggle-main)", "Toggle compiling main vs subfile", {buffer=true})
map.desc('n', '<LocalLeader>t', "TOC open")
map.desc('n', '<LocalLeader>T', "TOC toggle")
map.n("<LocalLeader>u", "<Plug>Latex2Unicode", "TeX -> unicode", {buffer=true})
map.nx("<LocalLeader>U", "<Plug>Unicode2Latex", "Unicode -> TeX", {buffer=true})
map.n('<leader>cv', '<Plug>(vimtex-view)', "View", {buffer=true})
map.n('<LocalLeader>r', "<Plug>(vimtex-reload)", "Reload", {buffer=true})
map.n('<LocalLeader>R', "<Plug>(vimtex-reload-state)", "Reload state", {buffer=true})
map.n('<leader>ck', '<Plug>(vimtex-stop)', "Stop", {buffer=true})
map.n('<leader>cK', '<Plug>(vimtex-stop-all)', "Stop all", {buffer=true})
map.n('<leader>cc', function()
    -- check for already running latexmk, since running it twice will break things.
    latexmk.is_running(function(is_running)
        if is_running then
            print("Looks like latexmk is already running, e.g. in another tab. Use :VimtexComileSS to force.")
        else
            vim.schedule(vim.fn["vimtex#compiler#compile_ss"])
        end
    end)
end, "Compile single shot", {buffer=true})
map.n('<leader>cC', function()
    latexmk.is_running(function(is_running)
        if is_running then
            print("Looks like latexmk is already running, e.g. in another tab. Use :VimtexComile to force.")
        else
            vim.schedule(vim.fn["vimtex#compiler#compile"])
        end
    end)
end, "Compile continuously", {buffer=true})
map.n('<leader>cl', '<Plug>(vimtex-compile-output)', "Output", {buffer=true})
map.n('<leader>cc', '<Plug>(vimtex-compile-selected)', "Compile selected", {buffer=true})
map.n('<CR>', "<plug>(vimtex-compile-selected)", "Compile motion", {buffer=true})
map.n('<CR>', "<plug>(vimtex-compile-selected)", "Compile selection", {buffer=true})

-- e.g. \section*{}
map.desc('n', 'tsc', "Cmd/Star")
map.desc('n', 'tse', "Env/Star")
-- e.g. with(out) \left
map.desc('n', 'tsd', "Delim")
map.n('<LocalLeader>)', "<plug>(vimtex-delim-toggle-modifier)", "Delim", {buffer=true})
map.n('<LocalLeader>(', "<plug>(vimtex-delim-toggle-modifier-reverse)", "Delim", {buffer=true})
-- same as d, but looks through g:vimtex_delim_toggle_mod_list in reverse
map.desc('n', 'tsD', "Delim rev")
-- toggle / <-> \frac
map.desc('n', 'tsf', "Fraction")
-- change surrounding ...
map.desc('n', 'csc', "Cmd")
map.desc('n', 'cse', "Env")
map.desc('n', 'csm', "Math")
-- in/around ...
map.desc({ 'o', 'x' }, 'id', "Delim")
map.desc({ 'o', 'x' }, 'iP', "Section")
map.n("ts$", "<Plug>(vimtex-env-toggle-math)", "Inline <-> display", {buffer=true})
map.n("ts4", "<Plug>(vimtex-env-toggle-math)", "Inline <-> display", {buffer=true})
map.n("tsm", "<plug>(vimtex-env-toggle-math)", "Inline <-> display", {buffer=true})
map.n("dsm", "<plug>(vimtex-env-delete-math)", "Delete math", {buffer=true})
map.n("csm", "<plug>(vimtex-env-change-math)", "Change math", {buffer=true})
map.n("xad", "yaddad", "Cut a delim", { remap = true, buffer=true })
map.n("xid", "yiddid", "Cut in delim", { remap = true, buffer=true })
-- item with i instead of m and math with m
map.o("ai", "<Plug>(vimtex-am)", "An item", {buffer=true})
map.o("ii", "<Plug>(vimtex-im)", "In item", {buffer=true})
map.o("am", "<Plug>(vimtex-a$)", "An eq", {buffer=true})
map.o("im", "<Plug>(vimtex-i$)", "In eq", {buffer=true})
-- shorthand to $ just using 4 ($ without shift)
map.o("a4", "<Plug>(vimtex-a$)", "An eq", {buffer=true})
map.o("i4", "<Plug>(vimtex-i$)", "In eq", {buffer=true})
-- next/prev start/end of ...
map.n("[m", "<Plug>(vimtex-[n)", "Math start", {buffer=true})
map.n("[M", "<Plug>(vimtex-[N)", "Math end", {buffer=true})
map.n("[4", "<Plug>(vimtex-[n)", "Math start", {buffer=true})
map.n("[$", "<Plug>(vimtex-[N)", "Math end", {buffer=true})
map.n("]m", "<Plug>(vimtex-]n)", "Math start", {buffer=true})
map.n("]M", "<Plug>(vimtex-]N)", "Math end", {buffer=true})
map.n("]4", "<Plug>(vimtex-]n)", "Math start", {buffer=true})
map.n("]$", "<Plug>(vimtex-]N)", "Math end", {buffer=true})

---Get latex root.
---@param main string? name of main tex file, by default "main".
---@return string? path to root directory.
local function get_root(main)
    if main == nil then
        main = "main.tex"
    else
        -- add .tex extension if not given
        main = main:gsub("%.%w+", "") .. ".tex"
    end
    local mains = vim.fs.find(main, { upward = true, limit = 5 })
    if #mains == 0 then
        return nil
    else
        return vim.fs.dirname(mains[1])
    end
end

---Get mappings of the project input tree.
---@param root string Root folder to search. E.g. use get_root()
---@return table input_mapping Keys are referenced files without extension and values are {parent,linenum},
---where parent is the file that includes them relative to root also without extension, and linenum is the linenumber within parent.
---Since keys are unique this assumes a file is only referenced in one place.
local function get_inputs(root)
    local stdouts = ""
    for _, cmd in ipairs {
        -- rg '^[^%]*\\input\{([\w./_-]+)\}' -t tex -r '$1' --vimgrep --no-column,
        -- rg '^[^%]*\\include\{([\w./_-]+)\}' -t tex -r '$1' --vimgrep --no-column,
        -- rg '^[^%]*\\import\{([\w./_-]+)\}\{([\w./_-]+)\}' -t tex -r '$1$2' --vimgrep --no-column,
        [[rg ^[^%]*\\input\{([\w./_-]+)\} -t tex -r $1 --vimgrep --no-column]],
        [[rg ^[^%]*\\include\{([\w./_-]+)\} -t tex -r $1 --vimgrep --no-column]],
        [[rg ^[^%]*\\import\{([\w./_-]+)\}\{([\w./_-]+)\} -t tex -r $1$2 --vimgrep --no-column]],
    } do
        stdouts = stdouts .. vim.system(vim.split(cmd, ' '), { text = true, cwd = root }):wait().stdout
    end
    local input_mapping = {}
    for filepath, linenum, input in stdouts:gmatch("([^:]+)%.tex:(%d+):([^:]+)\n") do
        input_mapping[input:gsub("%.tex$", "")] = { filepath, linenum }
    end
    return input_mapping
end

---Get the filename of a file that inputs, includes or imports a given file.
---@param filepath string?
---@return string? parent Absolute path of parent file with extension.
---@return integer? linenum Line number within parent where the given file is inputted.
local function get_file_reference(filepath)
    filepath = filepath or vim.api.nvim_buf_get_name(0)
    local root = get_root()
    if root == nil then return nil end
    -- filepath relative to root without extension
    filepath = filepath:gsub("^" .. root .. "/?", ""):gsub("%.tex", "")
    local input_mapping = get_inputs(root)
    local input = input_mapping[filepath]
    if input == nil then
        return nil
    else
        return vim.fs.joinpath(root, input[1]) .. ".tex", input[2]
    end
end

map.n('<LocalLeader>-', function()
    local parent, linenum = get_file_reference()
    if parent == nil then
        print("No file references found.")
    else
        vim.cmd("edit +" .. linenum .. " " .. parent)
    end
end, "Go up in latex structure", {buffer=true})

-- colorscheme aucmd to fix missing or inconsistent hl links
vim.api.nvim_create_autocmd("ColorScheme", {
    buffer = 0,
    group = vim.api.nvim_create_augroup("Tex", { clear = true }),
    callback = function()
        -- Pick a reduced colour for removing emphasis on things like \cite{...} where the body's color and underline gives it emphasis by itself.
        -- We want to differentiate from comment and nontext, and nontext is bold so the fg with italic should be enough differentiation, plus we would write comments more that using nontext.
        local gray = hi.fg("NonText")
        hi.link("texCmd", "@function.call")
        hi.link("texCmdEnv", "@keyword.function") -- italic instead of bold for begin end
        hi.link("texCmdRef", "@function.builtin") -- italic
        -- italic \section{...}, bold etc. Gray a bit since the "..." shows aesthetic
        hi.set("texCmdPart", { fg = gray, italic = true })
        hi.set("texCmdStyleBold", { fg = gray, italic = true })
        hi.set("texCmdStyleItal", { fg = gray, italic = true })
        hi.set("texTypeStyle", { fg = gray, italic = true })       -- e.g. \underline
        hi.set("texItalStyle", { italic = true })                  -- contents of \emph{...}
        hi.set("texCmdRefConcealed", { fg = gray, italic = true }) -- italic \cite
        hi.set("texCmdRef", { fg = gray, italic = true })
        hi.set("texCmdCRef", { fg = gray, italic = true })
        hi.set("texCmdAcro", { fg = gray })           -- custom cmd defined in after/syntax/tex.vim
        hi.link("texCmdPackage", "@function.builtin") -- italic \package
        hi.link("texCmdInput", "@function.builtin")   -- italic \inputgraphics
        hi.link("texCmdTitle", "@function.builtin")   -- italic \title
        hi.link("texCmdAuthor", "@function.builtin")  -- italic \author
        hi.link("texCmdLet", "@function.builtin")     -- italic \let
        hi.link("texStatement", "@function.builtin")  -- only seen for \mathrm so far
        hi.mod("texMatcher", { underline = true })    -- matched parenthesis, \underline body, etc.
        hi.link("texEnvArgName", "@method")           -- bold and shine instead of nothing
        hi.link("texCmdBeamer", "@function")
        hi.link("texOpt", "@parameter")
        hi.link("texBeamerOpt", "@parameter")
        hi.link("texOptEqual", "@operator")
        hi.link("texArg", "@parameter")
        hi.link("texFileArg", "@string")
        hi.link("texFilesArg", "@string")
        hi.link("texFileOpt", "@parameter")
        hi.link("TexBeamerDelim", "Delimiter")
        hi.link("superscript", "Type")                                               -- like \huge, \normalsize etc
        hi.link("subscript", "Type")                                                 -- like \huge, \normalsize etc
        hi.set("texRefConcealedArg", { fg = hi.fg("TexFileArg"), underline = true }) -- body of \cite{...}
        hi.link("texTitleArg", "Title")
        hi.link("texPartArgTitle", "Title")
        hi.link("texRefArg", "@tag")                                    -- body of \label
        hi.link("texSpecialChar", "@comment")                           -- unbreakable space ~, and \&
        hi.link("texMathZone", "@number")                               -- Most of tex math zone that isn't captured by anything else (such as math functions) is numbers and we don't use numbers much elsewhere.
        hi.set("texMathCmdText", { fg = gray, italic = true })          -- italic \text in math mode
        hi.set("texMathSymbol", { fg = hi.fg("@type"), italic = true }) -- type is similar colour to number
        hi.set("texMathSymbol", { fg = hi.fg("@type"), italic = true }) --
        hi.link("texSICmd", "@number")                                  -- not bold SI. Color like math mode
        hi.set("texLigature", { bold = true })                          -- bold instead of strong color to only give subtle focus to ``'', --, and the ' in don't
        hi.link("texCmdLigature", "@function.call")
        hi.mod("texCmdLigature", { italic = true })
        hi.link("texTabularChar", "Operator") -- & and \\ in tables. Could also use Delimiter but this makes them bold.
        hi.mod("texCmdClass", { italic = true, bold = true })
        hi.link("texOptSep", "Delimiter")
        hi.mod("texCmdDef", { bold = true, italic = true })    -- an actual function definition. \def. TeX primitive.
        hi.mod("texCmdNewcmd", { bold = true, italic = true }) -- an actual function definition. \newcommand. LaTeX wrapper on def.
        hi.link("texNewcmdArgName", "@parameter")
    end
})
