#!/usr/bin/env lua
require "tex.overleaf"
local tbl = require "tex.tables"
require "tex.cmds"
require "tex.textcolor"
local hi = require "utils/highlights"
local util = require "utils/init" -- string.contains and schedule_notify
local latexmk = require "tex.latexmk"

-- also defined in lua/plugins/mini.lua as shift+` but we can just use ` since 
-- I don't think it has use in latex, maybe except in some verbatim code block or something?
require'mini.surround'.config.custom_surroundings['`'] = {
    input = { "``().-()''" },
    output = { left = '``', right = "''" },
}

-- requires kana/vim-textobj-user, see lua/plugins/init.lua
vim.fn["textobj#user#plugin"]("tex", {
    ['latex-ticks'] = {
        pattern = {'``', "''"},
        ['select-a'] = '<buffer> a`',
        ['select-i'] = '<buffer> i`',
    },
})

-- keymaps assume vimtex is used.

local mapbuf = function (keys, rhs, desc, mode, opts)
    opts = opts or {}
    opts.desc = desc
    vim.keymap.set(mode or 'n', keys, rhs, opts)
end
local set_keymap_desc = function(...) pcall(require"mini.clue".set_keymap_desc, ...) end

set_keymap_desc('n', '<LocalLeader>a', "Context menu")
mapbuf('<LocalLeader>A', "<plug>TableAlign", "Align table")
mapbuf("<LocalLeader>y", "<plug>TableYank", "Yank as TSV")
mapbuf("<LocalLeader>p", "<plug>TablePaste", "Paste TSV", {'n', 'v'})
mapbuf('da|', "<Plug>TableDelCol", "Delete a table column")
mapbuf("dix", tbl.deleteInCell, "Delete in cell", 'n', {silent=true})
mapbuf("cix", tbl.changeInCell, "Change in cell", 'n', {silent=true})
mapbuf("xix", function () tbl.deleteInCell("+") end, "Cut in cell", 'n', {silent=true})
mapbuf("ix", "<Plug>TableSelInCell", "Select in cell", 'n', {silent=true,})
mapbuf('<LocalLeader>[', "<Plug>TableSwapLeft", "Swap table column left")
mapbuf('<LocalLeader>]', "<Plug>TableSwapRight", "Swap table column right")
mapbuf('<LocalLeader>{', "<Plug>TableAddColLeft", "Add new empty column to the left")
mapbuf('<LocalLeader>}', "<Plug>TableAddColRight", "Add new empty column to the right")
mapbuf('<leader>fc', "<Cmd>Telescope bibtex<CR>", "Cite")
mapbuf('<LocalLeader><left>', "<Plug>TableGoLeft", "Goto left cell")
mapbuf('<LocalLeader><right>', "<Plug>TableGoRight", "Goto right cell")
mapbuf('<LocalLeader><up>', "<Plug>TableGoUp", "Goto up cell")
mapbuf('<LocalLeader><down>', "<Plug>TableGoDown", "Goto down cell")
mapbuf('<LocalLeader><del>', "<Plug>(vimtex-clean)", "Clean (rm aux)")
mapbuf('<LocalLeader><S-del>', "<Plug>(vimtex-clean-all)", "Clean all (rm aux+out)")
set_keymap_desc('n', '<LocalLeader>e', "Errors")
-- hacky. VimtexErrors puts errors found by Vimtex in quickfix (should be
-- running, use <leader>Lb) then cclose closes quickfix, and then Telescope
-- opens the quickfix in a nicer view.
mapbuf('<space>E', "<Cmd>VimtexErrors<CR>|:cclose|<Cmd>Telescope quickfix<CR>", "Errors")
-- Avoid accidentally deleting aux files with default keymap that is very similar to the compile keymaps.
-- no-operation, instead of del since deleting throws error and this means we don't start a vim change motion etc.
vim.keymap.set('n', "<LocalLeader>c", "<nop>", {buffer=true})
vim.keymap.set('n', "<LocalLeader>C", "<nop>", {buffer=true})

mapbuf('<Leader>cg', function ()
    local auxs = vim.fs.find("aux", {upward=true, limit=5})
    if #auxs == 0 then
        vim.api.nvim_err_writeln("Couldn't makeglossaries. No aux/ dir found (searched upward)")
        return
    end
    local makeglossaries = function(on_exit)
        vim.system({"makeglossaries", "main"}, {text=true, cwd=auxs[1]}, on_exit)
    end
    makeglossaries(function (obj)
        if obj.code == 0 then
            print("makeglossaries complete")
            return
        end
        -- Retry once. Fixes a specific error.
        makeglossaries(function (obj)
            if obj.code == 0 then
                print("makeglossaries complete")
                return
            end
            local stderr = obj.stderr:gsub("\n$", "")
            vim.schedule(function () -- notify when we are ready
                vim.notify(stderr) -- vim.notify instead of print to see multiple lines
            end)
        end)
    end)
end, "Compile glossary")

---Run `biber --cache` to get biber cache dir, so far found in /var/folders/...
---Then call `rm -rf` on it.
---Then call the supplied on_exit function if given.
---@param on_exit function
local function biber_clear_cache(on_exit)
    vim.system({"biber", "--cache"}, {text=true}, function (obj)
        if obj.code ~= 0 then
            util.schedule_notify(obj)
        else
            local cachepath = obj.stdout:gsub('\n$', '')
            if cachepath:match("^/var/folders/") then
                vim.system({"rm", "-rf", cachepath}, {text=true}, function (obj)
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
vim.api.nvim_create_user_command("BiberClearCache", biber_clear_cache, {})

mapbuf('<leader>cb', function (main)
    if main == nil then main = "main" end
    local biber = function(on_exit)
        vim.system({"biber", main}, {text=true}, on_exit)
    end
    biber(function (obj)
        if obj.code == 0 then
            print("biber complete")
            return
        end
        -- Might have failed due to lack of pdflatex/lualatex/etc. run
        if obj.stdout:contains("ERROR - Cannot find '"..main..".bcf'!") then
            print("biber failed: no main.bcf")
            return
        end
        -- Otherwise clear cache then retry once.
        -- Clearing cache helped with cryptic error with message:
        -- Unicode::UCD: failed to find unicore/version in /var/folders/...
        biber_clear_cache(function ()
            biber(function (obj)
                if obj.code == 0 then
                    print("biber complete")
                    return
                end
                util.schedule_notify(obj)
            end)
        end)
    end)
end, "Compile bibliography")

mapbuf('gK', function ()
    local line = vim.api.nvim_get_current_line()
    local pac = line:match("\\usepackage.*{([%w_-]+)}")
    if pac == nil then
        return print("No package name found on line.")
    end
    local obj = vim.system({'texdoc', pac}):wait()
    if obj.code ~= 0 then
        print("texdoc " .. pac .. " failed.")
    end
end, "texdoc help")
-- gx can open ctan main site, see plugins/init.lua
mapbuf('gX', function ()
    local line = vim.api.nvim_get_current_line()
    local pac = line:match("\\usepackage.*{([%w_-]+)}")
    if pac ~= nil then
        local ctan = "https://ctan.org/pkg/"..pac.."?lang=en"
        -- open manual pdf directly.
        -- We read it from ctan site since the url may vary (e.g. font doc at https://au.mirrors.cicku.me/ctan/fonts/baskervillef/doc/baskervillef-doc.pdf)
        print("Opening manual pdf(s)...")
        return vim.system({'curl', ctan}, {text=true, stdout=function (err, data)
            if err ~= nil then
                print("Error opening ctan manual pdf.")
            else
                vim.system({'grep', '-o', [["http[^"]*\.pdf"]]}, {text=true, stdin=data}, function (obj)
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
        end})
    end
end, "Open CTAN manual(s) for package")

mapbuf('<LocalLeader>s', "<Plug>(vimtex-status)", "Status")
mapbuf('<LocalLeader>S', "<Plug>(vimtex-status-all)", "Status all")
set_keymap_desc('n', '<LocalLeader>i', "Info")
set_keymap_desc('n', '<LocalLeader>I', "Info full")
set_keymap_desc('n', '<LocalLeader>q', "Log of VimTeX actions")
-- j for jump. Not using p for preamble since I want to use it for pasting tables.
mapbuf("<LocalLeader>j", "<Plug>TableJumpPre", "goto/from preamble (table)")
mapbuf('<LocalLeader>m', "<Plug>(vimtex-toggle-main)", "Toggle compiling main vs subfile")
set_keymap_desc('n', '<LocalLeader>t', "TOC open")
set_keymap_desc('n', '<LocalLeader>T', "TOC toggle")
mapbuf("<LocalLeader>u", "<Plug>Latex2Unicode", "TeX -> unicode")
mapbuf("<LocalLeader>U", "<Plug>Unicode2Latex", "Unicode -> TeX", {'n', 'x'})
mapbuf('<leader>cv', '<Plug>(vimtex-view)', "View")
mapbuf('<LocalLeader>r', "<Plug>(vimtex-reload)", "Reload")
mapbuf('<LocalLeader>R', "<Plug>(vimtex-reload-state)", "Reload state")
mapbuf('<leader>ck', '<Plug>(vimtex-stop)', "Stop")
mapbuf('<leader>cK', '<Plug>(vimtex-stop-all)', "Stop all")
mapbuf('<leader>cc', function ()
    -- check for already running latexmk, since running it twice will break things.
    latexmk.is_running(function (is_running)
        if is_running then
            print("Looks like latexmk is already running, e.g. in another tab. Use :VimtexComileSS to force.")
        else
            vim.schedule(vim.fn["vimtex#compiler#compile_ss"])
        end
    end)
end, "Compile single shot")
mapbuf('<leader>cC', function ()
    latexmk.is_running(function (is_running)
        if is_running then
            print("Looks like latexmk is already running, e.g. in another tab. Use :VimtexComile to force.")
        else
            vim.schedule(vim.fn["vimtex#compiler#compile"])
        end
    end)
end, "Compile continuously")
mapbuf('<leader>cl', '<Plug>(vimtex-compile-output)', "Output")
mapbuf('<leader>cc', '<Plug>(vimtex-compile-selected)', "Compile selected")
mapbuf('<CR>', "<plug>(vimtex-compile-selected)", "Compile motion")
mapbuf('<CR>', "<plug>(vimtex-compile-selected)", "Compile selection")

-- e.g. \section*{}
set_keymap_desc('n', 'tsc', "Cmd/Star")
set_keymap_desc('n', 'tse', "Env/Star")
-- e.g. with(out) \left 
set_keymap_desc('n', 'tsd', "Delim")
mapbuf('<LocalLeader>)', "<plug>(vimtex-delim-toggle-modifier)", "Delim")
mapbuf('<LocalLeader>(', "<plug>(vimtex-delim-toggle-modifier-reverse)", "Delim")
-- same as d, but looks through g:vimtex_delim_toggle_mod_list in reverse
set_keymap_desc('n', 'tsD', "Delim rev")
-- toggle / <-> \frac
set_keymap_desc('n', 'tsf', "Fraction")
-- change surrounding ...
set_keymap_desc('n', 'csc', "Cmd")
set_keymap_desc('n', 'cse', "Env")
set_keymap_desc('n', 'csm', "Math")
-- in/around ...
set_keymap_desc({'o', 'x'}, 'id', "Delim")
set_keymap_desc({'o', 'x'}, 'iP', "Section")
mapbuf("ts$", "<Plug>(vimtex-env-toggle-math)", "Inline <-> display")
mapbuf("ts4", "<Plug>(vimtex-env-toggle-math)", "Inline <-> display")
mapbuf("tsm", "<plug>(vimtex-env-toggle-math)", "Inline <-> display")
mapbuf("dsm", "<plug>(vimtex-env-delete-math)", "Delete math")
mapbuf("csm", "<plug>(vimtex-env-change-math)", "Change math")
mapbuf("xad", "yaddad", "Cut a delim", 'n', {remap=true,})
mapbuf("xid", "yiddid", "Cut in delim", 'n', {remap=true,})
-- item with i instead of m and math with m
mapbuf("ai", "<Plug>(vimtex-am)", "An item")
mapbuf("ii", "<Plug>(vimtex-im)", "In item")
mapbuf("am", "<Plug>(vimtex-a$)", "An eq")
mapbuf("im", "<Plug>(vimtex-i$)", "In eq")
-- shorthand to $ just using 4 ($ without shift)
mapbuf("a4", "<Plug>(vimtex-a$)", "An eq")
mapbuf("i4", "<Plug>(vimtex-i$)", "In eq")
-- next/prev start/end of ...
mapbuf("[m", "<Plug>(vimtex-[n)", "Math start")
mapbuf("[M", "<Plug>(vimtex-[N)", "Math end")
mapbuf("[4", "<Plug>(vimtex-[n)", "Math start")
mapbuf("[$", "<Plug>(vimtex-[N)", "Math end")
mapbuf("]m", "<Plug>(vimtex-]n)", "Math start")
mapbuf("]M", "<Plug>(vimtex-]N)", "Math end")
mapbuf("]4", "<Plug>(vimtex-]n)", "Math start")
mapbuf("]$", "<Plug>(vimtex-]N)", "Math end")

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
    local mains = vim.fs.find(main, {upward=true, limit=5})
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
    for _, cmd in ipairs{
        -- rg '^[^%]*\\input\{([\w./_-]+)\}' -t tex -r '$1' --vimgrep --no-column,
        -- rg '^[^%]*\\include\{([\w./_-]+)\}' -t tex -r '$1' --vimgrep --no-column,
        -- rg '^[^%]*\\import\{([\w./_-]+)\}\{([\w./_-]+)\}' -t tex -r '$1$2' --vimgrep --no-column,
        [[rg ^[^%]*\\input\{([\w./_-]+)\} -t tex -r $1 --vimgrep --no-column]],
        [[rg ^[^%]*\\include\{([\w./_-]+)\} -t tex -r $1 --vimgrep --no-column]],
        [[rg ^[^%]*\\import\{([\w./_-]+)\}\{([\w./_-]+)\} -t tex -r $1$2 --vimgrep --no-column]],
    } do
        stdouts = stdouts .. vim.system(vim.split(cmd, ' '), {text=true, cwd=root}):wait().stdout
    end
    local input_mapping = {}
    for filepath, linenum, input in stdouts:gmatch("([^:]+)%.tex:(%d+):([^:]+)\n") do
        input_mapping[input:gsub("%.tex$", "")] = {filepath, linenum}
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
    filepath = filepath:gsub("^"..root.."/?", ""):gsub("%.tex", "")
    local input_mapping = get_inputs(root)
    local input = input_mapping[filepath]
    if input == nil then
        return nil
    else
        return vim.fs.joinpath(root, input[1]) .. ".tex", input[2]
    end
end

mapbuf('<LocalLeader>-', function ()
    local parent, linenum = get_file_reference()
    if parent == nil then
        print("No file references found.")
    else
        vim.cmd("edit +" .. linenum .. " " .. parent)
    end
end, "Go up in latex structure")

-- colorscheme aucmd to fix missing or inconsistent hl links
local grp = vim.api.nvim_create_augroup("Tex", {clear=true})
vim.api.nvim_create_autocmd("Colorscheme", {
    buffer = 0,
    group = grp,
    callback = function ()
        -- Pick a reduced colour for removing emphasis on things like \cite{...} where the body's color and underline gives it emphasis by itself.
        -- We want to differentiate from comment and nontext, and nontext is bold so the fg with italic should be enough differentiation, plus we would write comments more that using nontext.
        local gray = hi.getfg("NonText")
        hi.link("texCmd", "@function.call")
        hi.link("texCmdEnv", "@keyword.function") -- italic instead of bold for begin end
        hi.link("texCmdRef", "@function.builtin") -- italic
        -- italic \section{...}, bold etc. Gray a bit since the "..." shows aesthetic
        hi.set("texCmdPart", {fg=gray, italic=true})
        hi.set("texCmdStyleBold", {fg=gray, italic=true})
        hi.set("texCmdStyleItal", {fg=gray, italic=true})
        hi.set("texTypeStyle", {fg=gray, italic=true}) -- e.g. \underline
        hi.set("texItalStyle", {italic=true}) -- contents of \emph{...}
        hi.set("texCmdRefConcealed", {fg=gray, italic=true}) -- italic \cite
        hi.set("texCmdRef", {fg=gray, italic=true})
        hi.set("texCmdCRef", {fg=gray, italic=true})
        hi.set("texCmdAcro", {fg=gray}) -- custom cmd defined in after/syntax/tex.vim
        hi.link("texCmdPackage", "@function.builtin") -- italic \package
        hi.link("texCmdInput", "@function.builtin") -- italic \inputgraphics
        hi.link("texCmdTitle", "@function.builtin") -- italic \title
        hi.link("texCmdAuthor", "@function.builtin") -- italic \author
        hi.link("texCmdLet", "@function.builtin") -- italic \let
        hi.link("texStatement", "@function.builtin") -- only seen for \mathrm so far
        hi.mod("texMatcher", {underline=true}) -- matched parenthesis, \underline body, etc.
        hi.link("texEnvArgName", "@method") -- bold and shine instead of nothing
        hi.link("texCmdBeamer", "@function")
        hi.link("texOpt", "@parameter")
        hi.link("texBeamerOpt", "@parameter")
        hi.link("texOptEqual", "@operator")
        hi.link("texArg", "@parameter")
        hi.link("texFileArg", "@string")
        hi.link("texFilesArg", "@string")
        hi.link("texFileOpt", "@parameter")
        hi.link("TexBeamerDelim", "Delimiter")
        hi.link("superscript", "Type") -- like \huge, \normalsize etc
        hi.link("subscript", "Type") -- like \huge, \normalsize etc
        hi.set("texRefConcealedArg", {fg=hi.getfg("TexFileArg"), underline=true}) -- body of \cite{...}
        hi.link("texTitleArg", "Title")
        hi.link("texPartArgTitle", "Title")
        hi.link("texRefArg", "@tag") -- body of \label
        hi.link("texSpecialChar", "@comment") -- unbreakable space ~, and \&
        hi.link("texMathZone", "@number") -- Most of tex math zone that isn't captured by anything else (such as math functions) is numbers and we don't use numbers much elsewhere.
        hi.set("texMathCmdText", {fg=gray, italic=true}) -- italic \text in math mode
        hi.set("texMathSymbol", {fg=hi.getfg("@type"), italic=true}) -- type is similar colour to number
        hi.set("texMathSymbol", {fg=hi.getfg("@type"), italic=true}) -- 
        hi.link("texSICmd", "@number") -- not bold SI. Color like math mode
        hi.set("texLigature", {bold=true}) -- bold instead of strong color to only give subtle focus to ``'', --, and the ' in don't
        hi.link("texCmdLigature", "@function.call")
        hi.mod("texCmdLigature", {italic=true})
        hi.link("texTabularChar", "Operator") -- & and \\ in tables. Could also use Delimiter but this makes them bold.
        hi.mod("texCmdClass", {italic=true, bold=true})
        hi.link("texOptSep", "Delimiter")
        hi.mod("texCmdDef", {bold=true, italic=true}) -- an actual function definition. \def. TeX primitive.
        hi.mod("texCmdNewcmd", {bold=true, italic=true}) -- an actual function definition. \newcommand. LaTeX wrapper on def.
        hi.link("texNewcmdArgName", "@parameter")
    end
})

