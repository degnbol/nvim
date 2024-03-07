#!/usr/bin/env lua
require "tex.overleaf"
local tbl = require "tex.tables"
require "tex.cmds"
require "tex.textcolor"

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
require "utils/keymap"
set_keymap_desc('n', '<leader><leader>a', "Context menu")
vim.keymap.set('n', '<leader><leader>A', "<plug>TableAlign", { buffer=true, desc="Align table" })
vim.keymap.set('n', 'da|', "<Plug>TableDelCol", { buffer=true, desc="Delete a table column" })
vim.keymap.set("n", "dix", tbl.deleteInCell, { buffer=true, silent=true, desc="Delete in cell", })
vim.keymap.set("n", "cix", tbl.changeInCell, { buffer=true, silent=true, desc="Change in cell", })
vim.keymap.set("n", "xix", function () tbl.deleteInCell("+") end, { buffer=true, silent=true, desc="Cut in cell", })
vim.keymap.set("v", "ix", "<Plug>TableSelInCell", { buffer=true, silent=true, desc="Select in cell", })
vim.keymap.set('n', '<leader><leader>[', "<Plug>TableSwapLeft", { buffer=true, desc="Swap table column left"})
vim.keymap.set('n', '<leader><leader>]', "<Plug>TableSwapRight", { buffer=true, desc="Swap table column right"})
vim.keymap.set('n', '<leader><leader>{', "<Plug>TableAddColLeft", { buffer=true, desc="Add new empty column to the left"})
vim.keymap.set('n', '<leader><leader>}', "<Plug>TableAddColRight", { buffer=true, desc="Add new empty column to the right"})
vim.keymap.set("n", "<leader><leader>b", "<Cmd>Telescope bibtex<CR>", { buffer=true, desc="Cite" })
vim.keymap.set('n', '<leader><leader><left>', "<Plug>TableGoLeft", { buffer=true, desc="Goto left cell"})
vim.keymap.set('n', '<leader><leader><right>', "<Plug>TableGoRight", { buffer=true, desc="Goto right cell"})
vim.keymap.set('n', '<leader><leader><up>', "<Plug>TableGoUp", { buffer=true, desc="Goto up cell"})
vim.keymap.set('n', '<leader><leader><down>', "<Plug>TableGoDown", { buffer=true, desc="Goto down cell"})
set_keymap_desc('n', '<leader><leader>c', "Clean")
set_keymap_desc('n', '<leader><leader>C', "Clean all")
set_keymap_desc('n', '<leader><leader>e', "Errors")
-- hacky. VimtexErrors puts errors found by Vimtex in quickfix (should be
-- running, use <leader>Lb) then cclose closes quickfix, and then Telescope
-- opens the quickfix in a nicer view.
vim.keymap.set('n', '<space>E', "<Cmd>VimtexErrors<CR>|:cclose|<Cmd>Telescope quickfix<CR>", { buffer=true, desc="Errors"})
set_keymap_desc('n', '<leader><leader>g', "Status")
set_keymap_desc('n', '<leader><leader>G', "Status all")
set_keymap_desc('n', '<leader><leader>i', "Info")
set_keymap_desc('n', '<leader><leader>I', "Info full")
set_keymap_desc('n', '<leader><leader>q', "Log")
-- j for jump. Not using p for preamble since I want to use it for pasting tables.
vim.keymap.set("n", "<leader><leader>j", "<Plug>TableJumpPre", { buffer=true, desc="goto/from preamble (table)" })
set_keymap_desc('n', '<leader><leader>s', "Toggle main")
set_keymap_desc('n', '<leader><leader>t', "TOC open")
set_keymap_desc('n', '<leader><leader>T', "TOC toggle")
vim.keymap.set({"n", "x"}, "<leader><leader>u", "<Plug>Latex2Unicode", { buffer=true, desc="TeX -> unicode" })
vim.keymap.set({"n", "x"}, "<leader><leader>U", "<Plug>Unicode2Latex", { buffer=true, desc="Unicode -> TeX" })
vim.keymap.set('n', '<leader>cv', '<Plug>(vimtex-view)', {buffer=true, desc="View"})
set_keymap_desc('n', '<leader><leader>x', "Reload")
set_keymap_desc('n', '<leader><leader>X', "Reload state")
vim.keymap.set("n", "<leader><leader>y", "<plug>TableYank", { buffer=true, desc="Yank as TSV" })
vim.keymap.set('n', '<leader>ck', '<Plug>(vimtex-stop)', {buffer=true, desc="Stop"})
vim.keymap.set('n', '<leader>cK', '<Plug>(vimtex-stop-all)', {buffer=true, desc="Stop all"})
vim.keymap.set('n', '<leader>cc', '<Plug>(vimtex-compile)', {buffer=true, desc="Compile"})
vim.keymap.set('x', '<leader>cc', '<Plug>(vimtex-compile-selected)', {buffer=true, desc="Compile selected"})
-- TODO: maybe take inspo from these insert mode mappings and make 
-- snippet equivalents then disable them? Or use them if they are useful.
set_keymap_desc('n', '<leader><leader>m', "imaps list")
vim.keymap.set('n', '<leader>cl', '<Plug>(vimtex-compile-output)', {buffer=true, desc="Output"})

-- e.g. \section*{}
set_keymap_desc('n', 'tsc', "Cmd/Star")
set_keymap_desc('n', 'tse', "Env/Star")
-- e.g. with(out) \left 
set_keymap_desc('n', 'tsd', "Delim")
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
vim.keymap.set("n", "ts$", "<Plug>(vimtex-env-toggle-math)", { buffer=true, desc="Inline <-> display" })
vim.keymap.set("n", "ts4", "<Plug>(vimtex-env-toggle-math)", { buffer=true, desc="Inline <-> display" })
vim.keymap.set("n", "tsm", "<plug>(vimtex-env-toggle-math)", { buffer=true, desc="Inline <-> display"})
vim.keymap.set("n", "dsm", "<plug>(vimtex-env-delete-math)", { buffer=true, desc="Delete math"})
vim.keymap.set("n", "csm", "<plug>(vimtex-env-change-math)", { buffer=true, desc="Change math"})
vim.keymap.set("n", "xad", "yaddad", {remap=true, buffer=true, desc="Cut a delim"})
vim.keymap.set("n", "xid", "yiddid", {remap=true, buffer=true, desc="Cut in delim"})
-- item with i instead of m and math with m
vim.keymap.set({"o", "x"}, "ai", "<Plug>(vimtex-am)", {buffer=true, desc="An item"})
vim.keymap.set({"o", "x"}, "ii", "<Plug>(vimtex-im)", {buffer=true, desc="In item"})
vim.keymap.set({"o", "x"}, "am", "<Plug>(vimtex-a$)", {buffer=true, desc="An eq"})
vim.keymap.set({"o", "x"}, "im", "<Plug>(vimtex-i$)", {buffer=true, desc="In eq"})
-- shorthand to $ just using 4 ($ without shift)
vim.keymap.set({"o", "x"}, "a4", "<Plug>(vimtex-a$)", {buffer=true, desc="An eq"})
vim.keymap.set({"o", "x"}, "i4", "<Plug>(vimtex-i$)", {buffer=true, desc="In eq"})
-- next/prev start/end of ...
vim.keymap.set("n", "[m", "<Plug>(vimtex-[n)", { buffer=true, desc="Math start" })
vim.keymap.set("n", "[M", "<Plug>(vimtex-[N)", { buffer=true, desc="Math end" })
vim.keymap.set("n", "[4", "<Plug>(vimtex-[n)", { buffer=true, desc="Math start" })
vim.keymap.set("n", "[$", "<Plug>(vimtex-[N)", { buffer=true, desc="Math end" })
vim.keymap.set("n", "]m", "<Plug>(vimtex-]n)", { buffer=true, desc="Math start" })
vim.keymap.set("n", "]M", "<Plug>(vimtex-]N)", { buffer=true, desc="Math end" })
vim.keymap.set("n", "]4", "<Plug>(vimtex-]n)", { buffer=true, desc="Math start" })
vim.keymap.set("n", "]$", "<Plug>(vimtex-]N)", { buffer=true, desc="Math end" })



-- colorscheme aucmd to fix missing or inconsistent hl links
local grp = vim.api.nvim_create_augroup("Tex", {clear=true})
vim.api.nvim_create_autocmd("Colorscheme", {
    buffer = 0,
    group = grp,
    callback = function ()
        vim.api.nvim_set_hl(0, "texCmd",         {link="@function"})
        vim.api.nvim_set_hl(0, "texCmdEnv",      {link="@keyword.function"}) -- italic instead of bold for begin end
        vim.api.nvim_set_hl(0, "texEnvArgName",  {link="@method"}) -- bold and shine instead of nothing
        vim.api.nvim_set_hl(0, "texCmdBeamer",   {link="@function"})
        vim.api.nvim_set_hl(0, "texOpt",         {link="@parameter"})
        vim.api.nvim_set_hl(0, "texBeamerOpt",   {link="@parameter"})
        vim.api.nvim_set_hl(0, "texOptEqual",    {link="@operator"})
        vim.api.nvim_set_hl(0, "texArg",         {link="@parameter"})
        vim.api.nvim_set_hl(0, "texFileArg",     {link="@string"})
        vim.api.nvim_set_hl(0, "texFilesArg",    {link="@string"})
        vim.api.nvim_set_hl(0, "texFileOpt",     {link="@parameter"})
        vim.api.nvim_set_hl(0, "TexBeamerDelim", {link="Delimiter"})
        vim.api.nvim_set_hl(0, "superscript",    {link="Type"}) -- like \huge, \normalsize etc
        vim.api.nvim_set_hl(0, "subscript",      {link="Type"}) -- like \huge, \normalsize etc
    end
})

