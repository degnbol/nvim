local util = require "utils/init"
local md_table = require "utils.markdown_table"

local opt = vim.opt_local
-- markdown doesn't really have comments. 
-- The default is html comments <!-- ... -->.
-- For making tex documents with the markdown package it is possible to ignore latex comments, hence:
opt.commentstring = "% %s"
-- c=continue comment when wrapping. Should already be enabled.
-- r=with enter in insert mode. Should already be enabled.
-- o=with o/O in normal mode. Not sure if I want this.
-- a=auto format by default
-- v=inserted now, don't edit other lines.
opt.formatoptions:append "crv"
-- this also had to be set for the above to understand how a comment line looks.
-- the b means the % should always be followed by a blank, i.e. space
opt.comments = "b:%"
-- ai is necessary for lists but a problem if I try to indent start of sentence, then the next lines are also indented.
-- opt.autoindent = false
-- since we are actively using the textwidth in markdown then the moving of the window is just annoying for a slim window
opt.sidescrolloff = 0

-- conceal comment leader and the _ and * around emphasis.
-- Editing a file: level 1 conceals but keeps block markers visible (e.g. the
-- backslash of escapes stays — see after/queries/markdown_inline). Read-only
-- scratch markdown (LSP hover floats, file peeks): level 3 fully hides, so
-- escape backslashes disappear for a clean rendered look. Core sets floats to
-- 2 before setting filetype, so this ftplugin (run last) sets the final value.
opt.conceallevel = vim.bo.buftype == "" and 1 or 3

-- mini.hipatterns auto-enables only on normal buffers (buftype == ""), so
-- read-only scratch markdown (LSP hover floats, peeks) gets no hex swatch.
-- Enable it manually so colour definitions in hover floats get the same
-- end-of-line swatch as source buffers.
if vim.bo.buftype ~= "" then
    require("mini.hipatterns").enable(0)
end

opt.list = false

vim.keymap.set("n", "<localleader>a", function()
    local row = util.get_cursor()
    if not md_table.align_at(0, row, { max_width = 30, align_cells = true }) then
        vim.notify("not on a markdown table", vim.log.levels.INFO)
    end
end, { buffer = true, desc = "Align markdown table" })
