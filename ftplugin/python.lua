local map = require "utils/keymap"

-- remove o, we want to continue comments while editing them only (r).
-- no t and having c+a means only comments are autoformatted.
-- However, a made comment reformat slow, so don't use by default.
vim.opt.formatoptions = "jwcrql"
vim.opt.concealcursor = ""
vim.opt.list = false

map.buf('n', '<leader>cc', '<Cmd>!python %<CR>', "Run this script")

local hi = require "utils/highlights"

local function set_pymol_hl()
    -- Off-white fg from Normal, so pymol keywords show through green @string
    hi.set("@variable.builtin.pymol_select", { fg = hi.fg("Normal"), italic = true })
end

local grp = vim.api.nvim_create_augroup("colorscheme", { clear = true })
vim.api.nvim_create_autocmd("Colorscheme", {
    buffer = 0,
    group = grp,
    callback = function()
        hi.set("@cell", { reverse = true })
        if vim.g.loaded_pymol then set_pymol_hl() end
    end
})

local function load_pymol()
    -- Set up treesitter injection for pymol_select in Python strings.
    -- Scoped to function args, keyword args, and assignments — not docstrings.
    local read_query = require('utils/init').read_query
    local base = read_query('python', 'injections')
    local pymol_inject = read_query('pymol_select', 'python_injections')
    vim.treesitter.query.set('python', 'injections', base .. '\n' .. pymol_inject)
    set_pymol_hl()

    -- Force treesitter to re-evaluate injections with the new query
    local ok, parser = pcall(vim.treesitter.get_parser, 0)
    if ok and parser then
        parser:invalidate(true)
    end

    -- This global var is also used by blink to enable pymol_settings provider
    if not vim.g.loaded_pymol then
        vim.g.loaded_pymol = true
        require("luasnip").add_snippets("python", require 'luasnippets.python_pymol')
    end
end
-- manually load
map.n('<localleader>+', load_pymol, "Manually load pymol snippets+completion+syntax", { buffer = true, })
-- check if pymol is loaded by scanning first 10 lines
for _, line in ipairs(vim.api.nvim_buf_get_lines(0, 0, 10, false)) do
    -- might be using e.g. `from pymol_util import *`
    if line:match("import") and line:match("pymol") then
        return load_pymol()
    end
end

-- Filter the default goto references so we don't see
-- - "build/" references,
-- - The line we are calling from,
-- - Import statements.
map.n('grr', function()
    vim.lsp.buf.references(nil, map.filter_lsp_items(function(item)
        return not (
            map.qf_item_is_self(item) or
            item.filename:match("build/") or
            item.text:match("^import") or
            item.text:match("^from .* import")
        )
    end))
end, "Goto filtered references", { buffer = true })
