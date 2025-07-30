local util = require "utils/init"

-- remove o, we want to continue comments while editing them only (r).
-- no t and having c+a means only comments are autoformatted.
-- However, a made comment reformat slow, so don't use by default.
vim.opt.formatoptions = "jwcrql"
vim.opt.concealcursor = ""
vim.opt.list = false

local grp = vim.api.nvim_create_augroup("colorscheme", { clear = true })
vim.api.nvim_create_autocmd("Colorscheme", {
    buffer = 0,
    group = grp,
    callback = function()
        vim.api.nvim_set_hl(0, "@cell", { reverse = true })
    end
})

local function load_pymol()
    -- load additional pymol syntax hl
    local rtp = vim.opt.runtimepath:get()[1]
    vim.schedule(function()
        vim.cmd.source(rtp .. "/syntax/python_pymol.vim")
    end)
    -- The above needs to be done for each new buffer but the following only needs to be run once,
    -- hence the bool check.
    -- This global var is also used by blink to enable pymol_settings provider
    if not vim.g.loaded_pymol then
        vim.g.loaded_pymol = true
        -- load pymol snippets
        require("luasnip").add_snippets("python", require 'luasnippets.python_pymol')
        -- Load pymol settings completion source.
        -- Using blink.cmp instead.
        -- require "completion.pymol.cmp_pymol_settings".setup()
    end
end
-- manually load
vim.keymap.set('n', '<localleader>+', load_pymol, {
    buffer = true,
    desc =
    "Manually load pymol snippets+completion+syntax"
})
-- check if pymol is loaded by scanning first 10 lines
for _, line in ipairs(vim.api.nvim_buf_get_lines(0, 0, 10, false)) do
    -- might be using e.g. `from pymol_util import *`
    if line:match("import") and line:match("pymol") then
        return load_pymol()
    end
end

-- Filter the default goto references so we don't see the "build/" references.
vim.keymap.set('n', 'grr', function()
    vim.lsp.buf.references(nil, {
        on_list = function(options)
            local items = {}
            for _, item in ipairs(options.items) do
                -- :h setqflist-what
                if not item.filename:match("build/") then
                    table.insert(items, item)
                end
            end
            options.items = items
            vim.fn.setqflist({}, ' ', options)
            vim.cmd('botright copen')
        end
    })
end, { desc = "Goto references (excl build/)", buffer = true })
