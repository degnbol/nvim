-- have gf (goto file) work when writing the common $ROOT/PATH pattern.
vim.opt_local.includeexpr = [[substitute(v:fname,'\$ROOT/','','')]]

vim.cmd.iabbrev("edn", "end")

local map = require "utils/keymap"

-- Toggle a function definition between short (f(x) = …) and long form; conform
-- an overload set to the cursor def's current form.
map.buf('n', '<localleader>f', require("julia_fndef").toggle, "Toggle fn def form")
map.buf('n', '<localleader>F', require("julia_fndef").conform, "Conform fn def forms")

-- remove o, we want to continue comments while editing them only (r).
vim.opt_local.formatoptions = "jwcrql"

-- you have to be in insert mode to unconceal
vim.opt_local.concealcursor = "nvc"

vim.opt_local.list = false

map.buf('n', '<leader>cc', '<Cmd>!julia %<CR>', "Run this script")

-- using blink.cmp instead
-- require"completion.plotlyjs.cmp_plotlyjs".setup()

-- manually mark that plotlyjs is being used
map.buf('n', '<localleader>+', function() vim.g.loaded_plotly = true end, "Manually load plotlyjs completion")
-- ... or check if plotlyjs is loaded by scanning first 20 lines
for _, line in ipairs(vim.api.nvim_buf_get_lines(0, 0, 20, false)) do
    if line:lower():match("plotly") then
        vim.g.loaded_plotly = true
    end
end

-- Manual efforts. Install julia LSP as described on
-- https://github.com/neovim/nvim-lspconfig/blob/master/lsp/julials.lua
vim.lsp.enable("julials")
-- vim.lsp.enable("jetls")

-- Annoying bug where something keeps setting buftype=nofile
vim.api.nvim_create_autocmd("BufWritePre", {
    buffer = 0,
    group = vim.api.nvim_create_augroup("julia_bug", { clear = true }),
    callback = function() vim.bo.buftype = nil end
})

-- TEMP: there should be a better way to do this. syn must be getting reset somewhere.
vim.defer_fn(function ()
    vim.cmd [[syn keyword @variable.builtin stdin stdout stderr]]
end, 0)

