-- have gf (goto file) work when writing the common $ROOT/PATH pattern.
vim.opt_local.includeexpr = [[substitute(v:fname,'\$ROOT/','','')]]

-- quick macros for toggling between inline and not inline functions.
-- Maybe write full lua functions so we can also support begin/end notation.
-- @i: $[f go to beginning of function we are inside. $ so we don't go to
-- previous function in case we are on first char.
-- @f: Ifunction <Esc>f(%f= — we want to go to the = that defines the function
-- but there may be = inside the function args for default values and the
-- contents of the function could be e.g. arg == something so we find it by
-- finding open paren, jumping to matching paren to jump over args, then first
-- = should be it.
local nvim_code = require("utils/init").nvim_code
vim.fn.setreg("i", nvim_code("$[fdwA = <Esc>JJD"))
vim.fn.setreg("f", nvim_code("Ifunction <Esc>f(%f=caw<BS><CR><Esc>oend<Esc>[f=af"))

vim.cmd.iabbrev("edn", "end")

local map = require "utils/keymap"

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

