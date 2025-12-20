local util = require "utils/init"
local map = require "utils/keymap"

-- assume using typst.vim
map.n('<leader>cc', "<Cmd>TypstWatch<CR>", "Compile continuously", { buffer=true, })
local path_pdf = vim.api.nvim_buf_get_name(0):gsub(".typ", ".pdf")
if util.is_mac() then
    map.n('<leader>cv', "!open -a skim " .. path_pdf .. "<CR>", "Compile view", { buffer=true, silent=true, })
else
    map.n('<leader>cv', "!pdf " .. path_pdf .. "<CR>", "Compile view", { buffer=true, silent=true, })
end

-- toggle updating pdf on typing.
-- Doesn't work if running module file, hence we use TypstWatch instead by default.
map.n('<LocalLeader>t', function ()
    local lsp = require"lspconfig".typst_lsp
    local exportPdf
    if lsp.manager.config.settings.exportPdf == "onType" then
        exportPdf = "never" -- use TypstWatch to compile
    else
        exportPdf = "onType"
    end
    vim.cmd.LspStop()
    require"lspconfig".typst_lsp.setup { settings = { exportPdf = exportPdf } }
    vim.cmd.LspStart()
    print(exportPdf)
end, "Toggle update on type", { buffer=true})

vim.opt_local.wrap = true
vim.opt_local.sidescrolloff = 0
-- definitely don't break at @ in typst (when wrapping with linebreak). 
-- It is used right before citations and is not a place to break.
vim.opt_local.breakat:remove('@')
