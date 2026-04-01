vim.opt.list = false
-- to wrap comments using gw (default vim formatting)
vim.opt.formatoptions:append("t")

-- lz.n spec snippets for lua/plugins/ files (deferred until LuaSnip loads)
if vim.api.nvim_buf_get_name(0):match("lua/plugins/") and not vim.g.lzn_snippets then
    vim.g.lzn_snippets = true
    if package.loaded["luasnip"] then
        require("luasnip").add_snippets("lua", require 'luasnippets.lzn')
    else
        vim.api.nvim_create_autocmd("InsertEnter", {
            once = true,
            callback = function()
                vim.schedule(function()
                    if package.loaded["luasnip"] then
                        require("luasnip").add_snippets("lua", require 'luasnippets.lzn')
                    end
                end)
            end,
        })
    end
end

-- use :help instead of default Man for visual selection in lua
vim.keymap.set('x', 'K', [["hy:h <C-r>h<CR>]], { buffer=true, desc="Help" })

-- when doing gf or similar obviously we should look in the lua/ folder since 
-- this is where all scripts are required from.
local rtp = vim.opt.runtimepath:get()[1]
vim.opt_local.path:append(rtp .. "/lua")

vim.cmd.iabbrev("ture", "true")
vim.cmd.iabbrev("flase", "false")

