require "paths"
require "options"

-- bootstrap lazy.nvim
local lazypath = vim.fn.stdpath("data") .. "/lazy/lazy.nvim"
if not vim.loop.fs_stat(lazypath) then
  vim.fn.system {"git", "clone", "--filter=blob:none",
    "https://github.com/folke/lazy.nvim.git", "--branch=stable", lazypath, }
end
vim.opt.rtp:prepend(lazypath)
require'lazy'.setup('plugins', {
    -- https://github.com/folke/lazy.nvim
    change_detection = { notify = false, },
    dev = {
        -- directory where you store your local plugin projects
        path = "$XDG_CONFIG_HOME/nvim",
        ---@type string[] plugins that match these patterns will use your local versions instead of being fetched from GitHub
        patterns = {"degnbol"},
        fallback = true, -- Fallback to git when local plugin doesn't exist
    },
    install = {colorscheme = {"default"}}, -- Don't switch colorscheme after startup install
})

require "highlights"
require "keymap"
require "blockim"
require "autocmds"

