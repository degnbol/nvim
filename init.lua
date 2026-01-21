-- Check for large files before plugins load
require("largefile").check_argv()

-- map leaders has to be set before running lazy
require "options"

-- bootstrap lazy.nvim
local lazypath = vim.fn.stdpath("data") .. "/lazy/lazy.nvim"
if not vim.uv.fs_stat(lazypath) then
    vim.fn.system { "git", "clone", "--filter=blob:none",
        "https://github.com/folke/lazy.nvim.git", "--branch=stable", lazypath, }
end
local rtp = vim.opt.runtimepath:get()[1]
vim.opt.runtimepath:prepend(lazypath)
require 'lazy'.setup('plugins', {
    -- https://github.com/folke/lazy.nvim
    change_detection = { notify = false, },
    dev = {
        -- Directory with local plugin projects.
        path = rtp,
        ---@type string[] plugins that match these patterns will use your local versions instead of being fetched from GitHub
        patterns = { "degnbol" },
        fallback = true,                       -- Fallback to git when local plugin doesn't exist
    },
    install = { colorscheme = { "terafox" } }, -- Don't switch colorscheme after startup install
})

