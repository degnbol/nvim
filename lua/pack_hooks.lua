local M = {}

local hooks = {
    ["mason.nvim"] = function()
        vim.cmd.packadd("mason.nvim")
        require("mason").setup()
        vim.cmd("MasonUpdate")
    end,
    ["nvim-treesitter"] = function()
        vim.cmd.packadd("nvim-treesitter")
        pcall(vim.cmd, "TSUpdate")
    end,
    ["LuaSnip"] = function(ev)
        vim.system({ "make", "install_jsregexp" }, { cwd = ev.data.path }):wait()
    end,
    ["telescope-fzf-native.nvim"] = function(ev)
        vim.system({ "make" }, { cwd = ev.data.path }):wait()
    end,
    ["markdown-preview.nvim"] = function(ev)
        if not ev.data.active then vim.cmd.packadd("markdown-preview.nvim") end
        vim.fn["mkdp#util#install"]()
    end,
    ["blink-cmp-dictionary"] = function()
        vim.system({ "brew", "install", "wordnet" }):wait()
    end,
    ["vimtex"] = function()
        vim.system({ "brew", "install", "pstree" })
    end,
    ["math-conceal.nvim"] = function(ev)
        vim.system({ "make", "lua51" }, { cwd = ev.data.path }):wait()
    end,
    ["fff.nvim"] = function(ev)
        vim.system({ "cargo", "build", "--release" }, { cwd = ev.data.path }):wait()
    end,
}

function M.on_changed(ev)
    local name = ev.data.spec.name
    if (ev.data.kind == "install" or ev.data.kind == "update") and hooks[name] then
        hooks[name](ev)
    end
end

return M
