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
    ["markdown-preview.nvim"] = function(ev)
        if not ev.data.active then vim.cmd.packadd("markdown-preview.nvim") end
        vim.fn["mkdp#util#install"]()
    end,
    ["blink-cmp-dictionary"] = function()
        if vim.fn.executable("brew") == 1 then
            vim.system({ "brew", "install", "wordnet" }):wait()
        else
            vim.notify("brew not found, skipping wordnet install for blink-cmp-dictionary", vim.log.levels.WARN)
        end
    end,
    ["vimtex"] = function()
        if vim.fn.executable("brew") == 1 then
            vim.system({ "brew", "install", "pstree" })
        else
            vim.notify("brew not found, skipping pstree install for vimtex", vim.log.levels.WARN)
        end
    end,
    ["math-conceal.nvim"] = function(ev)
        vim.system({ "make", "lua51" }, { cwd = ev.data.path }):wait()
    end,
}

function M.on_changed(ev)
    local name = ev.data.spec.name
    if (ev.data.kind == "install" or ev.data.kind == "update") and hooks[name] then
        hooks[name](ev)
    end
end

return M
