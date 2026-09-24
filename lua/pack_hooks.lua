local util = require "utils/init"

local M = {}

---Run a command to completion, notifying if it fails.
---@param cmd string[]
---@param cwd string|nil nil runs in the current directory
local function run(cmd, cwd)
    util.notify_failure(cmd, vim.system(cmd, { cwd = cwd, text = true }):wait())
end

local hooks = {
    ["mason.nvim"] = function()
        vim.cmd.packadd("mason.nvim")
        require("mason").setup()
        vim.cmd("MasonUpdate")
    end,
    ["nvim-treesitter"] = function()
        vim.cmd.packadd("nvim-treesitter")
        require("nvim-treesitter").update()
    end,
    ["LuaSnip"] = function(ev)
        run({ "make", "install_jsregexp" }, ev.data.path)
    end,
    ["markdown-preview.nvim"] = function(ev)
        if not ev.data.active then vim.cmd.packadd("markdown-preview.nvim") end
        vim.fn["mkdp#util#install"]()
    end,
    ["blink-cmp-dictionary"] = function()
        if vim.fn.executable("brew") == 1 then
            run({ "brew", "install", "wordnet" })
        else
            vim.notify("brew not found, skipping wordnet install for blink-cmp-dictionary", vim.log.levels.WARN)
        end
    end,
    ["vimtex"] = function()
        if vim.fn.executable("brew") == 1 then
            local cmd = { "brew", "install", "pstree" }
            vim.system(cmd, { text = true }, function(obj) util.notify_failure(cmd, obj) end)
        else
            vim.notify("brew not found, skipping pstree install for vimtex", vim.log.levels.WARN)
        end
    end,
    ["math-conceal.nvim"] = function(ev)
        run({ "make", "lua51" }, ev.data.path)
    end,
}

function M.on_changed(ev)
    local name = ev.data.spec.name
    if (ev.data.kind == "install" or ev.data.kind == "update") and hooks[name] then
        hooks[name](ev)
    end
end

return M
