local map = require "utils/keymap"

return {
    {
        "nvim-dap",
        before = function ()
            map.n("<leader>Dn", "<Cmd>DapNew<CR>", "New DAP")
            map.n("<leader>Dr", "<Cmd>DapToggleRepl<CR>", "REPL toggle")
            map.n("<leader>Db", "<Cmd>DapToggleBreakpoint<CR>", "Breakpoint toggle")
            map.n("<leader><left>", "<Cmd>DapStepOut<CR>", "Step out")
            map.n("<leader><right>", "<Cmd>DapStepInto<CR>", "Step into")
            map.n("<leader><down>", "<Cmd>DapStepOver<CR>", "Step over")
        end,
    },
    {
        "mason-nvim-dap.nvim",
        before = function()
            require("lz.n").trigger_load("mason.nvim")
            vim.cmd.packadd("nvim-dap")
        end,
        after = function()
            require("mason-nvim-dap").setup {
                ensure_installed = {
                    "python"
                },
                handlers = {
                    python = function(config)
                        config.configurations[1].cwd = "/Users/cmadsen/Documents/enxyme-flow/src/"
                        require('mason-nvim-dap').default_setup(config)
                    end,
                },
            }
        end,
    },
    {
        "nvim-dap-virtual-text",
        after = function()
            vim.cmd.packadd("nvim-dap")
            require("nvim-dap-virtual-text").setup()
        end,
        before = function ()
            map.n("<leader>Dv", "<Cmd>DapVirtualTextToggle<CR>", "Virtual text toggle")
        end,
    },
    {
        "nvim-dap-ui",
        after = function()
            vim.cmd.packadd("nvim-nio")
            vim.cmd.packadd("nvim-dap")
            require("dapui").setup()
        end,
    }
}
