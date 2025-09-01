local map = require "utils/keymap"

return {
    {
        "mfussenegger/nvim-dap",
        init = function ()
            map.n("<leader>Dn", "<Cmd>DapNew<CR>", "New DAP")
            map.n("<leader>Dr", "<Cmd>DapToggleRepl<CR>", "REPL toggle")
            map.n("<leader>Db", "<Cmd>DapToggleBreakpoint<CR>", "Breakpoint toggle")
            map.n("<leader><left>", "<Cmd>DapStepOut<CR>", "Step out")
            map.n("<leader><right>", "<Cmd>DapStepInto<CR>", "Step into")
            map.n("<leader><down>", "<Cmd>DapStepOver<CR>", "Step over")
        end,
    },
    {
        "jay-babu/mason-nvim-dap.nvim",
        opts = {
            ensure_installed = {
                "python"
            },
            handlers = {
                python = function(config)
                    config.configurations[1].cwd = "/Users/cmadsen/Documents/enxyme-flow/src/"
                    require('mason-nvim-dap').default_setup(config)
                end,
            },
        },
        dependencies = {
            "williamboman/mason.nvim",
            "mfussenegger/nvim-dap",
        },
    },
    {
        "theHamsta/nvim-dap-virtual-text",
        config = true,
        init = function ()
            map.n("<leader>Dv", "<Cmd>DapVirtualTextToggle<CR>", "Virtual text toggle")
        end,
    },
    {
        "rcarriga/nvim-dap-ui",
        config = true,
    }
}
