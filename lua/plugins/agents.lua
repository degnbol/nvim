-- AI coding agent integrations
-- Toggle `enabled` to switch between plugins (only one should be active)

return {
    -- Claude Code via WebSocket MCP protocol (same as VS Code extension)
    -- https://github.com/coder/claudecode.nvim
    {
        "coder/claudecode.nvim",
        enabled = true,
        dependencies = { "folke/snacks.nvim" },
        keys = function()
            local function send_operator()
                vim.opt.operatorfunc = "v:lua.ClaudeCodeSendOperator"
                return "g@"
            end
            local function send_line()
                vim.opt.operatorfunc = "v:lua.ClaudeCodeSendOperator"
                return "g@_"
            end
            function ClaudeCodeSendOperator(type)
                if type == "char" then
                    vim.cmd('silent normal! `[v`]')
                else
                    vim.cmd('silent normal! `[V`]')
                end
                vim.cmd('ClaudeCodeSend')
            end
            return {
                { "<D-\\>", "<cmd>ClaudeCode<cr>", desc = "Agent toggle" },
                { "<M-CR>", send_operator, desc = "Agent send motion", expr = true },
                { "<M-CR><M-CR>", send_line, desc = "Agent send line", expr = true },
                { "<M-CR>", "<cmd>ClaudeCodeSend<cr>", desc = "Agent send selection", mode = "v" },
                { "<leader>i", nil, desc = "Intelligence/AI" },
                { "<leader>ii", "<cmd>ClaudeCode<cr>", desc = "Toggle agent" },
                { "<leader>if", "<cmd>ClaudeCodeFocus<cr>", desc = "Focus agent" },
                { "<leader>iq", "<cmd>ClaudeCodeClose<cr>", desc = "Close agent" },
                { "<leader>in", "<cmd>ClaudeCode<cr>", desc = "New session" },
                { "<leader>ib", "<cmd>ClaudeCodeAdd %<cr>", desc = "Add current buffer" },
                { "<leader>is", "<cmd>ClaudeCodeSend<cr>", mode = "v", desc = "Send selection" },
                { "<leader>ir", "<cmd>ClaudeCode --resume<cr>", desc = "Resume session" },
                { "<leader>ic", "<cmd>ClaudeCode --continue<cr>", desc = "Continue session" },
                { "<leader>im", "<cmd>ClaudeCodeSelectModel<cr>", desc = "Select model" },
                { "<leader>ia", "<cmd>ClaudeCodeDiffAccept<cr>", desc = "Accept diff" },
                { "<leader>id", "<cmd>ClaudeCodeDiffDeny<cr>", desc = "Deny diff" },
            }
        end,
        opts = {
            terminal_cmd = vim.fn.expand("~/.local/bin/claude"),
            terminal = {
                split_side = "right",
                split_width_percentage = 0.4,
            },
        },
    },
    -- Native chat UI via Agent Client Protocol (multi-provider)
    -- https://github.com/carlos-algms/agentic.nvim
    {
        "carlos-algms/agentic.nvim",
        enabled = false,
        build = "pnpm install",
        dependencies = { "nvim-lua/plenary.nvim" },
        keys = function()
            local agentic = function(fn)
                return function() require("agentic")[fn]() end
            end
            local function send_operator()
                vim.opt.operatorfunc = "v:lua.AgenticSendOperator"
                return "g@"
            end
            local function send_line()
                vim.opt.operatorfunc = "v:lua.AgenticSendOperator"
                return "g@_"
            end
            function AgenticSendOperator(type)
                if type == "char" then
                    vim.cmd('silent normal! `[v`]')
                else
                    vim.cmd('silent normal! `[V`]')
                end
                require("agentic").add_selection()
            end
            return {
                { "<D-\\>", agentic("toggle"), desc = "Agent toggle" },
                { "<M-CR>", send_operator, desc = "Agent send motion", expr = true },
                { "<M-CR><M-CR>", send_line, desc = "Agent send line", expr = true },
                { "<M-CR>", agentic("add_selection"), desc = "Agent send selection", mode = "v" },
                { "<leader>i", nil, desc = "Intelligence/AI" },
                { "<leader>ii", agentic("toggle"), desc = "Toggle agent" },
                { "<leader>if", agentic("open"), desc = "Focus agent" },
                { "<leader>iq", agentic("close"), desc = "Close agent" },
                { "<leader>in", agentic("new_session"), desc = "New session" },
                { "<leader>ix", agentic("stop_generation"), desc = "Stop generation" },
                { "<leader>ib", agentic("add_file"), desc = "Add current buffer" },
                { "<leader>is", agentic("add_selection"), mode = "v", desc = "Send selection" },
            }
        end,
        opts = {
            provider = "claude-acp",
        },
    },
}
