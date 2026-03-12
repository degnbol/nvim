-- AI coding agent integrations
-- Toggle `enabled` to switch between plugins (only one should be active)

-- Open agentic in a dedicated tab; close tab when toggling off.
-- Exposed globally so `nvim -c 'lua AgenticToggle()'` works from the shell.
function AgenticToggle()
    local a = require("agentic")
    local reg = require("agentic.session_registry")
    local tab = vim.api.nvim_get_current_tabpage()
    local session = reg.sessions[tab]
    if session and session.widget:is_open() then
        session.widget:hide()
        -- Close the tab if it was created for agentic (only has 1 window left)
        if #vim.api.nvim_tabpage_list_wins(tab) <= 1 then
            vim.cmd("tabclose")
        end
    else
        -- Reuse current tab if it's just a dashboard or empty buffer
        local wins = vim.fn.filter(
            vim.api.nvim_tabpage_list_wins(tab),
            function(_, w) return vim.api.nvim_win_get_config(w).relative == "" end
        )
        local buf = vim.api.nvim_win_get_buf(wins[1])
        local ft = vim.bo[buf].filetype
        local fresh = #wins == 1 and (ft == "dashboard" or (ft == "" and vim.fn.bufname(buf) == ""))
        if not fresh then
            vim.cmd("tabnew")
        end
        local empty_win = vim.api.nvim_get_current_win()
        a.open({ auto_add_to_context = false })
        -- Close the leftover empty window so agentic fills the tab
        if vim.api.nvim_win_is_valid(empty_win) then
            local ebuf = vim.api.nvim_win_get_buf(empty_win)
            pcall(vim.api.nvim_win_close, empty_win, true)
            pcall(vim.api.nvim_buf_delete, ebuf, { force = true })
        end
    end
end

return {
    -- Claude Code via WebSocket MCP protocol (same as VS Code extension)
    -- https://github.com/coder/claudecode.nvim
    {
        "coder/claudecode.nvim",
        enabled = false,
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
            _G.ClaudeCodeSendOperator = function(type)
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
    -- Local fork of agentic.nvim — native chat UI via Agent Client Protocol
    {
        dir = vim.fn.stdpath("config") .. "/agentic.nvim",
        name = "agentic.nvim",
        enabled = true,
        build = "npm i -g @zed-industries/claude-agent-acp",
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
            _G.AgenticSendOperator = function(type)
                if type == "char" then
                    vim.cmd('silent normal! `[v`]')
                else
                    vim.cmd('silent normal! `[V`]')
                end
                require("agentic").add_selection()
            end
            return {
                { "<D-\\>", AgenticToggle, desc = "Agent toggle" },
                { "<M-CR>", send_operator, desc = "Agent send motion", expr = true },
                { "<M-CR><M-CR>", send_line, desc = "Agent send line", expr = true },
                { "<M-CR>", agentic("add_selection"), desc = "Agent send selection", mode = "v" },
                { "<leader>i", nil, desc = "Intelligence/AI" },
                { "<leader>ii", AgenticToggle, desc = "Toggle agent" },
                { "<leader>if", agentic("open"), desc = "Focus agent" },
                { "<leader>iq", agentic("close"), desc = "Close agent" },
                { "<leader>in", agentic("new_session"), desc = "New session" },
                { "<leader>ix", agentic("stop_generation"), desc = "Stop generation" },
                { "<leader>ib", agentic("add_file"), desc = "Add current buffer" },
                { "<leader>is", agentic("add_selection"), mode = "v", desc = "Send selection" },
            }
        end,
        opts = {},
    },
}
