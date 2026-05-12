-- AI coding agent integrations

return {
    -- Local fork of agentic.nvim — native chat UI via Agent Client Protocol
    {
        "agentic.nvim",
        load = function() end,
        cmd = { "Agentic", "AgenticResume" },
        keys = {
            { "<D-\\>", "<Plug>(agentic-toggle-tab)", desc = "Agent toggle" },
            { "<M-CR>", "<Plug>(agentic-send)", desc = "Agent send motion" },
            { "<M-CR><M-CR>", "<Plug>(agentic-send-line)", desc = "Agent send line" },
            { "<M-CR>", "<Plug>(agentic-send)", desc = "Agent send selection", mode = "v" },
            { "<leader>i", nil, desc = "Intelligence/AI" },
            { "<leader>ii", "<Plug>(agentic-toggle-tab)", desc = "Toggle agent" },
            { "<leader>if", "<Plug>(agentic-open)", desc = "Focus agent" },
            { "<leader>iq", "<Plug>(agentic-close)", desc = "Close agent" },
            { "<leader>in", "<Plug>(agentic-new-session)", desc = "New session" },
            { "<leader>ib", "<Plug>(agentic-add-file)", desc = "Add current buffer" },
            { "<leader>is", "<Plug>(agentic-send)", mode = "v", desc = "Send selection" },
        },
        before = function()
            vim.api.nvim_set_hl(0, "AgenticSearchMatch", { link = "DiagnosticError" })
        end,
        after = function()
            require("agentic").setup {
                debug = false,
                log = true,
                winbar = false,
                headers = {
                    chat = { title = "Claude" },
                },
                auto_scroll = {
                    enabled = true,
                },
                notifications = {
                    bell = true,
                },
                acp_providers = {
                    ["opencode-acp"] = {
                        env = {
                            OPENAI_BASE_URL = "https://litellm.dev.xyme.cloud/v1",
                            -- set in untracked .env and available in apple Passwords.
                            -- Also needed ~/.config/opencode/config.json
                            OPENAI_API_KEY = os.getenv("OPENCODE_LITELLM_API_KEY"),
                        },
                    },
                },
            }

            vim.api.nvim_create_user_command("Agentic", function()
                require("agentic").toggle_tab()
            end, { desc = "Open agentic agent tab" })
        end,
    },
}
