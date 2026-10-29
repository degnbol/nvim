local hi = require "utils/highlights"
local map = require "utils/keymap"
-- AI coding agent integrations

return {
    -- Local fork of agentic.nvim — native chat UI via Agent Client Protocol
    {
        "agentic.nvim",
        load = function() end,
        before = function()
            hi.onColorScheme(function()
                hi.link("AgenticSearchMatch", "DiagnosticError")
                hi.setbg("AgenticUnapprovedCommand", hi.bg("DiffDelete"))
            end)
        end,
        after = function()
            require("agentic").setup {
                debug = false,
                log = true,
                winbar = false,
                headers = {
                    agent = { title = "Claude" },
                },
                follow = {
                    enabled = true,
                },
                notifications = {
                    bell = true,
                },
                permissions = {
                    structured = {
                        blender = {
                            safe_write = { { options = { "background" } }, { options = { "b" } } },
                            ask = {
                                {
                                    options = {
                                        "python",
                                        "python-expr",
                                        "python-text",
                                        "python-console",
                                        "P",
                                        "y",
                                        "enable-autoexec",
                                    },
                                },
                            },
                        },
                        pymol = { safe_write = { {} } },
                    },
                },
                keymaps = {
                    prompts = {
                        ["<localLeader>y"] = "Go ahead.",
                        ["<localLeader>P"] = "Spawn the plan reviewer.",
                        ["<localLeader>C"] = "Spawn the code and docstring reviewers on your implementation.",
                    },
                },
                acp_providers = {
                    ["claude-agent-acp"] = {
                        default_mode = "auto",
                        subagents = {
                            force_background = true,
                        },
                        env = {
                            -- Opt out of the bash-first experiment. Under auto
                            -- and bypass modes the SDK otherwise attaches a
                            -- system reminder telling the agent to read and
                            -- edit files through cat/sed/heredocs instead of
                            -- Read/Edit/Write. The flag is a tri-bool read
                            -- straight off the environment, ahead of the
                            -- cohort assignment.
                            CLAUDE_CODE_THRIFTY_SONIC = "0",
                        },
                    },

                    ["opencode-acp"] = {
                        env = {
                            OPENAI_BASE_URL = "https://litellm.svc.xyme.cloud/",
                            -- set in untracked .env and available in apple Passwords.
                            -- Also needed ~/.config/opencode/config.json
                            OPENAI_API_KEY = os.getenv("OPENCODE_LITELLM_API_KEY"),
                        },
                    },
                },
            }

            map.n("<M-CR>", "<Plug>(agentic-send)", "Agent send motion")
            map.n("<M-CR><M-CR>", "<Plug>(agentic-send-line)", "Agent send line")
            map.v("<M-CR>", "<Plug>(agentic-send)", "Agent send selection")
            map.n("<leader>if", "<Plug>(agentic-open)", "Focus agent")
            map.n("<leader>in", "<Plug>(agentic-new-session)", "New session")
            map.n("<leader>ib", "<Plug>(agentic-add-file)", "Add current buffer")
            map.v("<leader>is", "<Plug>(agentic-send)", "Send selection")
        end,
    },
}
