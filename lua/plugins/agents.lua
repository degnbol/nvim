-- AI coding agent integrations
-- Toggle `enabled` to switch between plugins (only one should be active)

-- Extract cwd from the first few lines of a JSONL session file.
-- Returns the cwd string or nil.
local function extract_cwd(jsonl_path)
    local f = io.open(jsonl_path, "r")
    if not f then return nil end
    for _ = 1, 10 do
        local line = f:read("*l")
        if not line then break end
        local cwd = line:match('"cwd":"([^"]+)"')
        if cwd then
            f:close()
            return cwd
        end
    end
    f:close()
    return nil
end

-- Resolve a UUID prefix to a full UUID by scanning Claude session files.
-- Returns (full_uuid, cwd) or nil if not found / ambiguous.
local function resolve_session_prefix(prefix)
    local claude_dir = vim.fn.expand("~/.claude/projects")
    if vim.fn.isdirectory(claude_dir) == 0 then return nil end
    local matches = {}
    for project_dir, dtype in vim.fs.dir(claude_dir) do
        if dtype == "directory" then
            local dir_path = vim.fs.joinpath(claude_dir, project_dir)
            for file, ftype in vim.fs.dir(dir_path) do
                if ftype == "file" and file:match("%.jsonl$") then
                    local uuid = file:gsub("%.jsonl$", "")
                    if uuid:sub(1, #prefix) == prefix then
                        table.insert(matches, {
                            uuid = uuid,
                            path = vim.fs.joinpath(dir_path, file),
                        })
                    end
                end
            end
        end
    end
    if #matches == 1 then
        return matches[1].uuid, extract_cwd(matches[1].path)
    elseif #matches > 1 then
        vim.notify("Ambiguous session prefix: " .. #matches .. " matches", vim.log.levels.WARN)
    end
    return nil
end

-- Open agentic and resume an ACP session by UUID prefix.
-- Usage: `nvim -c 'lua AgenticResume("8583a113")'`
function AgenticResume(prefix)
    local full_id, cwd = resolve_session_prefix(prefix)
    if not full_id then
        vim.notify("No session found matching: " .. prefix, vim.log.levels.ERROR)
        return
    end

    require("agentic").toggle_tab()
    -- Call immediately — load_acp_session defers internally if the agent
    -- isn't ready yet, avoiding the race where new_session() fires first.
    require("agentic").load_acp_session(full_id, cwd)
end

return {
    -- Claude Code via WebSocket MCP protocol (same as VS Code extension)
    -- https://github.com/coder/claudecode.nvim
    {
        "claudecode.nvim",
        enabled = false,
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
        after = function()
            require("claudecode").setup {
                terminal_cmd = vim.fn.expand("~/.local/bin/claude"),
                terminal = {
                    split_side = "right",
                    split_width_percentage = 0.4,
                },
            }
        end,
    },
    -- Local fork of agentic.nvim — native chat UI via Agent Client Protocol
    {
        "agentic.nvim",
        enabled = true,
        load = function() end,
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
        after = function()
            require("agentic").setup {
                headers = {
                    chat = { title = "Claude" },
                },
                notifications = {
                    bell = true,
                },
            }
        end,
        before = function()
            -- Red text like ripgrep's default match colour.
            -- terminal_color_1 (ANSI red) is nil — colorscheme doesn't set it.
            vim.api.nvim_set_hl(0, "AgenticSearchMatch", { link = "DiagnosticError" })
        end,
    },
}
