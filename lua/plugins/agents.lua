-- AI coding agent integrations

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

return {
    -- Local fork of agentic.nvim — native chat UI via Agent Client Protocol
    {
        "agentic.nvim",
        enabled = true,
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
            }

            vim.api.nvim_create_user_command("Agentic", function()
                require("agentic").toggle_tab()
            end, { desc = "Open agentic agent tab" })

            vim.api.nvim_create_user_command("AgenticResume", function(args)
                local full_id, cwd = resolve_session_prefix(args.args)
                if not full_id then
                    vim.notify("No session found matching: " .. args.args, vim.log.levels.ERROR)
                    return
                end
                require("agentic").toggle_tab()
                require("agentic").load_acp_session(full_id, cwd)
            end, { nargs = 1, desc = "Resume agentic session by UUID prefix" })
        end,
    },
}
