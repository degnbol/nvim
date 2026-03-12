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
    -- Native chat UI via Agent Client Protocol (multi-provider)
    -- https://github.com/carlos-algms/agentic.nvim
    {
        "carlos-algms/agentic.nvim",
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
        opts = {
            keymaps = {
                prompt = {
                    submit = { "<CR>" }, -- removed <C-s>; :w submits instead (see config)
                },
            },
            provider = "claude-agent-acp",
            windows = {
                width = "50%",
                chat = { win_opts = {
                    -- wrap = false, -- Do actually need wrap for long text outputs.
                } },
            },
            headers = {
                input = function() end,
                files = function() end,
                chat = function() end,
                code = { suffix = "" },
                diagnostics = { suffix = "" },
            },
        },
        config = function(_, opts)
            require("agentic").setup(opts)
            -- Show diffs in a separate tab instead of splitting the agentic tab
            local DiffPreview = require("agentic.ui.diff_preview")
            local AgConfig = require("agentic.config")
            local SM0 = require("agentic.session_manager")
            function SM0:_show_diff_in_buffer(tool_call_id)
                if not AgConfig.diff_preview.enabled then return end
                local tracker = tool_call_id
                    and self.message_writer.tool_call_blocks[tool_call_id]
                if not tracker or tracker.kind ~= "edit" or tracker.diff == nil then
                    return
                end
                DiffPreview.show_diff({
                    file_path = tracker.argument,
                    diff = tracker.diff,
                    get_winid = function(bufnr)
                        local agent_tab = vim.api.nvim_get_current_tabpage()
                        -- Reuse existing diff tab or create one
                        local ok, diff_tab = pcall(vim.api.nvim_tabpage_get_var, agent_tab, "_agentic_diff_tab")
                        if ok and vim.api.nvim_tabpage_is_valid(diff_tab) then
                            vim.api.nvim_set_current_tabpage(diff_tab)
                        else
                            vim.cmd("tabnew")
                            diff_tab = vim.api.nvim_get_current_tabpage()
                            vim.api.nvim_tabpage_set_var(agent_tab, "_agentic_diff_tab", diff_tab)
                            vim.api.nvim_set_current_tabpage(agent_tab)
                        end
                        local winid = vim.api.nvim_tabpage_list_wins(diff_tab)[1]
                        vim.api.nvim_win_set_buf(winid, bufnr)
                        return winid
                    end,
                })
            end
            -- Subtle status: just the state text, no spinner animation
            local SA = require("agentic.ui.status_animation")
            local NS = vim.api.nvim_create_namespace("agentic_animation")
            function SA:_render_frame()
                if not self._state or not vim.api.nvim_buf_is_valid(self._bufnr) then return end
                local lines = vim.api.nvim_buf_get_lines(self._bufnr, 0, -1, false)
                local line = math.max(0, #lines - 1)
                self._extmark_id = vim.api.nvim_buf_set_extmark(self._bufnr, NS, line, 0, {
                    id = self._extmark_id,
                    virt_lines = { { { " " .. self._state, "NonText" } } },
                    virt_lines_above = false,
                })
                -- No timer — static text, no animation loop
            end
            -- :w submits the prompt in the input buffer
            vim.api.nvim_create_autocmd("FileType", {
                pattern = "AgenticInput",
                callback = function(ev)
                    -- Schedule: _create_new_buf sets options in unordered loop,
                    -- so buftype=nofile may be set after filetype triggers this autocmd.
                    vim.schedule(function()
                        vim.api.nvim_buf_set_name(ev.buf, "agentic://prompt")
                        vim.bo[ev.buf].buftype = "acwrite"
                    end)
                    vim.api.nvim_create_autocmd("BufWriteCmd", {
                        buffer = ev.buf,
                        callback = function()
                            vim.bo[ev.buf].modified = false
                            vim.api.nvim_feedkeys(
                                vim.api.nvim_replace_termcodes("<CR>", true, false, true),
                                "m", false
                            )
                        end,
                    })
                end,
            })
            -- Defer prompt until agent is ready and session exists (avoids error on early submit)
            local SM = require("agentic.session_manager")
            local orig_submit = SM._handle_input_submit
            function SM:_handle_input_submit(...)
                if self.session_id and self.agent and self.agent.state == "ready" then
                    return orig_submit(self, ...)
                end
                local args = { ... }
                local timer = vim.uv.new_timer()
                local attempts = 0
                timer:start(100, 100, vim.schedule_wrap(function()
                    attempts = attempts + 1
                    if self.session_id and self.agent and self.agent.state == "ready" then
                        timer:stop()
                        timer:close()
                        orig_submit(self, unpack(args))
                    elseif attempts >= 100
                        or (self.agent and (self.agent.state == "error" or self.agent.state == "disconnected"))
                    then
                        timer:stop()
                        timer:close()
                    end
                end))
            end
            -- Mark user prompts with "❯" sign and [/] navigation
            local PROMPT_NS = vim.api.nvim_create_namespace("agentic_prompt_signs")
            vim.fn.sign_define("AgenticPrompt", { text = "❯", texthl = "NonText" })
            local function place_prompt_signs(bufnr)
                vim.fn.sign_unplace("AgenticPrompt", { buffer = bufnr })
                local lines = vim.api.nvim_buf_get_lines(bufnr, 0, -1, false)
                for i, line in ipairs(lines) do
                    if line == "##" then
                        vim.fn.sign_place(0, "AgenticPrompt", "AgenticPrompt", bufnr, { lnum = i })
                    end
                end
            end
            vim.api.nvim_create_autocmd("FileType", {
                pattern = "AgenticChat",
                callback = function(ev)
                    local bufnr = ev.buf
                    vim.api.nvim_create_autocmd("TextChanged", {
                        buffer = bufnr,
                        callback = function() place_prompt_signs(bufnr) end,
                    })
                    -- Jump between prompts
                    vim.keymap.set("n", "[[", function()
                        local row = vim.api.nvim_win_get_cursor(0)[1]
                        local lines = vim.api.nvim_buf_get_lines(bufnr, 0, row - 1, false)
                        for i = #lines, 1, -1 do
                            if lines[i] == "##" then
                                vim.api.nvim_win_set_cursor(0, { i, 0 })
                                return
                            end
                        end
                    end, { buffer = bufnr, desc = "Previous prompt" })
                    vim.keymap.set("n", "]]", function()
                        local row = vim.api.nvim_win_get_cursor(0)[1]
                        local lines = vim.api.nvim_buf_get_lines(bufnr, row, -1, false)
                        for i, line in ipairs(lines) do
                            if line == "##" then
                                vim.api.nvim_win_set_cursor(0, { row + i, 0 })
                                return
                            end
                        end
                    end, { buffer = bufnr, desc = "Next prompt" })
                end,
            })
            -- Cleaner chat messages: strip noise from headers and markers
            function SM._generate_welcome_header(_, session_id)
                local ts = os.date("%Y-%m-%d %H:%M")
                local short_id = session_id and session_id:sub(1, 8) or "unknown"
                return string.format("# %s · %s", ts, short_id)
            end
            local MW = require("agentic.ui.message_writer")
            local orig_write_msg = MW.write_message
            function MW:write_message(update, ...)
                if update and update.content and update.content.text then
                    local t = update.content.text
                    -- User header line → "##" (strip icon, name, timestamp)
                    t = t:gsub("\n*##[^\n]* User %- %d%d%d%d%-%d%d%-%d%d %d%d:%d%d:%d%d\n+", "\n\n##\n\n")
                    t = t:gsub("\n\n+###.* Agent %- .+$", "\n\n---")
                    t = t:gsub("^\n+", "")
                    -- Strip any finish marker line (🏁) and optional trailing "-----"
                    t = t:gsub("\n?###[^\n]*🏁[^\n]*\n?%-%-%-%-%-?", "")
                    t = t:gsub("\n?###[^\n]*🏁[^\n]*", "")
                    update.content.text = t
                end
                return orig_write_msg(self, update, ...)
            end
            -- Hide thinking chunks from chat
            local orig_write_chunk = MW.write_message_chunk
            function MW:write_message_chunk(update, ...)
                if update.sessionUpdate == "agent_thought_chunk" then return end
                return orig_write_chunk(self, update, ...)
            end
        end,
    },
}
