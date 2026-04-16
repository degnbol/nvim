local grp = vim.api.nvim_create_augroup("KittySetVar", { clear = true })

--- Set or clear a kitty user variable on the current window.
--- Skipped in headless mode: no kitty terminal is attached, and the OSC
--- bytes would contaminate captured stdout (e.g. `nvim --headless -c 'echo
--- $VIMRUNTIME'` used by Makefiles and LuaLS).
--- @param name string Variable name
--- @param value? string Value to set, or nil/omit to clear
local function set_user_var(name, value)
    if #vim.api.nvim_list_uis() == 0 then
        return
    end
    if value then
        local b64 = vim.base64.encode(value)
        io.stdout:write(("\x1b]1337;SetUserVar=%s=%s\007"):format(name, b64))
    else
        io.stdout:write(("\x1b]1337;SetUserVar=%s\007"):format(name))
    end
end

vim.api.nvim_create_autocmd({ "VimEnter", "VimResume" }, {
    group = grp,
    callback = function() set_user_var("nvim", "1") end,
})

-- Track agentic.nvim session ID for kitty session restore.
-- Fired by session_manager.lua on new session and session load (/new, /clear).
vim.api.nvim_create_autocmd("User", {
    group = grp,
    pattern = "AgenticSessionChanged",
    callback = function(ev)
        local id = ev.data and ev.data.session_id
        set_user_var("agentic_session", id)
    end,
})

-- Clear.
vim.api.nvim_create_autocmd({ "VimLeave", "VimSuspend" }, {
    group = grp,
    callback = function()
        set_user_var("nvim")
        set_user_var("agentic_session")
    end,
})

