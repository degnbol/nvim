local util = require "utils/init"

-- Theme switching commands and startup colorscheme.
-- The generated colorscheme (~/nvim/colors/generated.lua) is a hardlink
-- to ~/dotfiles/colors/generated/<theme>/colorscheme.lua, managed by switch.sh.

---Run a theme switch script, then load the colorscheme it generated. On a
---non-zero exit, notify instead and keep the current colorscheme.
---@param script string path, `~` allowed
local function switch(script)
    local cmd = { vim.fs.normalize(script) }
    local obj = vim.system(cmd, { text = true }):wait()
    if obj.code ~= 0 then
        util.notify_failure(cmd, obj)
        return
    end
    vim.cmd.colorscheme("generated")
end

vim.api.nvim_create_user_command("Dark", function() switch("~/dotfiles/colors/dark.sh") end, {})
vim.api.nvim_create_user_command("Light", function() switch("~/dotfiles/colors/light.sh") end, {})

-- Fall back to a bundled scheme when the generated one hasn't been linked in.
vim.api.nvim_create_autocmd("VimEnter", {
    callback = function()
        vim.schedule(function()
            local ok, err = pcall(vim.cmd.colorscheme, "generated")
            if ok then return end
            if #vim.api.nvim_get_runtime_file("colors/generated.*", false) > 0 then
                vim.notify(err, vim.log.levels.ERROR)
            end
            vim.cmd.colorscheme("dawnfox")
        end)
    end,
})
