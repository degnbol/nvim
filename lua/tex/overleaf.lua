local util = require "utils/init"

local M = {}

local gitsync = vim.fn.stdpath("config") .. "/tex/overleaf/gitsync.sh"
local group = vim.api.nvim_create_augroup("overleaf", { clear = true })

---Commit a file, then pull and push, in the background. Notifies on failure.
---@param path string file to commit, inside a git repo
local function sync(path)
    local cmd = { gitsync, path }
    -- gitsync.sh runs git in its working directory.
    vim.system(cmd, { cwd = vim.fs.root(path, ".git"), text = true }, function(obj)
        util.notify_failure(cmd, obj)
    end)
end

---Is the git repo containing buf cloned from Overleaf?
---@param buf integer
---@return boolean
local function is_overleaf(buf)
    local root = vim.fs.root(buf, ".git")
    if not root then return false end
    -- Non-zero exit: no origin remote.
    local obj = vim.system({ "git", "remote", "get-url", "origin" }, { cwd = root, text = true }):wait()
    return obj.code == 0 and vim.startswith(obj.stdout, "https://git.overleaf.com/")
end

---For the current buffer, if in an Overleaf clone: set 'wrap' locally, remove
---auto-format from 'formatoptions', and sync with Overleaf in the background on
---load, save and focus. Does nothing for other buffers.
function M.setup()
    local buf = vim.api.nvim_get_current_buf()
    if not is_overleaf(buf) then return end
    -- When collaborating you probably don't want to insert linebreaks all over the place.
    vim.opt_local.wrap = true
    vim.opt_local.formatoptions:remove("a") -- no autoformat
    vim.api.nvim_clear_autocmds { group = group, buffer = buf }
    vim.api.nvim_create_autocmd({ "BufRead", "BufWritePost", "FocusGained" }, {
        group = group,
        buffer = buf,
        callback = function(args) sync(vim.api.nvim_buf_get_name(args.buf)) end,
    })
end

return M
