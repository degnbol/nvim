-- Large file handling for both command line and session
local M = {}

M.threshold = 1024 * 1024              -- 1 MB: disable features after load
M.extreme_threshold = 50 * 1024 * 1024 -- 50 MB: defer/block events

-- Configure buffer for large file viewing
function M.configure(buf)
    vim.bo[buf].swapfile = false
    vim.bo[buf].undofile = false
    vim.bo[buf].syntax = ""
    vim.wo.foldmethod = "manual"
    vim.wo.foldenable = false
    vim.wo.spell = false
    pcall(vim.treesitter.stop, buf)
    for _, client in ipairs(vim.lsp.get_clients({ bufnr = buf })) do
        pcall(vim.lsp.buf_detach_client, buf, client.id)
    end
end

-- Check command line args for large files (call from init.lua before plugins)
function M.check_argv()
    for _, arg in ipairs(vim.fn.argv()) do
        local ok, stats = pcall(vim.uv.fs_stat, arg)
        if ok and stats and stats.size > M.extreme_threshold then
            local fullpath = vim.fn.fnamemodify(arg, ":p")
            local size = stats.size

            -- Remove from argument list and wipe buffer
            vim.cmd("silent! argdelete *")
            local bufnr = vim.fn.bufnr(fullpath)
            if bufnr ~= -1 then
                vim.cmd("silent! bwipeout! " .. bufnr)
            end

            -- Block events during startup
            vim.o.eventignore = "FileType,Syntax,BufReadPre,BufReadPost,BufEnter,BufWinEnter,BufAdd,BufNew"

            -- Open file after startup completes
            vim.api.nvim_create_autocmd("UIEnter", {
                once = true,
                callback = function()
                    vim.o.eventignore = ""
                    vim.schedule(function()
                        vim.schedule(function()
                            vim.cmd("edit " .. vim.fn.fnameescape(fullpath))
                            vim.b.largefile = true
                            M.configure(0)
                            vim.notify(
                                string.format("Large file (%.0f MB) — opened", size / 1024 / 1024),
                                vim.log.levels.WARN
                            )
                        end)
                    end)
                end,
            })
            return true
        end
    end
    return false
end

return M
