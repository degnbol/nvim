#!/usr/bin/env lua
local M = {}

---Check if latexmk is running.
---Checks if any process in the system matches the word latexmk.
---Maybe in future only look for latexmk running in the current project folder.
---@param callback function This function is asynchronously called with a boolean.
function M.is_running(callback)
    -- We need to grep in full name (-f), since the command is `perl /Library/TeX/texbin/latexmk ...`
    vim.system({"pgrep", "-f", "latexmk"}, {}, function (obj)
        -- returns success if process is found
        callback(obj.code == 0)
    end)
end

return M
