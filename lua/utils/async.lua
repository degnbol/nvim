-- Helpers for async functions run by vim._async.run.
local async = require "vim._async"

local M = {}

--- Run a command and wait for it to exit.
--- @async
--- @param cmd string[]
--- @param opts vim.SystemOpts|nil see `vim.system`
--- @return vim.SystemCompleted result
function M.system(cmd, opts)
    local result = async.await(3, vim.system, cmd, opts)
    -- vim.system calls back in a fast context, where most of the API is not allowed.
    async.await(1, vim.schedule)
    return result
end

--- Share runs of an async function: a call with the key of a run that has not yet
--- finished waits for that run, and returns its result or raises its error.
--- @param fn async fun(key: string): any
--- @return async fun(key: string): any
function M.shared(fn)
    --- @type table<string, fun(ok: boolean, result: any)[]>
    local waiting = {}
    return function(key)
        local ok, result
        if waiting[key] then
            ok, result = async.await(1, function(resume) table.insert(waiting[key], resume) end)
        else
            waiting[key] = {}
            ok, result = pcall(fn, key)
            local waiters = waiting[key]
            waiting[key] = nil
            for _, resume in ipairs(waiters) do
                vim.schedule(function() resume(ok, result) end)
            end
        end
        if not ok then error(result, 0) end
        return result
    end
end

return M
