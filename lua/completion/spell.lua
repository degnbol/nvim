#!/usr/bin/env lua
local async = require "blink.cmp.lib.async"
local Kind = vim.lsp.protocol.CompletionItemKind

local M = {}

-- Initialize a new ctags instance
function M.new(user_config)
    local self = setmetatable({}, { __index = M })
    self.is_cached = false
    self.cached_items = {}
    local dir = vim.opt.runtimepath:get()[1] .. "/spell"
    self.files = { dir .. "/en.dic", dir .. "/custom.utf8.add" }
    -- lazy load instead if manually invoked.
    -- self:_load_async()
    return self
end

function M._load_file(path, out)
    for line in io.lines(path) do
        if #line > 3 then
            table.insert(out,
                -- {
                --     label = "label",
                --     kind = lsp.CompletionItemKind.Text,
                --     insertText = "label",
                --     documentation = {
                --         kind = "markdown",
                --         value = "",
                --     },
                --     detail = "",
                --     source = "spell",
                -- }
                {
                    label=line,
                    insertText=line,
                    kind = Kind.Text,
                    source="spell",
                }
            )
        end
    end
    return out
end

function M:_lazy_load()
    if self.is_cached then return end
    for _, path in ipairs(self.files) do
        self._load_file(path, self.cached_items)
    end
    self.is_cached = true
end

function M:_load_async()
    local tasks = vim.tbl_map(function(path)
        return self:_load_file_async(path)
    end, self.files)

    async.task
        .await_all(tasks)
        :map(function(results)
            self.cached_items = vim.tbl_extend(results)
            self.is_cached = true
        end)
        :catch(function(err)
    end)
end

function M:_load_file_async(path)
    return async.task.new(function(resolve, reject)
        resolve(self._load_file(path, {}))
    end)
end

-- Get completions based on the current word prefix and cached items
function M:get_completions(context, callback)
    self:_lazy_load()
    if not self.is_cached then
        callback({ is_incomplete_forward = true, is_incomplete_backward = true, items = {} })
        return function() end
    end
	callback({
		is_incomplete_forward = false,
        is_incomplete_backward = false,
        items = self.cached_items,
	})
	return function() end
end

return M
