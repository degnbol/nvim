
local async = require "blink.cmp.lib.async"
local Kind = vim.lsp.protocol.CompletionItemKind
local ts = vim.treesitter

-- treesitter context functions

local function get_text(node)
    return ts.get_node_text(node, 0)
end

local function in_comment()
    -- for some reason the in_set function ignores it when we are in a comment, 
    -- and similarly vim.treesitter.get_captures_at_cursor and get_node doesn't 
    -- seem to see the comment here, even though they see it just fine when run 
    -- in normal mode. We fallback to using regex.
    local r, c = unpack(vim.api.nvim_win_get_cursor(0))
    local line = vim.api.nvim_get_current_line()
    local comment_match = line:match("^[ \t]+#")
    return comment_match ~= nil and c > #comment_match
end

local function in_set()
    local node = vim.treesitter.get_node()
    if node == nil or node:type() ~= "argument_list" then return false end
    local parent = node:parent()
    if parent == nil then return false end
    local func_name = get_text(parent:child(0))
    return func_name == "set" or func_name == "cmd.set"
end

local M = {}

-- Initialize a new ctags instance
function M.new(user_config)
    local self = setmetatable({}, { __index = M })
    self.is_cached = false
    self.cached_items = {}
    self:_load_async()
    return self
end

function M:_load()
    local rtp = vim.opt.runtimepath:get()[1]
    local fh = io.open(rtp .. "/lua/completion/pymol/pymol_settings.txt")
    if fh == nil then
        print("pymol settings file not found")
        return
    end
    for line in fh:lines() do
        table.insert(self.cached_items, {
            label=line,
            insertText=line.."=",
            kind=Kind.Property,
            source="pymol_settings",
        })
    end
    fh:close()

    -- add descriptions if available
    for _, item in ipairs(self.cached_items) do
        local filepath = rtp .. "/lua/completion/pymol/pymol_settings_descriptions/" .. item.label .. ".md"
        fh = io.open(filepath)
        if fh ~= nil then
            local content = fh:read("*a")
            fh:close()
            item["documentation"] = {
                kind = "markdown",
                value = content:gsub("<p>", ""):gsub("<\\p>", "")
            }
        end
    end
    self.is_cached = true
end

function M:_lazy_load()
    if self.is_cached then return end
    self:_load()
end

function M:_load_async()
    async.task.new(function(resolve, reject)
        resolve(self:_load())
    end)
end

-- Get completions based on the current word prefix and cached items
function M:get_completions(context, callback)
    -- self:_lazy_load()
    if not self.is_cached or in_comment() or not in_set() then
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
