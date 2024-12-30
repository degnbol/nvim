#!/usr/bin/env lua
local async = require "blink.cmp.lib.async"
local Kind = vim.lsp.protocol.CompletionItemKind
local ts = vim.treesitter

-- treesitter context functions

local function get_text(node)
    return ts.get_node_text(node, 0)
end
local function iscall(node)
    return node:type() == "call_expression"
end
local function isattr(node)
    local text = get_text(node)
    local attrNames = {a=true, attr=true, ["…"]=true}
    return iscall(node) and attrNames[text:match("^[^(]+")]
end
local function namedArg2item(node)
    local item = get_text(node):match("([%w_]+)=")
    if item == nil then return nil end
    -- Support arbitrary xaxis2_domain etc by ignoring the 2 when looking up completion items
    item = item:gsub("axis%d+", "axis")
    return item -- written in two lines, in order to only take the first gsub returned value
end
local function func2toplevel(node)
    local funcname = get_text(node):match("^[%w_]+")
    return ({
        relayout="layout",
        ["relayout!"]="layout",
        Layout="layout",
    })[funcname] or funcname
end


local M = {}

-- Initialize a new ctags instance
function M.new(user_config)
    local self = setmetatable({}, { __index = M })
    self.is_cached = false
    self.cached_items = {}
    self.tree = {}
    self:_load_async()
    return self
end

local function blink_format(obj)
    for k, v in pairs(obj) do
        if k == "items" then
            for _, item in ipairs(v) do
                item.kind = Kind.Property
                item.insertText = item.label .. "="
                item.source = "plotly"
            end
        else
            blink_format(v)
        end
    end
end

function M:_load()
    local rtp = vim.opt.runtimepath:get()[1]
    local fh = io.open(rtp .. "/lua/completion/plotlyjs/plotlyjs.json")
    if fh == nil then
        print("Error loading plotlyjs.json")
        return
    end
    self.tree = vim.json.decode(fh:read("*a"))
    fh:close()

    -- extend with shortcut entries such as marker_color=... instead of marker=attr(color=...).
    -- Only for top level for now. May extend if valid and useful.
    for kData, vData in pairs(self.tree) do
        for kGrp, vGrp in pairs(vData) do
            if vGrp.items then
                for i, item in ipairs(vGrp.items) do
                    table.insert(vData.items, {
                        label= kGrp .. "_" .. item.label,
                        documentation = item.documentation,
                    })
                end
            end
        end
    end

    -- format all items for blink
    blink_format(self.tree)

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

function M:istoplevel(node)
    return iscall(node) and self.tree[func2toplevel(node)] ~= nil
end

function M:get_func_parents()
    local parent = ts.get_node():parent()
    if parent == nil then return {} end
    if self:istoplevel(parent) then
        return {func2toplevel(parent)}
    else
        local parents = {}
        while isattr(parent) do
            local namedArg = parent:parent()
            -- prepend
            local namedArgItem = namedArg2item(namedArg)
            if namedArgItem ~= nil then table.insert(parents, 1, namedArgItem) end
            parent = namedArg:parent()
            if parent:type() == "argument_list" then
                parent = parent:parent()
            end
        end
        if self:istoplevel(parent) then
            -- prepend
            table.insert(parents, 1, func2toplevel(parent))
            return parents
        end
    end
    return {}
end

-- Get completions based on the current word prefix and cached items
function M:get_completions(context, callback)
    -- self:_lazy_load()
    if self.is_cached then
        -- treesitter context
        local func_parents = self:get_func_parents()
        if func_parents then
            local d = self.tree
            for _, k in ipairs(func_parents) do d = d[k] end
            if d ~= nil then
                callback({
                    is_incomplete_forward = false,
                    is_incomplete_backward = false,
                    items = d["items"] or {},
                })
                return function() end
            end
        end
    end

    callback({
        is_incomplete_forward = true,
        is_incomplete_backward = true,
        items = {},
    })
    return function() end
end

return M
