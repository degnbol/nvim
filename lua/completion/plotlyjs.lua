local ts = vim.treesitter

local M = {}

local registered = false

local rtp = vim.opt.runtimepath:get()[1]
local fh = io.open(rtp .. "/lua/completion/plotlyjs.json")
local tree = vim.json.decode(fh:read("*a"))
fh:close()
-- extend with shortcut entries such as marker_color=... instead of marker=attr(color=...).
-- Only for top level for now. May extend if valid and useful.
for kData, vData in pairs(tree) do
    for kGrp, vGrp in pairs(vData) do
        if vGrp["items"] then
            for i, item in ipairs(vGrp["items"]) do
                table.insert(vData["items"], {label=kGrp .. "_" .. item["label"], detail=item["detail"]})
            end
        end
    end
end

M.text = function (node)
    return ts.get_node_text(node,0)
end
M.iscall = function (node)
    return node:type() == "call_expression"
end
M.isattr = function (node)
    local text = M.text(node)
    return M.iscall(node) and (vim.startswith(text, "a(") or vim.startswith(text, "attr("))
end
local func2toplevel = {relayout="layout", ["relayout!"]="layout", Layout="layout"}
M.func2toplevel = function (node)
    local funcname = M.text(node):match("^[%w_]+")
    return func2toplevel[funcname] or funcname
end
M.istoplevel = function (node)
    return M.iscall(node) and tree[M.func2toplevel(node)] ~= nil
end

M.get_func_parents = function ()
    local parent = ts.get_node():parent()
    if parent == nil then return {} end
    if M.istoplevel(parent) then
        return {M.func2toplevel(parent)}
    else
        local parents = {}
        while M.isattr(parent) do
            local namedArg = parent:parent()
            -- prepend
            table.insert(parents, 1, M.text(namedArg):match("([%w_]+)="))
            parent = namedArg:parent()
            if parent:type() == "argument_list" then
                parent = parent:parent()
            end
        end
        if M.istoplevel(parent) then
            -- prepend
            table.insert(parents, 1, M.func2toplevel(parent))
            return parents
        end
    end
    return {}
end

M.setup = function()
    if registered then return end
    registered = true

    local has_cmp, cmp = pcall(require, 'cmp')
    if not has_cmp then return end

    -- the rest is modified from :help cmp-develop
    
    local source = {}

    ---Return the keyword pattern for triggering completion (optional).
    ---If this is ommited, nvim-cmp will use a default keyword pattern. See |cmp-config.completion.keyword_pattern|.
    ---@return string
    function source:get_keyword_pattern()
        return [[\k]]
    end

    -- main purpose of this is to separate these completion items from regular text items, 
    -- so regular text suggestions comes after due to my custom sorting in completion.lua
    local types = require('cmp.types')
    
    ---Invoke completion (required).
    ---@param params cmp.SourceCompletionApiParams
    ---@param callback fun(response: lsp.CompletionResponse|nil)
    function source:complete(params, callback)
        local func_parents = M.get_func_parents()
        if func_parents then
            local d = tree
            for _, k in ipairs(func_parents) do d = d[k] end
            if d ~= nil then
                local items = d["items"]
                if items ~= nil then
                    for _, item in ipairs(items) do
                        item["kind"] = types.lsp.CompletionItemKind.Property
                        item["insertText"] = item["label"] .. "="
                    end
                end
                callback {items=items, isIncomplete=true}
                return
            end
        end
        -- isIncomplete allows adding other completion items, such as snippets.
        callback {isIncomplete=true}
    end

    ---Resolve completion item (optional). This is called right before the completion is about to be displayed.
    ---Useful for setting the text shown in the documentation window (`completion_item.documentation`).
    ---@param completion_item lsp.CompletionItem
    ---@param callback fun(completion_item: lsp.CompletionItem|nil)
    function source:resolve(completion_item, callback)
        callback(completion_item)
    end

    ---Executed after the item was selected.
    ---@param completion_item lsp.CompletionItem
    ---@param callback fun(completion_item: lsp.CompletionItem|nil)
    function source:execute(completion_item, callback)
        callback(completion_item)
    end

    ---Register your source to nvim-cmp.
    require'cmp'.register_source('plotlyjs', source)
end

return M
