local ts = vim.treesitter

local M = {}

local registered = false

M.get_text = function(node)
    return ts.get_node_text(node, 0)
end

M.in_set = function()
    local node = vim.treesitter.get_node()
    if node == nil or node:type() ~= "argument_list" then return false end
    local parent = node:parent()
    if parent == nil then return false end
    local func_name = M.get_text(parent:child(0))
    return func_name == "set" or func_name == "cmd.set"
end

M.setup = function()
    if registered then return end
    registered = true

    local has_cmp, cmp = pcall(require, 'cmp')
    if not has_cmp then return end

    local rtp = vim.opt.runtimepath:get()[1]
    local fh = io.open(rtp .. "/lua/completion/pymol/pymol_settings.txt")
    local items = {}
    for line in fh:lines() do
        table.insert(items, {label=line})
    end
    fh:close()

    local types = require 'cmp.types'
    for _, item in ipairs(items) do
        item["kind"] = types.lsp.CompletionItemKind.Property
        item["insertText"] = item["label"] .. "="
    end

    -- the rest is modified from :help cmp-develop

    local source = {}

    ---Return the keyword pattern for triggering completion (optional).
    ---If this is ommited, nvim-cmp will use a default keyword pattern. See |cmp-config.completion.keyword_pattern|.
    ---@return string
    function source:get_keyword_pattern()
        return [[\k]]
    end

    -- Adding trigger chars so we get suggestions pop up automatically (otherwise it says **kwargs for my custom set function which is useless)
    function source:get_trigger_characters()
        return { ',', ' ' }
    end

    ---Invoke completion (required).
    ---@param params cmp.SourceCompletionApiParams
    ---@param callback fun(response: lsp.CompletionResponse|nil)
    function source:complete(params, callback)
        if M.in_set() then
            callback {items=items, isIncomplete=true}
            return
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
    require'cmp'.register_source('pymol_settings', source)
end

return M
