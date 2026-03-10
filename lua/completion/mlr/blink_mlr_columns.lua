local Kind = vim.lsp.protocol.CompletionItemKind

--- Cache: absolute_path -> { mtime = number, columns = string[] }
local cache = {}

--- Get text of a named child node.
local function node_text(node)
    return vim.treesitter.get_node_text(node, 0)
end

--- Get the command_name text for a command node.
local function command_name(cmd_node)
    for child in cmd_node:iter_children() do
        if child:type() == "command_name" then
            return node_text(child)
        end
    end
end

--- Check if a command node ends with `+` (mlr verb chain continuation).
local function ends_with_plus(cmd_node)
    local last
    for child in cmd_node:iter_children() do
        if child:named() then last = child end
    end
    return last and node_text(last) == "+"
end

--- Walk backwards from the current command node through `+`-chained siblings
--- to find the root `mlr` command. Returns the mlr command node and an ordered
--- list of all command nodes in the chain, or nil if not in an mlr chain.
local function find_mlr_chain(start_node)
    -- Walk up to the nearest `command` node
    local node = start_node
    while node and node:type() ~= "command" do
        node = node:parent()
    end
    -- Cursor in whitespace (e.g. after `-f `, or between `|` and next word)
    -- lands on a container node (program, pipeline, etc.) instead of command.
    -- Find the last command node whose end <= cursor position.
    if not node or node:type() ~= "command" then
        local row, col = unpack(vim.api.nvim_win_get_cursor(0))
        row = row - 1
        -- Walk up to a container that has command children
        local container = start_node
        while container and container:type() ~= "program" do
            container = container:parent()
        end
        if not container then return end
        node = nil
        -- Recursive search for the closest preceding command node
        local function scan(parent)
            for child in parent:iter_children() do
                if child:type() == "command" then
                    local _, _, er, ec = child:range()
                    if er < row or (er == row and ec <= col) then
                        node = child
                    end
                else
                    scan(child)
                end
            end
        end
        scan(container)
        if not node then return end
    end

    -- Collect chain going backwards
    local chain = { node }
    local cur = node
    while true do
        local prev = cur:prev_named_sibling()
        if not prev or prev:type() ~= "command" then break end
        if not ends_with_plus(prev) then break end
        table.insert(chain, 1, prev)
        cur = prev
    end

    -- First node in chain must be `mlr`
    if command_name(chain[1]) ~= "mlr" then return end

    -- Also collect forward continuations (cursor might be mid-chain)
    cur = node
    while ends_with_plus(cur) do
        local nxt = cur:next_named_sibling()
        if not nxt or nxt:type() ~= "command" then break end
        chain[#chain + 1] = nxt
        cur = nxt
    end

    return chain
end

--- Collect all word children of a command node as a flat list of strings.
local function collect_words(cmd_node)
    local words = {}
    for child in cmd_node:iter_children() do
        local t = child:type()
        if t == "word" or t == "raw_string" or t == "string" then
            local text = node_text(child)
            -- Strip quotes from strings
            text = text:gsub("^['\"]", ""):gsub("['\"]$", "")
            words[#words + 1] = text
        elseif t == "command_name" then
            words[#words + 1] = node_text(child)
        end
    end
    return words
end

--- Extract file paths from a chain of mlr command nodes.
local function collect_paths(chain)
    local paths = {}
    local all_words = {}

    for i, cmd_node in ipairs(chain) do
        local words = collect_words(cmd_node)
        all_words[i] = words
    end

    for i, words in ipairs(all_words) do
        local j = 1
        while j <= #words do
            local w = words[j]
            -- --from <path>
            if w == "--from" and words[j + 1] then
                paths[#paths + 1] = words[j + 1]
                j = j + 2
            -- join -f <path> (not the global -f format flag)
            elseif w == "-f" and j > 1 then
                -- Check if a preceding word in this command is "join"
                local after_join = false
                for k = 1, j - 1 do
                    if words[k] == "join" then after_join = true; break end
                end
                if after_join and words[j + 1] then
                    paths[#paths + 1] = words[j + 1]
                end
                j = j + 2
            else
                j = j + 1
            end
        end

        -- Trailing filenames from the last command in the chain
        if i == #all_words then
            for k = #words, 1, -1 do
                local w = words[k]
                if w == "+" then
                    -- skip
                elseif w:match("^%-") then
                    break
                elseif w:match("[%./]") then
                    paths[#paths + 1] = w
                else
                    break
                end
            end
        end
    end

    return paths
end

--- Resolve a path relative to the buffer's directory.
local function resolve_path(path)
    if path:sub(1, 1) == "/" or path:sub(1, 1) == "~" then
        return vim.fn.expand(path)
    end
    local buf_dir = vim.fn.fnamemodify(vim.api.nvim_buf_get_name(0), ":h")
    return buf_dir .. "/" .. path
end

--- Read column headers from a file, using cache.
local function read_headers(abs_path)
    local stat = vim.uv.fs_stat(abs_path)
    if not stat then return {} end

    local entry = cache[abs_path]
    if entry and entry.mtime == stat.mtime.sec then
        return entry.columns
    end

    -- Detect format from extension
    local is_csv = abs_path:match("%.csv") ~= nil
    local fmt = is_csv and "-c" or "-t"

    -- head -n 1 outputs header + 1 data row; first line is the header
    local result = vim.system(
        { "mlr", fmt, "head", "-n", "1", abs_path },
        { text = true }
    ):wait()

    if result.code ~= 0 or not result.stdout or result.stdout == "" then
        return {}
    end

    local header = result.stdout:match("^([^\n]+)")
    if not header then return {} end

    local sep = is_csv and "," or "\t"
    local columns = vim.split(header, sep, { plain = true })

    cache[abs_path] = { mtime = stat.mtime.sec, columns = columns }
    return columns
end

local M = {}

function M.new()
    return setmetatable({}, { __index = M })
end

function M:get_completions(_, callback)
    local empty = { is_incomplete_forward = true, is_incomplete_backward = true, items = {} }

    local node = vim.treesitter.get_node()
    if not node then
        callback(empty)
        return function() end
    end

    local chain = find_mlr_chain(node)
    if not chain then
        callback(empty)
        return function() end
    end

    local paths = collect_paths(chain)
    if #paths == 0 then
        callback(empty)
        return function() end
    end

    -- Deduplicate columns across files, tracking source
    local seen = {}
    local items = {}
    for _, path in ipairs(paths) do
        local abs = resolve_path(path)
        local columns = read_headers(abs)
        local basename = vim.fn.fnamemodify(path, ":t")
        for _, col in ipairs(columns) do
            if not seen[col] then
                seen[col] = true
                items[#items + 1] = {
                    label = col,
                    kind = Kind.Field,
                    labelDetails = { description = basename },
                    source = "mlr",
                }
            end
        end
    end

    callback({
        is_incomplete_forward = false,
        is_incomplete_backward = false,
        items = items,
    })
    return function() end
end

return M
