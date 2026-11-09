local M = {}

--- Send `method` to every running client that supports it, except those named
--- in `opts.exclude`, and call back once all have answered or failed, or after
--- `opts.timeout_ms` with the replies so far. With no eligible client, calls
--- back at once with `{}`. Errors and timeouts are reported with vim.notify at
--- WARN, naming the client, and those clients get no entry in the replies.
--- @param method string
--- @param params table
--- @param opts { exclude: string[]|nil, timeout_ms: integer }
--- @param on_done fun(replies: { client: vim.lsp.Client, result: any }[]) ordered by client id
function M.request_all(method, params, opts, on_done)
    local excluded = opts.exclude or {}
    local clients = vim.tbl_filter(function(client)
        return not vim.list_contains(excluded, client.name)
    end, vim.lsp.get_clients { method = method })
    table.sort(clients, function(a, b) return a.id < b.id end)

    --- @type table<integer, { client: vim.lsp.Client, result: any }> by index in `clients`
    local replies = {}
    local answered = {} --- @type table<integer, true>
    local n_pending = #clients
    local done = false
    local timer = assert(vim.uv.new_timer())

    local function finish()
        if done then return end
        done = true
        timer:stop()
        timer:close()
        local ordered = {}
        for i = 1, #clients do
            ordered[#ordered + 1] = replies[i]
        end
        on_done(ordered)
    end

    local function settle(i)
        answered[i] = true
        n_pending = n_pending - 1
        if n_pending == 0 then finish() end
    end

    if n_pending == 0 then return finish() end
    timer:start(opts.timeout_ms, 0, vim.schedule_wrap(function()
        if done then return end
        for i, client in ipairs(clients) do
            if not answered[i] then
                vim.notify(("%s (%s): timed out"):format(method, client.name), vim.log.levels.WARN)
            end
        end
        finish()
    end))
    for i, client in ipairs(clients) do
        local sent = client:request(method, params, function(err, result)
            if done then return end
            if err then
                vim.notify(("%s (%s): %s"):format(method, client.name, err.message), vim.log.levels.WARN)
            else
                replies[i] = { client = client, result = result }
            end
            settle(i)
        end)
        if not sent then
            vim.notify(("%s (%s): request not sent"):format(method, client.name), vim.log.levels.WARN)
            settle(i)
        end
    end
end

--- Lines `first` to `last` (0-indexed, inclusive) of the file at `uri`, from
--- its loaded buffer if any, else from disk.
--- @param uri string
--- @param first integer
--- @param last integer
--- @return string[] fewer than asked past the last line, empty for an unreadable file
local function uri_lines(uri, first, last)
    local fname = vim.uri_to_fname(uri)
    -- Buffer names are resolved paths; bufexists matches either. With the
    -- buffer existing, bufadd only returns its number.
    if vim.fn.bufexists(fname) == 1 then
        local bufnr = vim.fn.bufadd(fname)
        if vim.api.nvim_buf_is_loaded(bufnr) then
            return vim.api.nvim_buf_get_lines(bufnr, first, last + 1, false)
        end
    end
    if not vim.uv.fs_access(fname, "R") then return {} end
    return vim.list_slice(vim.fn.readfile(fname, "", last + 1), first + 1)
end

--- Convert an LSP range in the file at `uri` between position encodings.
--- Reads the file from disk when it has no loaded buffer. Lines it cannot
--- read convert as empty.
--- @param uri string
--- @param range lsp.Range
--- @param from "utf-8"|"utf-16"|"utf-32"
--- @param to "utf-8"|"utf-16"|"utf-32"
--- @return lsp.Range
function M.convert_range(uri, range, from, to)
    -- Saves reading the file.
    if from == to then return range end
    local first = range.start.line
    local lines = uri_lines(uri, first, range["end"].line)
    --- @param pos lsp.Position
    --- @return lsp.Position
    local function convert(pos)
        local line = lines[pos.line - first + 1] or ""
        local byte = vim.str_byteindex(line, from, pos.character, false)
        return { line = pos.line, character = vim.str_utfindex(line, to, byte, false) }
    end
    return { start = convert(range.start), ["end"] = convert(range["end"]) }
end

return M
