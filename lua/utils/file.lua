-- Bounded reads of file contents, gzip included: a read of a 10 GB file costs
-- no more than a read of a 10 kB one, so nothing here reads a whole file.

local M = {}

--- Bytes pulled from one file — after decompression — before the read is cut
--- short. Weighed once per chunk, so a read overshoots it by up to CHUNK.
M.MAX_BYTES = 4 * 2 ^ 20

local CHUNK = 2 ^ 16

-- gzip only has to keep up with MAX_BYTES; a wait this long means it is stuck.
local TIMEOUT_MS = 5000

---@alias utils.file.Sink fun(chunk: string): boolean true = enough, stop

--- Feed `path` to `sink` in chunks until it asks to stop or the file ends.
---@param path string
---@param sink utils.file.Sink
---@return boolean eof the file was read to its end
local function stream_plain(path, sink)
    local fd = assert(vim.uv.fs_open(path, "r", 438))
    local offset = 0
    while true do
        local chunk, err = vim.uv.fs_read(fd, CHUNK, offset)
        if type(chunk) ~= "string" then
            vim.uv.fs_close(fd)
            error(tostring(err))
        end
        if chunk == "" then break end
        offset = offset + #chunk
        if sink(chunk) then
            vim.uv.fs_close(fd)
            return false
        end
    end
    vim.uv.fs_close(fd)
    return true
end

--- Feed the decompressed contents of a gzipped `path` to `sink` in chunks.
--- gzip is killed the moment `sink` has had enough, so no more of the archive
--- is inflated than was asked for.
---@param path string
---@param sink utils.file.Sink
---@return boolean eof the file was read to its end
local function stream_gzip(path, sink)
    local stopped, pipe_err = false, nil
    local proc -- used in the callback, which only runs under the wait below
    proc = vim.system({ "gzip", "-cd", "--", path }, {
        stdout = function(err, data)
            if stopped then return end
            -- Throwing here would leave the pipe open and the wait below
            -- hanging until its timeout, so carry the error out instead. A nil
            -- `data` is the end of the stream, which gzip reports by exiting.
            if err then
                stopped, pipe_err = true, err
            elseif data and sink(data) then
                stopped = true
            end
            if stopped then proc:kill("sigkill") end
        end,
    })
    local res = proc:wait(TIMEOUT_MS)
    local function fail(reason)
        error(("gzip -cd %s: %s"):format(path, reason))
    end
    if pipe_err then fail(pipe_err) end
    if stopped then return false end
    -- wait() hands back nothing when even its own SIGKILL does not land.
    if not res or res.code == 124 then fail(("no output for %ds"):format(TIMEOUT_MS / 1000)) end
    if res.code ~= 0 then fail(vim.trim(res.stderr or "")) end
    return true
end

--- Read a window of lines out of a file, decompressing a `.gz` path on the fly.
--- At most MAX_BYTES are read to reach the window, so the result is empty both
--- when `first` is past the end of the file and when it is past what the cap
--- reached. A line the cap cut in half is dropped.
---@param path string
---@param first integer 1-based line to start at
---@param count integer maximum number of lines returned
---@return string[] lines
function M.read_lines(path, first, count)
    local last = first + count - 1
    local chunks = {}
    local size, newlines = 0, 0
    local function sink(chunk)
        chunks[#chunks + 1] = chunk
        size = size + #chunk
        local _, found = chunk:gsub("\n", "")
        newlines = newlines + found
        return size >= M.MAX_BYTES or newlines >= last
    end

    local stream = vim.endswith(path, ".gz") and stream_gzip or stream_plain
    local eof = stream(path, sink)

    local lines = vim.split(table.concat(chunks), "\n", { plain = true })
    -- What follows the final newline is either nothing, or a line left
    -- incomplete by the stop — a full last line only exists at a real EOF.
    if not eof or lines[#lines] == "" then table.remove(lines) end
    return vim.list_slice(lines, first, last)
end

return M
