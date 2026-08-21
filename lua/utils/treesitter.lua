local M = {}

---Climb from `start` (default: node under the cursor) to the nearest ancestor
---satisfying `match`. `match` is either a node-type string, or a predicate
---returning a truthy value; the returned value is that predicate's result (a
---bare type string yields the matching node), so this doubles as an extractor.
---@param match string|fun(node: TSNode): any
---@param start TSNode|nil
---@return any
function M.ancestor(match, start)
    local pred = match
    if type(match) == "string" then
        pred = function(n) return n:type() == match and n or nil end
    end
    local node = start or vim.treesitter.get_node()
    while node do
        local result = pred(node)
        if result then return result end
        node = node:parent()
    end
end

---Advance `n` bytes forward from a start position through `text`, tracking
---(row, col, byte). Newline-aware: `\n` increments the row and resets the
---byte-column to 0. All coordinates are 0-indexed, matching `TSNode:range(true)`.
---@param text string
---@param n integer bytes to advance
---@param row integer
---@param col integer byte-column
---@param byte integer byte offset
---@return integer row
---@return integer col
---@return integer byte
local function advance_bytes(text, n, row, col, byte)
    for i = 1, n do
        if text:byte(i) == 0x0a then
            row = row + 1
            col = 0
        else
            col = col + 1
        end
        byte = byte + 1
    end
    return row, col, byte
end

---Query directive `(#trim! @cap PREFIX_BYTES SUFFIX_BYTES)`.
---Trims PREFIX_BYTES from the start and SUFFIX_BYTES from the end of @cap, then
---sets metadata.range to a (row, col, byte) triple consistent across all three
---coordinates. Both ends are computed by a newline-aware forward walk, so the
---range stays valid when the node spans lines — or when error recovery of an
---unterminated string ends the node at column 0, where `#offset!`'s naive
---(row+drow, col+dcol) arithmetic would emit a negative column and crash
---`set_included_ranges` ("Range value out of bounds"). Both walks are clamped to
---the node's own bytes, so trimming more than the node holds yields a
---zero-length range at the trim start — an empty injected region — rather than a
---range reaching past the node.
---@param match table<integer, TSNode[]>
---@param _ integer pattern index (unused)
---@param source integer|string buffer or string
---@param pred any[]
---@param metadata vim.treesitter.query.TSMetadata
function M.trim_directive(match, _, source, pred, metadata)
    local capture_id = pred[2]
    local prefix_bytes = tonumber(pred[3]) or 0
    local suffix_bytes = tonumber(pred[4]) or 0
    local nodes = match[capture_id]
    if not nodes or #nodes == 0 then return end
    local node = nodes[1]
    local sr, sc, sb = node:range(true)
    local text = vim.treesitter.get_node_text(node, source)
    local prefix = math.min(prefix_bytes, #text)
    local keep = math.max(#text - suffix_bytes, prefix)
    local new_sr, new_sc, new_sb = advance_bytes(text, prefix, sr, sc, sb)
    local new_er, new_ec, new_eb = advance_bytes(text, keep, sr, sc, sb)
    if not metadata[capture_id] then metadata[capture_id] = {} end
    metadata[capture_id].range = { new_sr, new_sc, new_sb, new_er, new_ec, new_eb }
end

---Query directive `(#inject-by-ext! @dest)` — set `injection.language` from the
---file extension of @dest (a heredoc redirect destination). Resolves extension
---→ filetype (`vim.filetype.match`) → parser language
---(`vim.treesitter.language.get_lang`, which honours the sh/bash → zsh
---registrations). Surrounding quotes on the destination are stripped first.
---Unknown extensions leave `injection.language` unset, so the capture is
---ignored and the base heredoc-tag injection (`<<LUA ... LUA`) still applies.
---@param match table<integer, TSNode[]>
---@param _ integer pattern index (unused)
---@param source integer|string buffer or string
---@param pred any[]
---@param metadata vim.treesitter.query.TSMetadata
function M.inject_by_ext_directive(match, _, source, pred, metadata)
    local capture_id = pred[2]
    local nodes = match[capture_id]
    if not nodes or #nodes == 0 then return end
    local dest = vim.treesitter.get_node_text(nodes[1], source)
    dest = dest:gsub("^['\"]", ""):gsub("['\"]$", "")
    local ft = vim.filetype.match({ filename = dest })
    if not ft then return end
    metadata["injection.language"] = vim.treesitter.language.get_lang(ft) or ft
end

---Basename of `node`'s text — any leading `dir/` components stripped, so
---`.venv/bin/python` yields `python`.
---@param node TSNode
---@param source integer|string buffer or string
---@return string
local function basename(node, source)
    return (vim.treesitter.get_node_text(node, source):gsub(".*/", ""))
end

---Interpreter basename → { command-flag char, injection language }. The `char`
---is the getopt command flag that makes the interpreter read its next arg as
---code (`-c` for shells/python, `-e` for julia/node/lua/R); it gates against a
---flag that merely ends in the same letter on an unrelated command.
---@type table<string, { char: string, lang: string }>
local INTERPRETERS = {
    python  = { char = "c", lang = "python" },
    python3 = { char = "c", lang = "python" },
    sh      = { char = "c", lang = "zsh" },
    bash    = { char = "c", lang = "zsh" },
    zsh     = { char = "c", lang = "zsh" },
    julia   = { char = "e", lang = "julia" },
    node    = { char = "e", lang = "javascript" },
    lua     = { char = "e", lang = "lua" },
    R       = { char = "e", lang = "r" },
    Rscript = { char = "e", lang = "r" },
    gnuplot = { char = "e", lang = "gnuplot" },
}

---Walk backward from `node` over `-`-prefixed tokens (flags and the bare `-`
---stdin marker) to the interpreter token, and return its INTERPRETERS entry (or
---nil if no non-flag token is reached or its basename is off-table). `node` is
---the token adjacent to the code: a `-c`/`-e` flag's previous sibling, or a
---command's last argument.
---@param node TSNode|nil
---@param source integer|string buffer or string
---@return { char: string, lang: string }|nil
local function resolve_interp(node, source)
    while node and vim.treesitter.get_node_text(node, source):sub(1, 1) == "-" do
        node = node:prev_named_sibling()
    end
    if not node then return end
    return INTERPRETERS[basename(node, source)]
end

---Commands that run another command named by their first non-option argument.
---@type table<string, true>
local WRAPPERS = {
    command = true, builtin = true, exec = true, eval = true,
    nohup = true, setsid = true, time = true,
    sudo = true, doas = true, env = true, xargs = true,
    timeout = true, nice = true, ionice = true, stdbuf = true, unbuffer = true,
}

---True when `node` could name a command — false for the tokens a wrapper puts
---before the wrapped name: options, a `timeout`/`nice` numeric argument, and
---`env`-style `NAME=value` assignments.
---@param node TSNode
---@param source integer|string buffer or string
---@return boolean
local function may_be_command(node, source)
    if node:type() == "number" then return false end
    local text = vim.treesitter.get_node_text(node, source)
    return text:sub(1, 1) ~= "-" and text:match("^[%w_]+=") == nil
end

---Resolve the command a token actually invokes, seeing through wrapper
---prefixes: `command grep`, `sudo grep`, `env A=1 grep` and `timeout 5 grep`
---all resolve to the `grep` node. This is the forward counterpart of
---`resolve_interp`, which walks *backward* to the interpreter.
---
---tree-sitter-zsh keeps wrapped commands flat — the wrapper is the
---`command_name` and the wrapped command an ordinary `word` argument — so the
---walk advances over named siblings to the first that `may_be_command`, and
---repeats while that token is itself a wrapper.
---
---May return a token that is no command at all: the walk holds no per-wrapper
---knowledge of which options take a value, so `sudo -u me grep` yields `me` and
---`timeout 5s grep` yields `5s`. Such a token matches no command name, so
---name-based tests fail closed rather than matching the wrong command. Preserve
---that property when extending WRAPPERS.
---@param node TSNode a `command_name` or argument token
---@param source integer|string buffer or string
---@return TSNode
local function effective_command(node, source)
    while WRAPPERS[basename(node, source)] do
        local arg = node:next_named_sibling()
        while arg and not may_be_command(arg, source) do
            arg = arg:next_named_sibling()
        end
        if not arg then return node end
        node = arg
    end
    return node
end

---Query directive `(#inject-interp! @flag)` — resolve an interpreter injection
---from the command flag adjacent to the injected string, in O(1) concurrent
---partial matches. A floating `@_interp` capture is O(command length) and
---silently starves the injection out past tree-sitter's match_limit (256) on
---long commands.
---
---1. Flag-shape gate: @flag must be a single-dash cluster (`-…`, not `--long`)
---   whose last char is a command char (`c`/`e`). This subsumes plain `-c`/`-e`
---   and clusters where the command flag is last (`-lc`, `-uc`).
---2. Backward walk: from @flag, skip `-`-prefixed flags and take the first
---   non-flag token (an argument word or the command_name) as the interpreter.
---   Handles every wrapper/prefix uniformly (`uv run python -c`,
---   `timeout 180 python3 -u -c`, `env -i bash -c`).
---3. Table lookup: set `injection.language` iff the interpreter basename is in
---   INTERPRETERS and its `char` matches the flag's command char. Otherwise
---   leave it unset — the capture is ignored (same contract as
---   `#inject-by-ext!`), letting `nvim -c` / `grep -e` fall through.
---@param match table<integer, TSNode[]>
---@param _ integer pattern index (unused)
---@param source integer|string buffer or string
---@param pred any[]
---@param metadata vim.treesitter.query.TSMetadata
function M.inject_interp_directive(match, _, source, pred, metadata)
    local nodes = match[pred[2]]
    if not nodes or #nodes == 0 then return end
    local flag = vim.treesitter.get_node_text(nodes[1], source)
    if flag:sub(1, 1) ~= "-" or flag:sub(2, 2) == "-" then return end
    local cmd_char = flag:sub(-1)
    if cmd_char ~= "c" and cmd_char ~= "e" then return end
    local interp = resolve_interp(nodes[1]:prev_named_sibling(), source)
    if interp and interp.char == cmd_char then
        metadata["injection.language"] = interp.lang
    end
end

---Query directive `(#inject-interp-cmd! @cmd)` — like `#inject-interp!` but for
---an interpreter that reads its code from stdin via a heredoc, with no flag
---(`gnuplot <<GP … GP`, `python <<EOF … EOF`, `uv run python - <<EOF … EOF`).
---@cmd is the whole interpreter (command) node. The interpreter may be wrapped
---(`uv run python`) and/or trailed by a bare `-` stdin marker and flags
---(`python -u -`), so walk backward from the last child skipping `-`-prefixed
---tokens (flags and the bare `-`) to the interpreter token — the same backward
---walk as `#inject-interp!`, but anchored on the command's last argument rather
---than a `-c`/`-e` flag. Resolves that token's basename in the INTERPRETERS
---table and sets `injection.language` to its `lang` (the `char` command flag is
---irrelevant here). An off-table interpreter leaves the language unset, so the
---capture is ignored (same contract as `#inject-by-ext!`).
---@param match table<integer, TSNode[]>
---@param _ integer pattern index (unused)
---@param source integer|string buffer or string
---@param pred any[]
---@param metadata vim.treesitter.query.TSMetadata
function M.inject_interp_cmd_directive(match, _, source, pred, metadata)
    local nodes = match[pred[2]]
    if not nodes or #nodes == 0 then return end
    local cmd = nodes[1]
    local interp = resolve_interp(cmd:named_child(cmd:named_child_count() - 1), source)
    if interp then metadata["injection.language"] = interp.lang end
end

---Query predicate `(#command-is? @cap "name" ...)` — true when the command @cap
---invokes is one of the listed names. @cap is a `(command_name)` node. Wrapper
---prefixes are resolved first (`effective_command`, so `sudo grep` counts as
---`grep`), then the basename is compared, so a path-prefixed command such as
---`.venv/bin/python` or `/usr/bin/python3` matches the bare name.
---@param match table<integer, TSNode[]>
---@param _ integer pattern index (unused)
---@param source integer|string buffer or string
---@param pred any[]
---@return boolean
function M.command_is(match, _, source, pred)
    for _, node in ipairs(match[pred[2]] or {}) do
        local name = basename(effective_command(node, source), source)
        for i = 3, #pred do
            if name == pred[i] then return true end
        end
    end
    return false
end

---Query predicate `(#arg-after? @cap @cmd N)` — true when @cap is at least the
---Nth argument (1-indexed) of the command @cmd invokes. @cmd is a
---`(command_name)` node; wrapper prefixes are resolved first
---(`effective_command`), so the wrapper's own tokens do not fill argument slots
---and `sudo sqlite3 db 'sql'` counts `db` as argument 1 exactly as the
---unwrapped form does. Options count as arguments — the point is to require
---*some* preceding argument, not to model any command's option grammar.
---@param match table<integer, TSNode[]>
---@param _ integer pattern index (unused)
---@param source integer|string buffer or string
---@param pred any[]
---@return boolean
function M.arg_after(match, _, source, pred)
    local caps, cmds = match[pred[2]], match[pred[3]]
    if not caps or not caps[1] or not cmds or not cmds[1] then return false end
    local arg = effective_command(cmds[1], source):next_named_sibling()
    for _ = 2, tonumber(pred[4]) or 1 do
        if not arg then return false end
        arg = arg:next_named_sibling()
    end
    if not arg then return false end
    local _, _, cap_byte = caps[1]:start()
    local _, _, arg_byte = arg:start()
    return cap_byte >= arg_byte
end

---Query directive `(#head! @cap N)` — narrow @cap to its first N bytes.
---Assumes the first N bytes do not span a newline.
---@param match table<integer, TSNode[]>
---@param _ integer pattern index (unused)
---@param _source integer|string buffer or string (unused)
---@param pred any[]
---@param metadata vim.treesitter.query.TSMetadata
function M.head_directive(match, _, _source, pred, metadata)
    local capture_id = pred[2]
    local n = tonumber(pred[3]) or 1
    local nodes = match[capture_id]
    if not nodes or #nodes == 0 then return end
    local sr, sc, sb = nodes[1]:range(true)
    if not metadata[capture_id] then metadata[capture_id] = {} end
    metadata[capture_id].range = { sr, sc, sb, sr, sc + n, sb + n }
end

---Query directive `(#tail! @cap N)` — narrow @cap to its last N bytes.
---Assumes the last N bytes do not span a newline.
---@param match table<integer, TSNode[]>
---@param _ integer pattern index (unused)
---@param _source integer|string buffer or string (unused)
---@param pred any[]
---@param metadata vim.treesitter.query.TSMetadata
function M.tail_directive(match, _, _source, pred, metadata)
    local capture_id = pred[2]
    local n = tonumber(pred[3]) or 1
    local nodes = match[capture_id]
    if not nodes or #nodes == 0 then return end
    local _, _, _, er, ec, eb = nodes[1]:range(true)
    if not metadata[capture_id] then metadata[capture_id] = {} end
    metadata[capture_id].range = { er, ec - n, eb - n, er, ec, eb }
end

return M
