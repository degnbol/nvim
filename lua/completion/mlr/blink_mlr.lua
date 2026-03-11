local Kind = vim.lsp.protocol.CompletionItemKind

--- Load mlr verb/flag data from generated JSON.
local data_path = vim.fn.fnamemodify(debug.getinfo(1, "S").source:sub(2), ":h") .. "/mlr_verbs.json"
local data = vim.json.decode(table.concat(vim.fn.readfile(data_path), "\n"))

-- Convert positional_field_verbs array to set for O(1) lookup
local positional_field_verbs = {}
for _, v in ipairs(data.positional_field_verbs) do
    positional_field_verbs[v] = true
end

-- Convert field_flags per-verb arrays to sets
local field_flags = {}
for verb, flags in pairs(data.field_flags) do
    field_flags[verb] = {}
    for _, f in ipairs(flags) do
        field_flags[verb][f] = true
    end
end

--- Cache: absolute_path -> { mtime = number, columns = string[] }
local cache = {}

local function node_text(node)
    return vim.treesitter.get_node_text(node, 0)
end

local function command_name(cmd_node)
    for child in cmd_node:iter_children() do
        if child:type() == "command_name" then
            return node_text(child)
        end
    end
end

local function ends_with_plus(cmd_node)
    local last
    for child in cmd_node:iter_children() do
        if child:named() then last = child end
    end
    return last and node_text(last) == "+"
end

--- Walk backwards from the current command node through `+`-chained siblings
--- to find the root `mlr` command. Returns an ordered list of all command
--- nodes in the chain, or nil if not in an mlr chain.
local function find_mlr_chain(start_node)
    local node = start_node
    while node and node:type() ~= "command" do
        node = node:parent()
    end
    -- Cursor in whitespace lands on a container node.
    -- Find the last command node whose end <= cursor position.
    if not node or node:type() ~= "command" then
        local row, col = unpack(vim.api.nvim_win_get_cursor(0))
        row = row - 1
        local container = start_node
        while container and container:type() ~= "program" do
            container = container:parent()
        end
        if not container then return end
        node = nil
        local function scan(parent)
            for child in parent:iter_children() do
                if child:type() == "command" then
                    local sr, sc, er, ec = child:range()
                    if (sr < row or (sr == row and sc <= col))
                        and (er > row or (er == row and ec >= col)) then
                        node = child
                    elseif er < row or (er == row and ec <= col) then
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

    local chain = { node }
    local cur = node
    while true do
        local prev = cur:prev_named_sibling()
        if not prev or prev:type() ~= "command" then break end
        if not ends_with_plus(prev) then break end
        table.insert(chain, 1, prev)
        cur = prev
    end

    if command_name(chain[1]) ~= "mlr" then return end

    cur = node
    while ends_with_plus(cur) do
        local nxt = cur:next_named_sibling()
        if not nxt or nxt:type() ~= "command" then break end
        chain[#chain + 1] = nxt
        cur = nxt
    end

    return chain
end

--- Collect all word-like children across the chain as tokens with positions.
--- Returns {{ text, sr, sc, er, ec }, ...} ordered by position.
local function collect_tokens(chain)
    local tokens = {}
    for _, cmd_node in ipairs(chain) do
        for child in cmd_node:iter_children() do
            local t = child:type()
            local text
            if t == "word" or t == "raw_string" or t == "string" then
                text = node_text(child)
                text = text:gsub("^['\"]", ""):gsub("['\"]$", "")
            elseif t == "command_name" then
                text = node_text(child)
            end
            if text then
                local sr, sc, er, ec = child:range()
                tokens[#tokens + 1] = { text = text, sr = sr, sc = sc, er = er, ec = ec }
            end
        end
    end
    return tokens
end

--- Determine cursor context within the mlr chain.
--- Returns { type = "verb" | "field" | "flag", verb = string|nil }
local function cursor_context(chain)
    local cursor_row, cursor_col = unpack(vim.api.nvim_win_get_cursor(0))
    cursor_row = cursor_row - 1

    local tokens = collect_tokens(chain)

    -- Treesitter ranges are [start, end) — ec is exclusive. In insert mode,
    -- cursor sits at ec when typing a word at EOL. Distinguish from cursor
    -- on whitespace after a token (where ec also == cursor_col) by checking
    -- the actual character: EOL or non-whitespace → user is typing; space → past it.
    local cur_line = vim.api.nvim_buf_get_lines(0, cursor_row, cursor_row + 1, false)[1] or ""
    local typing_at_boundary = cursor_col >= #cur_line
        or not cur_line:sub(cursor_col + 1, cursor_col + 1):match("%s")

    local at_cursor_idx, before_cursor_idx
    for i, tok in ipairs(tokens) do
        local starts_before = tok.sr < cursor_row or (tok.sr == cursor_row and tok.sc <= cursor_col)
        local ends_after = tok.er > cursor_row or (tok.er == cursor_row and tok.ec > cursor_col)
        local ends_at = tok.er == cursor_row and tok.ec == cursor_col

        if starts_before and ends_after then
            at_cursor_idx = i
        elseif starts_before and ends_at and typing_at_boundary then
            at_cursor_idx = i
        end
        if tok.er < cursor_row or (tok.er == cursor_row and tok.ec <= cursor_col) then
            before_cursor_idx = i
        end
    end

    -- For context, parse tokens up to (but not including) the word being typed
    local context_end = at_cursor_idx and (at_cursor_idx - 1) or (before_cursor_idx or 0)

    -- Walk tokens to determine state
    local in_mlr_flags = true
    local current_verb = nil
    local prev_text = nil
    local skip_next = false

    for i = 1, context_end do
        local text = tokens[i].text

        if skip_next then
            skip_next = false
            prev_text = text
            goto continue
        end

        if text == "+" or text == "then" then
            current_verb = nil
            in_mlr_flags = false
            prev_text = text
            goto continue
        end

        if in_mlr_flags then
            if text == "--from" or text == "--mfrom" then
                skip_next = true
            elseif not text:match("^%-") and text ~= "mlr" then
                in_mlr_flags = false
                current_verb = text
            end
        elseif not current_verb then
            current_verb = text
        end

        prev_text = text
        ::continue::
    end

    -- Verb position: no verb yet, or last token was +/then
    if not current_verb then
        return { type = "verb" }
    end
    if prev_text == "+" or prev_text == "then" then
        return { type = "verb" }
    end

    -- Field position: previous token is a field-taking flag for this verb
    local verb_ff = field_flags[current_verb]
    if prev_text and verb_ff and verb_ff[prev_text] then
        return { type = "field", verb = current_verb }
    end

    -- Positional field verbs: args after verb name are field names
    -- (unless preceded by a non-field flag that takes an argument)
    if positional_field_verbs[current_verb] then
        -- If prev_text is the verb itself, or a non-flag token, we're in field position
        if prev_text == current_verb
            or (prev_text and not prev_text:match("^%-")) then
            return { type = "field", verb = current_verb }
        end
    end

    -- Flag position: inside a verb's arguments
    return { type = "flag", verb = current_verb }
end

--- Collect all word children of a command node as a flat list of strings.
local function collect_words(cmd_node)
    local words = {}
    for child in cmd_node:iter_children() do
        local t = child:type()
        if t == "word" or t == "raw_string" or t == "string" then
            local text = node_text(child)
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
            if w == "--from" and words[j + 1] then
                paths[#paths + 1] = words[j + 1]
                j = j + 2
            elseif w == "-f" and j > 1 then
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

local function resolve_path(path)
    if path:sub(1, 1) == "/" or path:sub(1, 1) == "~" then
        return vim.fn.expand(path)
    end
    local buf_dir = vim.fn.fnamemodify(vim.api.nvim_buf_get_name(0), ":h")
    return buf_dir .. "/" .. path
end

local function read_headers(abs_path)
    local stat = vim.uv.fs_stat(abs_path)
    if not stat then return {} end

    local entry = cache[abs_path]
    if entry and entry.mtime == stat.mtime.sec then
        return entry.columns
    end

    local is_csv = abs_path:match("%.csv") ~= nil
    local fmt = is_csv and "-c" or "-t"

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

--- Extract verb segments from a command node's word children.
--- Each segment: { verb = string, args = {string...}, end_row, end_col }.
local function verb_segments(cmd_node)
    local segments = {}
    local cur_verb, cur_args, cur_end_row, cur_end_col
    local is_root = command_name(cmd_node) == "mlr"
    local past_flags = not is_root
    local skip_next = false

    for child in cmd_node:iter_children() do
        local t = child:type()
        local text
        if t == "word" or t == "raw_string" or t == "string" then
            text = node_text(child)
            text = text:gsub("^['\"]", ""):gsub("['\"]$", "")
        elseif t == "command_name" then
            text = node_text(child)
        end
        if not text then goto continue end

        local _, _, er, ec = child:range()

        if skip_next then
            skip_next = false
        elseif text == "+" then
            if cur_verb then
                segments[#segments + 1] = {
                    verb = cur_verb, args = cur_args,
                    end_row = cur_end_row, end_col = cur_end_col,
                }
            end
            cur_verb, cur_args = nil, nil
            past_flags = true
        elseif not past_flags then
            if text == "--from" then
                skip_next = true
            elseif not text:match("^%-") and text ~= "mlr" then
                past_flags = true
                cur_verb = text
                cur_args = {}
            end
        elseif not cur_verb then
            cur_verb = text
            cur_args = {}
        else
            cur_args[#cur_args + 1] = text
        end

        cur_end_row, cur_end_col = er, ec
        ::continue::
    end
    if cur_verb then
        segments[#segments + 1] = {
            verb = cur_verb, args = cur_args,
            end_row = cur_end_row, end_col = cur_end_col,
        }
    end
    return segments
end

--- Parse rename pairs from a single rename verb segment's args.
local function parse_rename_pairs(args)
    local renames = {}
    for _, arg in ipairs(args) do
        if arg:match("^%-") then break end
        local parts = vim.split(arg, ",", { plain = true })
        for k = 1, #parts - 1, 2 do
            if parts[k] ~= "" and parts[k + 1] ~= "" then
                renames[parts[k]] = parts[k + 1]
            end
        end
    end
    return renames
end

--- Analyse rename verbs in the chain relative to the cursor.
--- Returns:
---   renames: { old = new } mapping from rename verbs BEFORE cursor's verb
---   suppress: true if cursor is in a rename verb at a "to" (new name) position
local function analyse_renames(chain)
    local cursor_row, cursor_col = unpack(vim.api.nvim_win_get_cursor(0))
    cursor_row = cursor_row - 1

    local renames = {}
    local suppress = false

    for _, cmd_node in ipairs(chain) do
        local segs = verb_segments(cmd_node)
        for si, seg in ipairs(segs) do
            local seg_before_cursor = seg.end_row < cursor_row
                or (seg.end_row == cursor_row and seg.end_col <= cursor_col)
            local next_seg = segs[si + 1]
            local seg_contains_cursor = not seg_before_cursor
                and (not next_seg
                    or cursor_row < next_seg.end_row
                    or (cursor_row == next_seg.end_row and cursor_col <= next_seg.end_col))

            if seg.verb == "rename" then
                if seg_before_cursor then
                    local rename_map = parse_rename_pairs(seg.args)
                    for old, new in pairs(rename_map) do
                        renames[old] = new
                    end
                elseif seg_contains_cursor then
                    local sr, sc = cmd_node:range()
                    local lines = vim.api.nvim_buf_get_text(
                        0, sr, sc, cursor_row, cursor_col, {})
                    local text = table.concat(lines, "\n")
                    local last_pos = 1
                    while true do
                        local s = text:find("rename%s", last_pos)
                        if not s then break end
                        last_pos = s + 1
                    end
                    local after = text:sub(last_pos - 1):match("rename%s+(.*)")
                    if after then
                        local before_plus = after:match("^(.-)%s+%+") or after
                        local commas = select(2, before_plus:gsub(",", ""))
                        suppress = commas % 2 == 1
                    end
                end
            end
        end
    end

    return renames, suppress
end

--- Build column name completion items from file headers.
local function column_items(chain)
    local paths = collect_paths(chain)
    if #paths == 0 then return {} end

    local renames, suppress = analyse_renames(chain)
    if suppress then return {} end

    local seen = {}
    local items = {}
    for _, path in ipairs(paths) do
        local abs = resolve_path(path)
        local columns = read_headers(abs)
        local basename = vim.fn.fnamemodify(path, ":t")
        for _, col in ipairs(columns) do
            local display = renames[col] or col
            if not seen[display] then
                seen[display] = true
                items[#items + 1] = {
                    label = display,
                    kind = Kind.Field,
                    labelDetails = { description = basename },
                    source = "mlr",
                }
            end
        end
    end
    return items
end

--- Build verb name completion items.
local function verb_items()
    local items = {}
    for name, desc in pairs(data.verbs) do
        items[#items + 1] = {
            label = name,
            kind = Kind.Keyword,
            labelDetails = { description = desc },
            source = "mlr",
        }
    end
    return items
end

--- Build flag completion items for a specific verb.
--- Blink's keyword can't start with `-`, so typing `-f` gives keyword `f`.
--- Without textEdit, accepting `-f` replaces only `f` → `--f`.
--- Use textEdit to replace from the start of the typed prefix (including dashes).
local function flag_items(verb)
    local cursor_row, cursor_col = unpack(vim.api.nvim_win_get_cursor(0))
    local row = cursor_row - 1
    local line_text = vim.api.nvim_buf_get_lines(0, row, row + 1, false)[1] or ""
    local before = line_text:sub(1, cursor_col)
    local prefix = before:match("%-[-a-zA-Z]*$") or ""
    local edit_start = cursor_col - #prefix

    local items = {}
    local flags = data.verb_flags[verb]
    if flags then
        for _, entry in ipairs(flags) do
            local flag, desc = entry[1], entry[2]
            if flag ~= "-h" and flag ~= "--help" then
                items[#items + 1] = {
                    label = flag,
                    kind = Kind.Property,
                    labelDetails = { description = desc or "" },
                    textEdit = {
                        range = {
                            start = { line = row, character = edit_start },
                            ["end"] = { line = row, character = cursor_col },
                        },
                        newText = flag,
                    },
                    source = "mlr",
                }
            end
        end
    end
    items[#items + 1] = {
        label = "+",
        kind = Kind.Operator,
        labelDetails = { description = "Chain next verb" },
        sortText = "~~~+",
        source = "mlr",
    }
    return items
end

local M = {}

function M.new()
    return setmetatable({}, { __index = M })
end

function M:get_trigger_characters()
    return { "," }
end

function M:get_completions(_, callback)
    local empty = { is_incomplete_forward = true, is_incomplete_backward = true, items = {} }

    local node = vim.treesitter.get_node()
    if not node then
        callback(empty)
        return function() end
    end

    local nt = node:type()
    local items

    -- Inside a string: offer column names if it's a filter/put DSL expression
    if nt == "raw_string" or nt == "string" or nt == "string_content" then
        local str_node = (nt == "string_content") and node:parent() or node
        local parent = str_node:parent()
        if not parent or parent:type() ~= "command" then
            callback(empty)
            return function() end
        end
        -- Check if a preceding sibling word is "filter" or "put"
        local is_dsl = false
        for child in parent:iter_children() do
            if child == str_node then break end
            local ct = child:type()
            if ct == "word" or ct == "command_name" then
                local txt = node_text(child)
                if txt == "filter" or txt == "put" then
                    is_dsl = true
                end
            end
        end
        if not is_dsl then
            callback(empty)
            return function() end
        end
        local chain = find_mlr_chain(str_node)
        if not chain then
            callback(empty)
            return function() end
        end
        items = column_items(chain)
    else
        local chain = find_mlr_chain(node)
        if not chain then
            callback(empty)
            return function() end
        end

        local ctx = cursor_context(chain)

        if ctx.type == "verb" then
            items = verb_items()
        elseif ctx.type == "field" then
            items = column_items(chain)
        elseif ctx.type == "flag" then
            items = flag_items(ctx.verb)
        else
            callback(empty)
            return function() end
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
