--- Lua replacement for vim-abolish's :Abolish iabbrev generator.
--- Mirrors s:expand_braces + s:create_dictionary semantics: brace expansion
--- with parallel pairing, recursive multi-brace handling, and 3-way case
--- variants (lowercase, MixedCase, UPPERCASE) plus the literal lhs.
--- See plugin/abolish.vim in tpope/vim-abolish for the original.

local M = {}

--- MixedCase: uppercase first char of camelcase form. Matches vim-abolish
--- s:mixedcase, which wraps s:camelcase (snake_case → camelCase, then
--- uppercase first char).
--- @param word string
--- @return string
function M.mixedcase(word)
    word = word:gsub("%-", "_")
    if not word:find("_") and word:find("%l") then
        return word:sub(1, 1):upper() .. word:sub(2)
    end
    local out = {}
    local upper_next = false
    for i = 1, #word do
        local c = word:sub(i, i)
        if c == "_" then
            upper_next = true
        elseif upper_next then
            out[#out + 1] = c:upper()
            upper_next = false
        else
            out[#out + 1] = c:lower()
        end
    end
    local s = table.concat(out)
    if #s > 0 then s = s:sub(1, 1):upper() .. s:sub(2) end
    return s
end

local function split_commas(s)
    local parts = {}
    for part in (s .. ","):gmatch("([^,]*),") do
        parts[#parts + 1] = part
    end
    return parts
end

--- Expand {a,b,c} braces in (lhs, rhs) into a list of concrete pairs.
--- Rules (from vim-abolish):
---   - First {...} in lhs is expanded; remaining braces handled recursively.
---   - rhs's first {...} pairs by index (mod len) with lhs targets.
---   - rhs's empty `{}` copies the lhs target verbatim.
---   - rhs without any `{...}` maps every expansion to the literal rhs.
--- @param lhs string
--- @param rhs string
--- @return string[][] pairs list of {lhs, rhs}
function M.expand_braces(lhs, rhs)
    local kbefore, kmiddle, kafter = lhs:match("^(.-){([^}]*)}(.*)$")
    if not kbefore then
        return { { lhs, rhs } }
    end

    local vbefore, vmiddle, vafter = rhs:match("^(.-){([^}]*)}(.*)$")
    if not vbefore then
        vbefore, vmiddle, vafter = rhs, ",", ""
    end

    local targets = split_commas(kmiddle)
    local replacements = split_commas(vmiddle)
    if #replacements == 1 and replacements[1] == "" then
        replacements = targets
    end

    local out = {}
    for i, target in ipairs(targets) do
        local rep = replacements[((i - 1) % #replacements) + 1]
        local new_lhs = kbefore .. target .. kafter
        local new_rhs = vbefore .. rep .. vafter
        for _, pair in ipairs(M.expand_braces(new_lhs, new_rhs)) do
            out[#out + 1] = pair
        end
    end
    return out
end

--- Case variants for one (lhs, rhs): mixed/lower/upper plus the literal.
--- Duplicates collapse via dict semantics (e.g. all-lower input yields 3 keys).
--- @param lhs string
--- @param rhs string
--- @return table<string, string> variants
function M.case_variants(lhs, rhs)
    local out = {}
    out[M.mixedcase(lhs)] = M.mixedcase(rhs)
    out[lhs:lower()] = rhs:lower()
    out[lhs:upper()] = rhs:upper()
    out[lhs] = rhs
    return out
end

-- Registry of (expanded) lhs → {rhs, predicate} entries for the `<expr>` form
-- that every abbrev uses. Buffer-local abbrevs live under their bufnr (looked
-- up first, so they shadow globals — matching vim); global abbrevs live in
-- `_dispatch_global`. Per-buffer entries go stale when a buf is wiped but are
-- only reachable through that buffer's iabbrevs, cleared at the same time, so
-- the leak is harmless.
M._dispatch = {}
M._dispatch_global = {}

--- Whether expanding `lhs` here would be at a genuine word start, reading the
--- actual buffer rather than trusting vim's full-id trigger. vim expands a
--- keyword (full-id) abbrev when the char before it is a non-keyword char *or
--- where the current insertion started* (`:help abbreviations`). That second
--- clause is a loophole: anything that resets the insert anchor mid-word — a
--- cursor move, or blink.cmp `auto_insert` rewriting the buffer as you type —
--- makes `id`→`I'd` fire inside "undid", `im`→`I'm` inside "nvim". Re-checking
--- the real preceding char closes it. Only constrains keyword-initial (full-id)
--- lhs; end-id/non-id abbrevs keep vim's own front rule.
--- @param lhs string
--- @return boolean
local function at_word_boundary(lhs)
    if vim.fn.match(lhs:sub(1, 1), "\\k") < 0 then return true end  -- not full-id
    local col = vim.api.nvim_win_get_cursor(0)[2]
    local prev = vim.api.nvim_get_current_line():sub(1, col - #lhs):sub(-1)
    return prev == "" or vim.fn.match(prev, "\\k") < 0
end

--- Called from the `<expr>` iabbrevs that `iabbrev` registers. Returns the rhs
--- when it would expand here (at a word boundary, and the predicate passes if
--- any), else the lhs (no expansion). Buffer-local entries shadow globals.
--- @param lhs string
--- @return string
function M._dispatch_lookup(lhs)
    local entry = (M._dispatch[vim.api.nvim_get_current_buf()] or {})[lhs]
        or M._dispatch_global[lhs]
    if not entry then return lhs end
    if entry.predicate and not entry.predicate() then return lhs end
    if not at_word_boundary(lhs) then return lhs end
    return entry.rhs
end

--- Check that `lhs` matches one of vim's three accepted abbreviation shapes
--- (`:help abbreviations`): full-id (all keyword chars), end-id (last char
--- keyword, others non-keyword), or non-id (last char non-keyword, others
--- any non-whitespace). Returns `true` on success, or `false, msg`.
--- @param lhs string
--- @return boolean, string?
function M.validate_lhs(lhs)
    if lhs == "" then return false, "empty lhs" end
    if vim.fn.match(lhs, "[ \t]") >= 0 then
        return false, ("lhs %q contains space or tab"):format(lhs)
    end
    local chars = vim.fn.split(lhs, "\\zs")
    local kwd = {}
    for i, c in ipairs(chars) do
        kwd[i] = vim.fn.match(c, "\\k") == 0
    end
    local n = #chars
    if not kwd[n] then return true end  -- non-id
    local has_kwd, has_non_kwd = false, false
    for i = 1, n - 1 do
        if kwd[i] then has_kwd = true else has_non_kwd = true end
    end
    if not (has_kwd and has_non_kwd) then return true end  -- full-id or end-id
    return false, ("lhs %q is not full-id, end-id, or non-id "
        .. "(see :help abbreviations) — keyword chars before a final keyword "
        .. "char must be all-keyword or all-non-keyword"):format(lhs)
end

--- Define insert-mode abbreviation(s). With `cases` (default true), also
--- expands {a,b} braces and emits lower/Title/UPPER variants — i.e. the
--- :Abolish behaviour. With `cases = false`, emits a single literal iabbrev.
--- With `buf_local = true`, registers as `<buffer>` (current buffer only).
--- With `predicate`, expansion is additionally gated on it (called at expansion
--- time, no args; expands only when truthy).
---
--- Every variant is registered as an `<expr>` iabbrev dispatching through
--- `_dispatch_lookup`, so the word-boundary guard (and any predicate) applies
--- uniformly. A plain `iabbrev` would let a full-id abbrev expand mid-word
--- whenever the insert anchor resets (see `at_word_boundary`).
---
--- Each generated lhs variant must satisfy vim's abbreviation rules
--- (`:help abbreviations`):
---   - full-id: all keyword chars (e.g. "foo", "bar123")
---   - end-id:  last char keyword, all others non-keyword (e.g. "'s", "#i")
---   - non-id:  last char non-keyword, others any non-whitespace (e.g. "def#")
--- Mixed keyword/non-keyword chars before a final keyword char (e.g. "it;s")
--- match none of the three and raise an error at registration.
--- @param lhs string
--- @param rhs string
--- @param cases? boolean default true
--- @param buf_local? boolean default false
--- @param predicate? fun(): boolean called at expansion time
function M.iabbrev(lhs, rhs, cases, buf_local, predicate)
    if cases == nil then cases = true end
    local buf_tag = buf_local and "<buffer> " or ""
    local bucket = M._dispatch_global
    if buf_local then
        local buf = vim.api.nvim_get_current_buf()
        M._dispatch[buf] = M._dispatch[buf] or {}
        bucket = M._dispatch[buf]
    end
    for _, pair in ipairs(M.expand_braces(lhs, rhs)) do
        local l, r = pair[1], pair[2]
        local entries = cases and M.case_variants(l, r) or { [l] = r }
        for k, v in pairs(entries) do
            local ok, err = M.validate_lhs(k)
            if not ok then
                error(("iabbrev(%q, %q): %s"):format(lhs, rhs, err), 2)
            end
            bucket[k] = { rhs = v, predicate = predicate }
            vim.cmd(string.format(
                "iabbrev <expr> %s%s v:lua.require'utils.iabbrev'._dispatch_lookup(%q)",
                buf_tag, k, k))
        end
    end
end

return M
