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

--- Define insert-mode abbreviation(s). With `cases` (default true), also
--- expands {a,b} braces and emits lower/Title/UPPER variants — i.e. the
--- :Abolish behaviour. With `cases = false`, emits a single literal iabbrev.
--- @param lhs string
--- @param rhs string
--- @param cases? boolean default true
function M.iabbrev(lhs, rhs, cases)
    if cases == nil then cases = true end
    for _, pair in ipairs(M.expand_braces(lhs, rhs)) do
        local l, r = pair[1], pair[2]
        local entries = cases and M.case_variants(l, r) or { [l] = r }
        for k, v in pairs(entries) do
            vim.cmd("iabbrev " .. k .. " " .. v)
        end
    end
end

return M
