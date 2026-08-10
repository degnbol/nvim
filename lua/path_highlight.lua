-- Underline filepath tokens in markdown prose that resolve to an existing path
-- (@string.special.path is underline in the colorscheme), so the underline
-- doubles as an on-disk existence indicator. Driven by a decoration provider:
-- on_range fires per redraw for visible lines only, setting ephemeral extmarks,
-- so highlights track typing and scrolling with no change-event wiring.
local M = {}

local paths = require "utils.paths"
local ts = require "utils/treesitter"

local ns = vim.api.nvim_create_namespace("path_highlight")

-- resolved abspath -> exists (bool). Pure single-candidate resolution makes the
-- abspath a valid key (a :cd changes the key for a relative token, so a stale
-- hit is impossible). Cleared on BufWritePost / FocusGained in setup().
local exists_cache = {}

-- A maximal run of path characters is a candidate iff it starts with a real
-- anchor (/, ./, ../, ~/, $VAR/) followed by at least one more char — so a
-- lone "/" or "and/or" in prose is not a path.
local PATH_RUN = "[%w/%._%-+~=@#%%$]+"
local function is_anchored(tok)
    return tok:match("^%.?%.?/.") or tok:match("^~/.") or tok:match("^%$[%w_]+/.")
end

--- Anchored path-token candidates in a single line.
---@param line string
---@return { col: integer, token: string }[] col is a 0-indexed byte column
function M.scan(line)
    local out = {}
    local start = 1
    while true do
        local s, e = line:find(PATH_RUN, start)
        if not s then break end
        local tok = line:sub(s, e)
        if is_anchored(tok) then
            out[#out + 1] = { col = s - 1, token = tok }
        end
        start = e + 1
    end
    return out
end

--- Whether byte position (row, col) is markdown prose, excluding code (spans,
--- fenced and indented blocks) and link destinations. Two-step because markdown
--- inline content is itself an injection: step 1's language check catches only
--- a fence whose injected language has a parser (it reports e.g. "lua"); a bare
--- or unknown-language fence, an indented block, a code span, and a link
--- destination all stay "markdown"/"markdown_inline", so step 2's ancestor
--- climb (with injections un-ignored, or code_span/link_destination are
--- unreachable) rejects them by node type.
---@param bufnr integer
---@param row integer 0-indexed
---@param col integer 0-indexed byte column
---@return boolean
function M.in_prose(bufnr, row, col)
    local parser = vim.treesitter.get_parser(bufnr)
    if not parser then return false end
    local lang = parser:language_for_range({ row, col, row, col }):lang()
    if lang ~= "markdown" and lang ~= "markdown_inline" then return false end
    local node = vim.treesitter.get_node({
        bufnr = bufnr,
        pos = { row, col },
        ignore_injections = false,
    })
    if not node then return true end
    local excluded = ts.ancestor(function(n)
        local t = n:type()
        return t == "code_span"
            or t == "link_destination"
            or t == "fenced_code_block"
            or t == "indented_code_block"
    end, node)
    return not excluded
end

--- Byte-length of `token`, or of its punctuation-stripped prefix, whichever
--- resolves to an existing path; nil if neither does. The raw token is tried
--- first so a real "../.." dir wins over stripping its trailing dots, while a
--- sentence-final "./foo.md." still resolves via the stripped candidate.
---@param token string
---@param bufnr integer
---@return integer|nil
function M.existing_len(token, bufnr)
    for _, cand in ipairs({ token, (token:gsub("%.+$", "")) }) do
        if cand ~= "" then
            local abspath = paths.resolve_path(cand, bufnr)
            if abspath then
                local exists = exists_cache[abspath]
                if exists == nil then
                    exists = vim.uv.fs_stat(abspath) ~= nil
                    exists_cache[abspath] = exists
                end
                if exists then return #cand end
            end
        end
    end
end

local function on_range(_, _, bufnr, brow, _, erow, ecol)
    -- End (lnum, 0) means EOL of the preceding line is the last included point.
    local last = ecol == 0 and erow - 1 or erow
    if last < brow then return end
    local lines = vim.api.nvim_buf_get_lines(bufnr, brow, last + 1, false)
    for i, line in ipairs(lines) do
        local row = brow + i - 1
        for _, m in ipairs(M.scan(line)) do
            local len = M.existing_len(m.token, bufnr)
            if len and M.in_prose(bufnr, row, m.col) then
                vim.api.nvim_buf_set_extmark(bufnr, ns, row, m.col, {
                    end_col = m.col + len,
                    hl_group = "@string.special.path",
                    ephemeral = true,
                })
            end
        end
    end
end

--- Register the decoration provider (markdown-parsed windows only) and the
--- autocmds that clear the existence cache. Idempotent.
function M.setup()
    vim.api.nvim_set_decoration_provider(ns, {
        on_win = function(_, _, bufnr)
            -- Gate on the treesitter language, not the raw filetype, so buffers
            -- that render as markdown under a different filetype (e.g.
            -- AgenticInput, which registers markdown) are included too.
            return vim.treesitter.language.get_lang(vim.bo[bufnr].filetype)
                == "markdown"
        end,
        on_range = on_range,
    })
    vim.api.nvim_create_autocmd({ "BufWritePost", "FocusGained" }, {
        group = vim.api.nvim_create_augroup("path_highlight", { clear = true }),
        desc = "path_highlight: invalidate filepath existence cache",
        callback = function() exists_cache = {} end,
    })
end

return M
