local M = {}

--- Whether a delimited LaTeX math source is inline (`$…$`, `$`…`$`, `\(…\)`)
--- rather than display (`$$…$$`, `\[…\]`).
--- @param src string math source including its delimiters, surrounding whitespace allowed
--- @return boolean
function M.is_inline(src)
    local s = vim.trim(src)
    return s:find("^%$[^$]") ~= nil or s:find("^\\%(") ~= nil
end

--- The body of a delimited LaTeX math source (`$…$`, `$`…`$`, `$$…$$`,
--- `\(…\)`, `\[…\]`). Whitespace outside the delimiters is trimmed, whitespace
--- inside is kept.
--- @param src string math source including its delimiters
--- @return string body
function M.math_body(src)
    local body = vim.trim(src):gsub("^%$+`?", ""):gsub("`?%$+$", "")
        :gsub("^\\[%[%(]", ""):gsub("\\[%]%)]$", "")
    return body
end

return M
