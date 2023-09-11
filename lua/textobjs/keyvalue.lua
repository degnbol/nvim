#!/usr/bin/env lua
local M = {}

---@param inner boolean inner value excludes trailing commas or semicolons, outer includes them. Both exclude trailing comments.
---@param lookForwL integer number of lines to look forward for the textobj
function M.value(inner, lookForwL)
	-- captures value till the end of the line
	-- negative sets to not find equality comparators == or css pseudo-elements ::
	local pattern = "([^=:][=:] ?)[^=:].*()"

	local valueFound = searchTextobj(pattern, true, lookForwL)
	if not valueFound then return end

	-- if value found, remove trailing comment from it
	local curRow = fn.line(".")
	local lineContent = u.getline(curRow)
	if bo.commentstring ~= "" then -- JSON has empty commentstring
		local commentPat = bo.commentstring:gsub(" ?%%s.*", "") -- remove placeholder and backside of commentstring
		commentPat = vim.pesc(commentPat) -- escape lua pattern
		commentPat = " *" .. commentPat .. ".*" -- to match till end of line
		lineContent = lineContent:gsub(commentPat, "") -- remove commentstring
	end
	local valueEndCol = #lineContent - 1

	-- inner value = exclude trailing comma/semicolon
	if inner and lineContent:find("[,;]$") then valueEndCol = valueEndCol - 1 end

	u.setCursor(0, { curRow, valueEndCol })
end

---@param inner boolean outer key includes the `:` or `=` after the key
---@param lookForwL integer number of lines to look forward for the textobj
function M.key(inner, lookForwL)
	local pattern = "(%s*).-( ?[:=] ?)"

	local valueFound = searchTextobj(pattern, inner, lookForwL)
	if not valueFound then return end

	-- 1st capture is included for the outer obj, but we don't want it
	if not inner then
		local curRow = fn.line(".")
		local leadingWhitespace = u.getline(curRow):find("[^%s]") - 1
		u.normal("o")
		u.setCursor(0, { curRow, leadingWhitespace })
	end
end

return M
