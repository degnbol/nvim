-- Shell out to the `astfix` tree-sitter formatter on a temp copy of the buffer.
-- nvim holds zero ast-grep/rule knowledge; astfix keys rules by file extension.
local M = {}

--- Apply an astfix rule to the whole buffer (one undo step).
---@param rule string|nil rule name (rules/<ext>/<rule>.yml); nil → the ext's default
function M.run(rule)
  local ext = vim.fn.expand "%:e"
  if ext == "" then
    vim.notify("astfix: buffer has no file extension", vim.log.levels.ERROR)
    return
  end
  local tmp = vim.fn.tempname() .. "." .. ext
  vim.fn.writefile(vim.api.nvim_buf_get_lines(0, 0, -1, false), tmp)

  -- Build argv without embedded nils (table literals truncate at the first nil).
  local cmd = { "astfix" }
  if rule then
    cmd[#cmd + 1] = "-r"
    cmd[#cmd + 1] = rule
  end
  cmd[#cmd + 1] = "-U"
  cmd[#cmd + 1] = tmp

  local res = vim.system(cmd):wait()
  if res.code ~= 0 then
    vim.notify("astfix: " .. (res.stderr or "failed"), vim.log.levels.ERROR)
    return
  end
  vim.api.nvim_buf_set_lines(0, 0, -1, false, vim.fn.readfile(tmp))
end

--- Rule-name candidates for the current buffer's extension (for :AstFix completion).
---@return string[]
function M.complete()
  local res = vim.system({ "astfix", "--list", vim.fn.expand "%" }):wait()
  return vim.split(res.stdout or "", "\n", { trimempty = true })
end

return M
