-- Shell out to `canon` (the source-form canonicalizer) over the whole buffer.
-- canon dispatches by file extension; nvim only picks the mode. Streams the
-- buffer through --stdin-filepath so no temp file is needed and the on-disk
-- file is never touched.
local M = {}

--- Canonicalize the whole buffer (one undo step).
---@param mode "pretty"|"ascii"|nil transform mode; nil → canon's default (--pretty)
function M.run(mode)
  local name = vim.api.nvim_buf_get_name(0)
  if name == "" then
    vim.notify("canon: buffer has no file name (for extension detection)", vim.log.levels.ERROR)
    return
  end

  local cmd = { "canon" }
  if mode then
    cmd[#cmd + 1] = "--" .. mode
  end
  cmd[#cmd + 1] = "--stdin-filepath"
  cmd[#cmd + 1] = name

  local input = table.concat(vim.api.nvim_buf_get_lines(0, 0, -1, false), "\n") .. "\n"
  local res = vim.system(cmd, { stdin = input }):wait()
  if res.code ~= 0 then
    -- --ascii tier-3 (no total-ASCII form) exits non-zero with the sites on stderr.
    vim.notify("canon: " .. (res.stderr ~= "" and res.stderr or "failed"), vim.log.levels.ERROR)
    return
  end

  local lines = vim.split(res.stdout, "\n", { plain = true })
  if lines[#lines] == "" then
    lines[#lines] = nil -- drop the trailing newline's empty tail
  end
  vim.api.nvim_buf_set_lines(0, 0, -1, false, lines)
end

--- Mode candidates for :Canon completion.
---@return string[]
function M.complete()
  return { "pretty", "ascii" }
end

return M
