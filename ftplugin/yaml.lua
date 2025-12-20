local map = require "utils/keymap"

-- Sort glossary entries by their `short` field
-- Useful for typst glossy YAML files
local function sort_glossary_by_short()
  local lines = vim.api.nvim_buf_get_lines(0, 0, -1, false)

  -- Parse YAML entries (top-level keys with nested short/long/description)
  local entries = {}
  local current_key = nil
  local current_lines = {}

  for _, line in ipairs(lines) do
    -- Skip empty lines and comments at the top
    if line:match("^%s*$") or line:match("^#") then
      if current_key then
        table.insert(current_lines, line)
      end
    -- Top-level key (no leading whitespace, ends with colon)
    elseif line:match("^%w[%w_-]*:%s*$") then
      -- Save previous entry
      if current_key then
        entries[current_key] = current_lines
      end
      current_key = line:match("^(%w[%w_-]*):")
      current_lines = { line }
    -- Nested content (has leading whitespace)
    elseif line:match("^%s+") and current_key then
      table.insert(current_lines, line)
    end
  end
  -- Save last entry
  if current_key then
    entries[current_key] = current_lines
  end

  -- Extract short names for sorting
  local keys_with_short = {}
  for key, entry_lines in pairs(entries) do
    local short_name = key -- fallback to key if no short field
    for _, line in ipairs(entry_lines) do
      local short = line:match("^%s+short:%s*(.+)%s*$")
      if short then
        short_name = short
        break
      end
    end
    table.insert(keys_with_short, { key = key, short = short_name:lower() })
  end

  -- Sort by short name
  table.sort(keys_with_short, function(a, b)
    return a.short < b.short
  end)

  -- Rebuild buffer
  local new_lines = {}
  for _, item in ipairs(keys_with_short) do
    for _, line in ipairs(entries[item.key]) do
      table.insert(new_lines, line)
    end
  end

  vim.api.nvim_buf_set_lines(0, 0, -1, false, new_lines)
  vim.notify("Sorted " .. #keys_with_short .. " glossary entries by short name", vim.log.levels.INFO)
end

vim.api.nvim_buf_create_user_command(0, "SortGlossary", sort_glossary_by_short, {
  desc = "Sort glossary entries by short name",
})

map.n("<localLeader>s", sort_glossary_by_short, "Sort typst glossy entries by \"short\"")

