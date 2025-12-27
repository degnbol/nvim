-- Custom treesitter highlighting for karabiner.json
-- Uses json parser with queries from queries/karabiner/*.scm

-- Read query file, concatenating all matches from runtimepath
local function read_query(lang, query_name)
  local files = vim.api.nvim_get_runtime_file('queries/' .. lang .. '/' .. query_name .. '.scm', true)
  local sources = {}
  for i = #files, 1, -1 do
    local f = io.open(files[i], 'r')
    if f then
      table.insert(sources, f:read('*a'))
      f:close()
    end
  end
  return table.concat(sources, '\n')
end

local base_highlights = read_query('json', 'highlights')
local karabiner_highlights = read_query('karabiner', 'highlights')
local karabiner_injections = read_query('karabiner', 'injections')

-- Stop any existing highlighter
vim.treesitter.stop()

-- Create parser and highlighter with combined queries
local bufnr = vim.api.nvim_get_current_buf()
local parser = vim.treesitter.get_parser(bufnr, 'json')
vim.treesitter.highlighter.new(parser, {
  queries = {
    json = base_highlights .. '\n' .. karabiner_highlights,
  }
})

-- Set up injections (global, but pattern is karabiner-specific)
if karabiner_injections ~= '' then
  vim.treesitter.query.set('json', 'injections', karabiner_injections)
end
