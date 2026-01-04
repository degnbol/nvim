-- Custom treesitter highlighting for karabiner.json
-- Uses json parser with queries from queries/karabiner/*.scm

local read_query = require('utils').read_query

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
