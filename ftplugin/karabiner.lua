-- Custom treesitter highlighting for karabiner.json
-- Uses json parser with extended highlights for karabiner-specific properties

-- Read base json highlights from runtime
local function get_query_source(lang, query_name)
  local files = vim.api.nvim_get_runtime_file('queries/' .. lang .. '/' .. query_name .. '.scm', true)
  local sources = {}
  -- Read in reverse order (last = highest priority in rtp)
  for i = #files, 1, -1 do
    local f = io.open(files[i], 'r')
    if f then
      table.insert(sources, f:read('*a'))
      f:close()
    end
  end
  return table.concat(sources, '\n')
end

local base_query = get_query_source('json', 'highlights')

local karabiner_patterns = [[
; Highlight karabiner-specific property names as @property.builtin
((pair
  key: (string
    (string_content) @_key) @property.builtin)
  (#any-of? @_key
    "global" "show_in_menu_bar"
    "rules" "description" "manipulators" "type" "from" "to"
    "key_code" "modifiers" "mandatory" "optional" "conditions" "is_keyboard"
    "frontmost_application_if" "frontmost_application_unless"
    "device_if" "device_unless" "keyboard_type_if" "keyboard_type_unless"
    "input_source_if" "input_sources" "input_source_id" "input_source_unless" "variable_if" "variable_unless"
    "event_changed_if" "event_changed_unless"
    "bundle_identifiers" "file_paths" "identifiers" "vendor_id" "product_id"
    "set_variable" "name" "value" "shell_command" "select_input_source"
    "mouse_key" "pointing_button" "consumer_key_code"
    "profiles" "selected" "parameters" "complex_modifications"
    "simple_modifications" "fn_function_keys" "devices"
    "virtual_hid_keyboard" "country_code" "keyboard_type_v2")
  (#set! priority 95))
]]

local extended_query = base_query .. '\n' .. karabiner_patterns

-- Stop any existing highlighter
vim.treesitter.stop()

-- Create parser and highlighter with custom query
local bufnr = vim.api.nvim_get_current_buf()
local parser = vim.treesitter.get_parser(bufnr, 'json')
vim.treesitter.highlighter.new(parser, {
  queries = { json = extended_query }
})
