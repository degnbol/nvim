#!/usr/bin/env lua
local make_default_opts = function()
  return {
    final_separator = ",",
    final_end_line = true,
    skip = {},
  }
end

local make_no_final_sep_opts = function()
  return {
    final_separator = false,
    final_end_line = true,
  }
end

require'trevj'.setup {
    containers = {
        julia = {
            parameters = make_default_opts(),
            argument_list = make_default_opts(),
            list = make_default_opts(),
            tuple = make_default_opts(),
            dictionary = make_default_opts(),
            set = make_default_opts(),
            list_comprehension = make_no_final_sep_opts(),
            generator_expression = make_no_final_sep_opts(),
            dictionary_comprehension = make_no_final_sep_opts(),
        }
    }
}
