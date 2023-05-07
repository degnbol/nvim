#!/usr/bin/env lua
-- treesitter colors in as @keyword.operator, but I would rather it be considered @repeat to match 'for'
-- note the wonderful ability to add ".julia" at the end and it only changes locally to julia and 
-- not other languages as opposed to standard hightlight groups
-- vim.api.nvim_set_hl(0, "@keyword.operator.julia", {link="@repeat"})
-- not needed since we change "in" to @repeat with the after/queries/julia/hightlights.scm that is copied into queries/julia/hightlights.scm

