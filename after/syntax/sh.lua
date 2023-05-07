#!/usr/bin/env lua
-- clear variable color inside string
vim.api.nvim_set_hl(0, "@variable.bash", {fg="white"})

-- waste of limited set of colors to destinguish between different function 
-- calls based on how core they are, and it's still bold in bash so there's still a hint anyways.
vim.api.nvim_set_hl(0, "@function.builtin.bash", {link="@function.call"})

vim.api.nvim_set_hl(0, "@path", {underline=true})
vim.api.nvim_set_hl(0, "@flag", {link="Special"})

