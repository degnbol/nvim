#!/usr/bin/env lua
-- conversion settings
vim.cmd 'let g:jupytext_fmt = "py:percent"'
-- highlight blocks
vim.cmd 'match Block /^# %%/'
vim.cmd 'hi link Block LineNr'
