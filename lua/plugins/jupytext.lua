#!/usr/bin/env lua
-- edit jupyter notebook. requires `pip install jupytext`
return {"goerz/jupytext.vim", config=function()
    -- conversion settings
    vim.cmd 'let g:jupytext_fmt = "py:percent"'
    -- highlight blocks
    vim.cmd 'match Block /^# %%/'
    vim.cmd 'hi link Block LineNr'
end}
