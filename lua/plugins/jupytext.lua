#!/usr/bin/env lua
-- edit jupyter notebook. requires `pip install jupytext`
return {
    "goerz/jupytext.vim",
    config=function()
        -- conversion settings
        vim.g.jupytext_fmt = "py:percent"
        -- highlight blocks.
        -- Doesn't work since treesitter trumps. Need to add a custom query.
        vim.cmd 'syn match Block /^# %%/'
        vim.api.nvim_set_hl(0, "Block", {link="LineNr"})
    end,
}
