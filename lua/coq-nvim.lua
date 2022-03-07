#!/usr/bin/env lua

vim.g.coq_settings = {
    ["keymap.recommended"] = false,
    ["clients.buffers.enabled"] = false,
    limits = {
        completion_auto_timeout = 1.0,
        completion_manual_timeout = 3.0
    },
    ["signature.target"] = "echo" -- disable signature help popup by default https://github.com/neoclide/coc.nvim/discussions/2778
}

vim.defer_fn(function() vim.opt.completeopt = 'menuone,noinsert' end, 1000) -- has to be deferred since COQ changes this setting

