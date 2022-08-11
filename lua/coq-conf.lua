#!/usr/bin/env lua
vim.g.coq_settings = {
    ["keymap.recommended"] = false,
    ["clients.buffers.enabled"] = false,
    limits = {
        completion_auto_timeout = 1.0,
        completion_manual_timeout = 3.0
    }
}

