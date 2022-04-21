#!/usr/bin/env lua
require("scrollbar").setup {
    show = false, -- enable with :ScrollbarToggle etc.
    marks = {
        Search = {
            color = "yellow"
        }
    },
    handlers = {
        diagnostics = false,
        search = true
    }
}
