#!/usr/bin/env lua
-- corner print what LSP is running
-- https://github.com/j-hui/fidget.nvim/blob/main/doc/fidget.md
return {
    "j-hui/fidget.nvim",
    opts = {
        text = {
            spinner = "star",
            commenced = "",
            completed = "",
        },
        fmt = {
            -- function to format fidget title
            fidget = function(fidget_name, spinner)
                -- hide which lsp is being loaded.
                return spinner
                -- default:
                -- return string.format("%s %s", spinner, fidget_name)
            end,
            -- function to format each task line
            task = function(task_name, message, percentage)
                -- by hiding task info we only see a minimal overall indication (the "fidget title" part)
                return ""
                -- default:
                -- return string.format("%s%s [%s]",
                -- message,
                -- percentage and string.format(" (%s%%)", percentage) or "",
                -- task_name)
            end,
        },
        -- disable for specific sources that are annoying and verbose
        sources = {
            -- ltex = {ignore = true},
        }
    }
}
