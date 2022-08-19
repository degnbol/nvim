#!/usr/bin/env lua
require"fidget".setup {
    text = {
        spinner = "star"
    },
    -- disable for specific sources that are annoying and verbose
    sources = {
        ltex = {ignore = true}
    }
}
