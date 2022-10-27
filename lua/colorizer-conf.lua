#!/usr/bin/env lua
require("colorizer").setup {
    filetypes = {
        julia = {
            RGB = false, -- #RGB hex codes are too simple and may show up in julia in rare cases unnamed variables
        }
    }
}

