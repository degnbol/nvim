#!/usr/bin/env lua
-- https://github.com/DaikyXendo/nvim-material-icon
return {{
    "DaikyXendo/nvim-material-icon",
    opts = {
        -- your personnal icons can go here (to override)
        -- you can specify color or cterm_color instead of specifying both of them
        -- DevIcon will be appended to `name`
        override = {
            tsv = {
                icon = "",
                color = "#89e051",
                cterm_color = "113",
                name = "TabSeparatedValues"
            },
            tab = {
                icon = "",
                color = "#89e051",
                cterm_color = "113",
                name = "Table"
            },
            csv = {
                icon = "",
                color = "#89e051",
                cterm_color = "113",
                name = "CommaSeparatedValues"
            },
            mat = {
                icon = "󰘨",
                color = "#89e051",
                cterm_color = "113",
                name = "Matrix"
            },
            -- for some reason Markdown and md are different by default, copied here is making md the same as Markdown.
            md = {
                icon = "",
                color = "#519aba",
                cterm_color = "67",
                name = "Markdown",
            },
            adoc = {
                icon = "",
                color = "#519aba",
                cterm_color = "67",
                name = "Asciidoc",
            },
            crd = {
                icon = "󰝱",
                color = "#519aba",
                cterm_color = "67",
                name = "Chordpro",
            },
        };
        -- globally enable different highlight colors per icon (default to true)
        -- if set to false all icons will have the default icon's color
        -- color_icons = true;
        -- globally enable default icons (default to false)
        -- will get overriden by `get_icons` option
        -- default = true;
    }
},
{"nvim-tree/nvim-web-devicons", opts = {
    override = require'nvim-material-icon'.get_icons()
},}}
