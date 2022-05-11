base16 = require'base16'
theme = base16.theme_from_array(require("themes/gigavoltArray"))
-- theme = base16.themes["unikitty-dark"]
base16(theme, true)
