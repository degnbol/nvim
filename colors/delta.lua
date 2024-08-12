local hi = require "utils/highlights"

require "fluoromachine".setup {
    glow = true,
    theme = 'delta',
    overrides = {
        ['Comment'] = { italic = false },
    }
}
require('fluoromachine.config').load()

-- visual selection is not clear enough.
hi.bg("Visual", "#6a3b6a")
-- default is the same color as functions, but @string.special is the only 
-- special hl with a different color.
hi.link("Special", "@string.special")

-- It's green by default which I associate with line added
hi.set("CursorLineNr", {fg="white", bold=true})

