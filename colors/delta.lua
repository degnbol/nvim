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
hi.setbg("Visual", "#6a3b6a")

-- It's green by default which I associate with line added
hi.set("CursorLineNr", { fg = "white", bold = true })
