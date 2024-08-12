require "fluoromachine".setup {
    glow = true,
    theme = 'fluoromachine',
    overrides = {
        ['Comment'] = { italic = false },
    }
}
require('fluoromachine.config').load()
