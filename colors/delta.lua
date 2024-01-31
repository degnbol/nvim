require "fluoromachine".setup {
    glow = true,
    theme = 'delta',
    overrides = {
        ['Comment'] = { italic = false },
    }
}
-- visual selection is not clear enough.
vim.api.nvim_set_hl(0, "Visual", {bg="#6a3b6a"})
