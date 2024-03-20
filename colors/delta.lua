require "fluoromachine".setup {
    glow = true,
    theme = 'delta',
    overrides = {
        ['Comment'] = { italic = false },
    }
}
-- visual selection is not clear enough.
vim.api.nvim_set_hl(0, "Visual", {bg="#6a3b6a"})
-- default is the same color as functions, but @string.special is the only 
-- special hl with a different color.
vim.api.nvim_set_hl(0, "Special", {link="@string.special"})
