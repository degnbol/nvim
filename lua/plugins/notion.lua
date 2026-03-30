
return {
    {
        "notion.nvim",
        enabled = false,
        after = function()
            require"notion".setup()
        end,
    },
}
