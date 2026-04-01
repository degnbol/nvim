return {
    {
        "nightfox.nvim",
        colorscheme = { "nightfox", "dayfox", "dawnfox", "duskfox", "nordfox", "terafox", "carbonfox" },
        after = function()
            require("nightfox").setup({
                options = {
                    styles = {
                        keywords = "italic",
                        operators = "bold",
                    },
                },
                groups = {
                    -- In carbonfox there's a lot of shades of blue being used and builtins were colored "red" from the palette.
                    -- We instead are using italic to signify builtin and matched colour on the builtins to the non-builtin equivalent group.
                    -- This means we have red from the palette freed up to show
                    -- something else so we change all the special hl groups to
                    -- red from a blue that is used for function.
                    carbonfox = {
                        Special = { fg = "palette.red" },
                    },
                    -- For terafox the red looks too aggressive and there is an underutilized pink.
                    terafox = {
                        Special = { fg = "palette.pink" },
                        -- @tag has color like keyword. We link it to Tag which is linked to Special which makes it pink.
                        ["@tag"] = { link = "Tag" },
                    }
                }
            })
        end
    },
    {
        "fluoromachine.nvim",
        colorscheme = { "fluoromachine", "delta", "retrowave" },
    },
}
