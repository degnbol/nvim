return {
    {
        "codedocs.nvim",
        after = function()
            -- Strip the `(type)` from Google-style Args lines: the types are
            -- already visible in the annotated signature.
            local blocks =
                vim.deepcopy(require("codedocs.config").opts.languages.python.styles.Google.annots.func.blocks)
            blocks[2].items[1].layout[1] = "%>%>%item_name: ${%snip_idx:description}"

            require"codedocs".setup {
                languages = {
                    python = {
                        default_style = "Google",
                        styles = {
                            Google = { annots = { func = { blocks = blocks } } },
                        },
                    },
                },
            }
        end,
    }
}
