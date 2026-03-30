return {
    "conform.nvim",
    event = "BufWritePre",
    keys = {
        { "grf", function()
            require("conform").format({ async = true, lsp_fallback = true })
        end, mode = { "n", "x" }, desc = "Format" },
    },
    after = function()
        require("conform").setup({
            formatters_by_ft = {
                sh = { "shfmt" },
                bash = { "shfmt" },
                zsh = { "shfmt" },
                ["sh.zsh"] = { "shfmt" },
            },
            formatters = {
                shfmt = {
                    prepend_args = function(_, ctx)
                        local args = { "-i", "4", "-bn", "-ci", "-sr" }
                        if ctx.filename and ctx.filename:match("%.zsh$") or
                            vim.bo[ctx.buf].filetype:find("zsh") then
                            table.insert(args, 1, "-ln")
                            table.insert(args, 2, "zsh")
                        end
                        return args
                    end,
                },
            },
        })
    end,
}
