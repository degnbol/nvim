return {
    "conform.nvim",
    event = "BufWritePre",
    keys = {
        { "grf", function()
            local buf = vim.api.nvim_get_current_buf()
            require("conform").format(
                { async = true, lsp_format = "fallback", bufnr = buf },
                function(err) require("conform_error").report(err, buf) end
            )
        end, mode = { "n", "x" }, desc = "Format" },
},
    after = function()
        require("conform").setup({
            -- Execution-error stderr is surfaced by conform_error's callback
            -- (real message + inline diagnostic) instead of conform's generic
            -- "See :ConformInfo" notification.
            notify_on_error = false,
            formatters_by_ft = {
                sh = { "shfmt" },
                bash = { "shfmt" },
                zsh = { "shfmt", "fmt_zsh_mlr" },
                ["sh.zsh"] = { "shfmt", "fmt_zsh_mlr" },
                typst = { "typstyle" },
            },
            formatters = {
                -- Cosmetic mlr-invocation formatter (config/fmt-zsh-mlr). Runs
                -- after shfmt; the two are order-independent — shfmt leaves the
                -- `+ \` reflow byte-identical. Only active if the binary is on PATH.
                fmt_zsh_mlr = {
                    command = "fmt-zsh-mlr",
                    stdin = true,
                },
                shfmt = {
                    prepend_args = function(_, ctx)
                        local args = { "-i", "4", "-ci", "-sr" }
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
