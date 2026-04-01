vim.filetype.add({
    pattern = {
        [".*%.blend%.py"] = "python.blender",
        [".*%.py"] = {
            function(_path, bufnr)
                local first = vim.api.nvim_buf_get_lines(bufnr, 0, 1, false)[1] or ""
                if first:match("blender") then
                    return "python.blender"
                end
            end,
            { priority = 10 },
        },
    },
})
