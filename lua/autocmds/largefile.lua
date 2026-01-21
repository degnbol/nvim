-- Large file handling for files opened within a session
local lf = require("largefile")
local group = vim.api.nvim_create_augroup("largefile", { clear = true })
local saved_eventignore = nil

-- Intercept at BufAdd (earliest point where filename is known)
vim.api.nvim_create_autocmd("BufAdd", {
    group = group,
    callback = function(args)
        local file = vim.api.nvim_buf_get_name(args.buf)
        if file == "" then return end

        local ok, stats = pcall(vim.uv.fs_stat, file)
        if not ok or not stats or stats.size <= lf.extreme_threshold then return end

        vim.b[args.buf].largefile = true
        vim.b[args.buf].largefile_size = stats.size
        saved_eventignore = vim.o.eventignore
        vim.o.eventignore = "FileType,Syntax,BufReadPost,BufEnter"
    end,
})

-- Restore eventignore after file is displayed
vim.api.nvim_create_autocmd("BufWinEnter", {
    group = group,
    callback = function(args)
        if saved_eventignore == nil then return end
        if not vim.b[args.buf].largefile then return end

        vim.o.eventignore = saved_eventignore
        saved_eventignore = nil

        lf.configure(args.buf)
        local size = vim.b[args.buf].largefile_size or 0
        vim.notify(
            string.format("Large file (%.0f MB) — loaded", size / 1024 / 1024),
            vim.log.levels.WARN
        )
    end,
})

-- Fallback for moderate files (1-50 MB)
vim.api.nvim_create_autocmd("BufReadPre", {
    group = group,
    callback = function(args)
        if vim.b[args.buf].largefile then return end

        local file = args.file
        local ok, stats = pcall(vim.uv.fs_stat, file)
        if not ok or not stats or stats.size <= lf.threshold then return end

        vim.b[args.buf].largefile = true
        vim.bo[args.buf].swapfile = false
        vim.bo[args.buf].undofile = false

        if stats.size > lf.extreme_threshold and saved_eventignore == nil then
            saved_eventignore = vim.o.eventignore
            vim.o.eventignore = "FileType,Syntax,BufReadPost,BufEnter"
            vim.b[args.buf].largefile_size = stats.size
        end
    end,
})

-- For moderate large files: disable features after load
vim.api.nvim_create_autocmd("BufReadPost", {
    group = group,
    callback = function(args)
        if not vim.b[args.buf].largefile then return end

        local file = args.file
        local ok, stats = pcall(vim.uv.fs_stat, file)
        if ok and stats and stats.size > lf.extreme_threshold then return end

        lf.configure(args.buf)
        if ok and stats then
            vim.notify(
                string.format("Large file (%.1f MB) — disabled features", stats.size / 1024 / 1024),
                vim.log.levels.WARN
            )
        end
    end,
})

-- Prevent LSP from attaching to large files
vim.api.nvim_create_autocmd("LspAttach", {
    group = group,
    callback = function(args)
        if vim.b[args.buf].largefile then
            vim.schedule(function()
                vim.lsp.buf_detach_client(args.buf, args.data.client_id)
            end)
        end
    end,
})
