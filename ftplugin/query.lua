#!/usr/bin/env lua

-- https://github.com/ribru17/ts_query_ls
local function LSP_start(cmdpath)
    if vim.bo.buftype == 'nofile' then return end
    vim.lsp.start {
        name = 'ts_query_ls',
        cmd = { cmdpath },
        root_dir = vim.fs.root(0, { 'queries' }),
        settings = {
            parser_install_directories = {
                -- If using nvim-treesitter with lazy.nvim
                vim.fs.joinpath(
                    vim.fn.stdpath('data'),
                    '/lazy/nvim-treesitter/parser/'
                ),
            },
            parser_aliases = {
                ecma = 'javascript',
            },
            language_retrieval_patterns = {
                'languages/src/([^/]+)/[^/]+\\.scm$',
            },
        },
    }
end

local lsppath = vim.opt.runtimepath:get()[1] .. 'lsp/ts_query_ls'
if vim.fn.filereadable(lsppath) then
    LSP_start(lsppath)
end

