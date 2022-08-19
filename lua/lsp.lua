#!/usr/bin/env lua
-- default config copied from https://github.com/neovim/nvim-lspconfig

-- Mappings.
-- See `:help vim.diagnostic.*` for documentation on any of the below functions
local opts = { noremap=true, silent=true }
-- vim.keymap.set('n', '<space>E', vim.diagnostic.open_float, opts)
vim.keymap.set('n', '<space>E', "<cmd>Telescope diagnostics<cr>", opts)
vim.keymap.set('n', '[d', vim.diagnostic.goto_prev, opts)
vim.keymap.set('n', ']d', vim.diagnostic.goto_next, opts)
vim.keymap.set('n', '<space>q', vim.diagnostic.setloclist, opts)

-- Use an on_attach function to only map the following keys
-- after the language server attaches to the current buffer
local on_attach = function(client, bufnr)
  -- Enable completion triggered by <c-x><c-o>
  vim.api.nvim_buf_set_option(bufnr, 'omnifunc', 'v:lua.vim.lsp.omnifunc')

  -- Mappings.
  -- See `:help vim.lsp.*` for documentation on any of the below functions
  local bufopts = { noremap=true, silent=true, buffer=bufnr }
  vim.keymap.set('n', 'gD', vim.lsp.buf.declaration, bufopts)
  vim.keymap.set('n', 'gd', vim.lsp.buf.definition, bufopts)
  vim.keymap.set('n', 'K', vim.lsp.buf.hover, bufopts)
  vim.keymap.set('n', 'gi', vim.lsp.buf.implementation, bufopts)
  vim.keymap.set('n', '<C-k>', vim.lsp.buf.signature_help, bufopts)
  vim.keymap.set('n', '<space>wa', vim.lsp.buf.add_workspace_folder, bufopts)
  vim.keymap.set('n', '<space>wr', vim.lsp.buf.remove_workspace_folder, bufopts)
  vim.keymap.set('n', '<space>wl', function() print(vim.inspect(vim.lsp.buf.list_workspace_folders())) end, bufopts)
  vim.keymap.set('n', '<space>D', vim.lsp.buf.type_definition, bufopts)
  vim.keymap.set('n', '<space>rn', vim.lsp.buf.rename, bufopts)
  vim.keymap.set('n', '<space>ca', vim.lsp.buf.code_action, bufopts)
  -- vim.keymap.set('n', 'gr', vim.lsp.buf.references, bufopts)
  vim.keymap.set('n', 'gr', "<cmd>Telescope lsp_references<CR>", bufopts)
  vim.keymap.set('n', '<space>F', vim.lsp.buf.format, bufopts)
end

-- replace the default lsp diagnostic letters with prettier symbols
vim.fn.sign_define("LspDiagnosticsSignError", {text = "", numhl = "LspDiagnosticsDefaultError"})
vim.fn.sign_define("LspDiagnosticsSignWarning", {text = "", numhl = "LspDiagnosticsDefaultWarning"})
vim.fn.sign_define("LspDiagnosticsSignInformation", {text = "", numhl = "LspDiagnosticsDefaultInformation"})
vim.fn.sign_define("LspDiagnosticsSignHint", {text = "", numhl = "LspDiagnosticsDefaultHint"})

-- hide diagnostics for hints and information.
vim.lsp.handlers["textDocument/publishDiagnostics"] = vim.lsp.with(
  vim.lsp.diagnostic.on_publish_diagnostics, {
    signs = {
      severity_limit = 'Warning',
    },
    underline = false,
    update_in_insert = false,
    virtual_text = {
      spacing = 4,
      severity_limit = 'Error',
    },
  }
)

-- color kinds
vim.api.nvim_set_hl(0, 'CmpItemAbbrDeprecated', {fg="#808080", strikethrough=true})
vim.api.nvim_set_hl(0, 'CmpItemAbbrMatch', {fg="#569CD6"})
vim.api.nvim_set_hl(0, 'CmpItemAbbrMatchFuzzy', {fg="#569CD6"})
vim.api.nvim_set_hl(0, 'CmpItemKindVariable', {fg="#9CDCFE"})
vim.api.nvim_set_hl(0, 'CmpItemKindInterface', {fg="#9CDCFE"})
vim.api.nvim_set_hl(0, 'CmpItemKindText', {fg="#9CDCFE"})
vim.api.nvim_set_hl(0, 'CmpItemKindFunction', {fg="#C586C0"})
vim.api.nvim_set_hl(0, 'CmpItemKindMethod', {fg="#C586C0"})
vim.api.nvim_set_hl(0, 'CmpItemKindKeyword', {fg="#D4D4D4"})
vim.api.nvim_set_hl(0, 'CmpItemKindProperty', {fg="#D4D4D4"})
vim.api.nvim_set_hl(0, 'CmpItemKindUnit', {fg="#D4D4D4"})

-- now for adding the language servers
local lsp = require "lspconfig"
-- coq for speed
-- local coq = require "coq"
-- cmp for functional, customizable and easy to control
local capabilities = require('cmp_nvim_lsp').update_capabilities(vim.lsp.protocol.make_client_capabilities())

-- naming: https://github.com/williamboman/mason-lspconfig.nvim/blob/main/doc/server-mapping.md
-- config help: https://github.com/neovim/nvim-lspconfig/blob/master/doc/server_configurations.md

lsp.bashls.setup { filetypes = {"sh", "bash", "zsh"}, on_attach=on_attach, capabilities=capabilities }

-- lsp.pyright.setup { on_attach=on_attach, capabilities=capabilities }
-- lsp.pylsp.setup { on_attach=on_attach, capabilities=capabilities }
lsp.jedi_language_server.setup { on_attach=on_attach, capabilities=capabilities }
-- lsp.jedi_language_server.setup(coq.lsp_ensure_capabilities { on_attach=on_attach })

lsp.julials.setup { on_attach=on_attach, capabilities=capabilities }

lsp.ltex.setup { on_attach=function(client, bufnr)
    on_attach(client, bufnr)
    -- hacky. VimtexErrors puts errors found by Vimtex in quickfix (should be 
    -- running, use <leader>Lb) then cclose closes quickfix, and then Telescope 
    -- opens the quickfix in a nicer view.
    vim.keymap.set('n', '<space>E', "<cmd>VimtexErrors<cr>|:cclose|<cmd>Telescope quickfix<cr>", opts)
end, capabilities=capabilities }

lsp.r_language_server.setup { on_attach=on_attach, capabilities=capabilities }

lsp.vimls.setup { on_attach=on_attach, capabilities=capabilities }

lsp.sumneko_lua.setup { on_attach=on_attach, capabilities=capabilities }

