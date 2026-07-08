local util = require "utils/init"
local map = require "utils/keymap"

-- assume using typst.vim
map.n('<leader>cc', "<Cmd>TypstWatch<CR>", "Compile continuously", { buffer=true, })
local path_pdf = vim.api.nvim_buf_get_name(0):gsub(".typ", ".pdf")
if util.is_mac() then
    map.n('<leader>cv', "!open -a skim " .. path_pdf .. "<CR>", "Compile view", { buffer=true, silent=true, })
else
    map.n('<leader>cv', "!pdf " .. path_pdf .. "<CR>", "Compile view", { buffer=true, silent=true, })
end

-- grd: jump to a glossary entry's definition for glossy `@term` refs; defer to
-- the LSP for everything else (`@fig:`/`@sec:`/`@eq:`, which it resolves right).
map.n('grd', function()
    local items = require"typst_glossary".resolve(0)
    if items then map.qf_mini { items = items } else map.lsp_definition() end
end, "Definition (glossary-aware)", { buffer=true })

-- K: our concise glossary hover for glossy `@term` refs; defer to the LSP
-- otherwise. tinymist doesn't strip `:pl`/`:cap`/… ref modifiers, so its hover
-- for a glossy ref is noisier and less useful — ours wins outright rather than
-- both showing (vim.lsp.buf.hover merges every client with no priority knob).
map.n('K', function()
    local row, col = unpack(vim.api.nvim_win_get_cursor(0))
    local markdown = require"typst_glossary".hover(0, row - 1, col)
    if markdown then
        vim.lsp.util.open_floating_preview(
            vim.split(markdown, "\n"), "markdown", { focus_id = "textDocument/hover" })
    else
        vim.lsp.buf.hover()
    end
end, "Hover (glossary-aware)", { buffer=true })

-- transliterate math symbol names → unicode glyphs, whole buffer (default rule).
map.n('<leader>tf', function() require("astfix").run() end, "astfix: names→glyphs", { buffer=true })

vim.opt_local.wrap = true
vim.opt_local.sidescrolloff = 0
-- definitely don't break at @ in typst (when wrapping with linebreak). 
-- It is used right before citations and is not a place to break.
vim.opt_local.breakat:remove('@')
