-- In-process LSP providing hover for glossy `@term` refs in typst.
-- No external server — cmd is a Lua function (cf. modules/kitty-conf.nvim).
-- tinymist resolves glossy refs to the package's runtime `label()` placeholder
-- and returns no hover for them, so this joins the standard multi-client hover
-- (`vim.lsp.buf.hover`) as the sole provider on glossary keys. All resolution
-- lives in lua/typst_glossary.lua, shared with the `grd` goto-def path.

local gloss = require "typst_glossary"

return {
    cmd = function()
        return {
            request = function(method, params, callback)
                if method == "initialize" then
                    callback(nil, { capabilities = { hoverProvider = true } })
                elseif method == "textDocument/hover" then
                    local bufnr = vim.uri_to_bufnr(params.textDocument.uri)
                    local markdown = gloss.hover(bufnr, params.position.line, params.position.character)
                    if markdown then
                        callback(nil, { contents = { kind = "markdown", value = markdown } })
                    else
                        callback(nil, nil)
                    end
                elseif method == "shutdown" then
                    callback(nil, nil)
                end
            end,
            notify = function() end,
            is_closing = function() return false end,
            terminate = function() end,
        }
    end,
    filetypes = { "typst" },
}
