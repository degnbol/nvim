-- :Canon [pretty|ascii] — canonicalize the buffer's source form (symbol
-- name↔glyph, punctuation). Global so it's reachable in any buffer; canon
-- errors cleanly on an unsupported extension. No argument → --pretty (default).
vim.api.nvim_create_user_command("Canon", function(opts)
  require("canon").run(opts.args ~= "" and opts.args or nil)
end, {
  nargs = "?",
  complete = function(arglead)
    return vim.tbl_filter(function(m)
      return m:find(arglead, 1, true) == 1
    end, require("canon").complete())
  end,
  desc = "Canonicalize the buffer (canon)",
})
