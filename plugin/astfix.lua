-- :AstFix [rule] — apply an astfix tree-sitter-formatter rule to the buffer.
-- Global so the entry point is discoverable everywhere; in a buffer whose
-- extension has no rules, completion is empty and an action errors cleanly.
vim.api.nvim_create_user_command("AstFix", function(opts)
  require("astfix").run(opts.args ~= "" and opts.args or nil)
end, {
  nargs = "?",
  complete = function()
    return require("astfix").complete()
  end,
  desc = "Apply an astfix rule to the buffer",
})
