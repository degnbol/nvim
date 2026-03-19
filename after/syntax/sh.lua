-- Disable built-in sh*Error groups — treesitter handles syntax errors now.
for _, g in ipairs{'shArithRegion', 'shCaseError', 'shCondError', 'shCurlyError',
    'shDerefOpError', 'shDerefWordError', 'shDoError', 'shEsacError',
    'shIfError', 'shInError', 'shParenError', 'shTestError', 'shDTestError'} do
  if vim.fn.hlexists(g) == 1 then
    vim.api.nvim_set_hl(0, g, {})
  end
end
