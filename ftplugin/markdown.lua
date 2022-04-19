-- markdown doesn't really have comments. 
-- The default is html comments <!-- ... -->.
-- For making tex documents with the markdown package it is possible to ignore latex comments, hence:
vim.api.nvim_buf_set_option(0, "commentstring", "% %s")
-- c=continue comment when wrapping. Should already be enabled.
-- r=with enter in insert mode
-- o=with o/O in normal mode
vim.opt.formatoptions:append "cro"
-- this also had to be set for the above to understand how a comment line looks.
vim.opt.comments = ":%"
