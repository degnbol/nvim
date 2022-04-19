-- markdown doesn't really have comments. 
-- The default is html comments <!-- ... -->.
-- For making tex documents with the markdown package it is possible to ignore latex comments, hence:
vim.api.nvim_buf_set_option(0, "commentstring", "% %s")
