-- markdown doesn't really have comments. 
-- The default is html comments <!-- ... -->.
-- For making tex documents with the markdown package it is possible to ignore latex comments, hence:
vim.api.nvim_buf_set_option(0, "commentstring", "% %s")
-- c=continue comment when wrapping. Should already be enabled.
-- r=with enter in insert mode
-- o=with o/O in normal mode
-- a=auto format by default
vim.opt.formatoptions:append "croa"
-- this also had to be set for the above to understand how a comment line looks.
-- the b means the % should always be followed by a blank, i.e. space
vim.opt.comments = "b:%"
-- ai is useful for lists but a problem if I try to indent start of sentence, then the next lines are also indented.
vim.opt.autoindent = false
