opt = vim.opt
-- markdown doesn't really have comments. 
-- The default is html comments <!-- ... -->.
-- For making tex documents with the markdown package it is possible to ignore latex comments, hence:
vim.api.nvim_buf_set_option(0, "commentstring", "% %s")
-- c=continue comment when wrapping. Should already be enabled.
-- r=with enter in insert mode
-- o=with o/O in normal mode
-- a=auto format by default
opt.formatoptions:append "croa"
-- this also had to be set for the above to understand how a comment line looks.
-- the b means the % should always be followed by a blank, i.e. space
opt.comments = "b:%"
-- ai is useful for lists but a problem if I try to indent start of sentence, then the next lines are also indented.
opt.autoindent = false
-- since we are actively using the textwidth in markdown then the moving of the window is just annoying for a slim window
opt.sidescrolloff = 0
