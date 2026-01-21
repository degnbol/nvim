
local opt = vim.opt_local
-- markdown doesn't really have comments. 
-- The default is html comments <!-- ... -->.
-- For making tex documents with the markdown package it is possible to ignore latex comments, hence:
opt.commentstring = "% %s"
-- c=continue comment when wrapping. Should already be enabled.
-- r=with enter in insert mode. Should already be enabled.
-- o=with o/O in normal mode. Not sure if I want this.
-- a=auto format by default
-- v=inserted now, don't edit other lines.
opt.formatoptions:append "crv"
-- this also had to be set for the above to understand how a comment line looks.
-- the b means the % should always be followed by a blank, i.e. space
opt.comments = "b:%"
-- ai is necessary for lists but a problem if I try to indent start of sentence, then the next lines are also indented.
-- opt.autoindent = false
-- since we are actively using the textwidth in markdown then the moving of the window is just annoying for a slim window
opt.sidescrolloff = 0

-- conceal comment leader and the _ and * around emphasis
-- level=1 -> conceal but don't remove block
opt.conceallevel = 1

opt.list = false
