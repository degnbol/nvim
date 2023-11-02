# Icons showing up as chinese characters
Terminal font is specified by kitty.conf.
It is currently referencing a nerd font patched Operator font.
Make sure it was installed with fonts/install.sh
Try reapplying kitty preferences (ctrl+cmd+,)
Make sure nerdfont is installed (install.sh)

# highlight groups acting weird when switching between different filetypes
The :hi command is global and a lot of highlight groups are shared, so editing 
a highlight group with one filetype in mind might cause unexpected changes for 
other filetypes. A solution is to not touch the main syntax groups and just 
link filetype specific groups to the right ones and make new ones. Another 
solution is to make :au BufEnter and :au BufLeave commands to disable changes 
whenever we leave the buffer that imposes them.

