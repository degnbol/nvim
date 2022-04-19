" Let's hide the comment leader by matching on it alone.
syntax match texCommentLeader '%' conceal
" define the syntax identifier texComment with regex for lines starting with %
syntax match texComment /%.*/ contains=texCommentLeader
" color tex comment as whatever color we have set for comments in general
highlight link texComment Comment
highlight link texCommentLeader Comment

" conceal comment leader and the _ and * around emphasis
" level=1 -> conceal but don't remove block
set conceallevel=1
" nvc=normal,visual,command -> only unconceal in insert mode with cursor on 
" the line.
set concealcursor=nvc

