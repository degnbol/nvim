" Let's hide the comment leader by matching on it alone.
syntax match texCommentLeader '%' conceal
" define the syntax identifier texComment with regex for lines starting with %
syntax match texComment /%.*/ contains=texCommentLeader

highlight link texComment Comment
highlight link texCommentLeader Comment

