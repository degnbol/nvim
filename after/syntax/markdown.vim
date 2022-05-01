" Let's hide the comment leader by matching on it alone.
syntax match texCommentLeader '%' conceal
" define the syntax identifier texComment with regex for lines starting with %
syntax match texComment /%.*/ contains=texCommentLeader

" color tex comment as whatever color we have set for comments in general
highlight link texComment Comment
highlight link texCommentLeader Comment
