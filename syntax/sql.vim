" postgres specific custom hls in this file. If we ever use another sql then 
" separate things.

syn match Statement '^\\copy'
syn keyword @keyword.modifier.sql header HEADER delimiter DELIMITER

" match string with escaped chars such as E'\t'
syn match @string.special /E\ze'/

