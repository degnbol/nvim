" Highlights doesn't seem to take effect if written in syntax/ or after/syntax/.
" It might be due to some double loading of some parts of config.
" Placing `hi ...` here works.

" In vim9 comments can be prefixed by "#", but neovim has locked its vimscript 
" for config before this was allowed. Since we only use vimscript to configure 
" neovim we show it as error.
hi default link vim9comment Error

