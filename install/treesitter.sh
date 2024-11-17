#!/usr/bin/env zsh
cd $0:h
# https://tree-sitter.github.io/tree-sitter/creating-parsers
# https://crates.io/crates/tree-sitter-cli
# also available on mason but might be the smaller version?
cargo install tree-sitter-cli

# LSP for tree-sitter query files
# https://github.com/ribru17/ts_query_ls
cd ../lsp/
wget 'https://github.com/ribru17/ts_query_ls/releases/download/v1.2.4/ts_query_ls-aarch64-apple-darwin.tar.gz'
tar xzf ts_query_ls*.tar.gz && rm ts_query_ls*.tar.gz
# looks like there is a bug in the code expecting the ts_query_ls to be in ~/.config/nvimlsp/ts_query_ls/ so we symlink it there as a quickfix
mkdir ~/.config/nvimlsp
ln -s $PWD/ts_query_ls ~/.config/nvimlsp/ts_query_ls
