local map = require "utils/keymap"

map.n('<leader>cc',     ":Crun\n",   "Cargo run", {buf=0} )
map.n('<LocalLeader>r', ":Crun\n",   "Cargo run", {buf=0} )
map.n('<LocalLeader>b', ":Cbuild\n", "Cargo build", {buf=0} )
map.n('<LocalLeader>c', ":Ccheck\n", "Cargo check", {buf=0} )
map.n('<LocalLeader>C', ":Cclean\n", "Cargo clean", {buf=0} )
map.n('<LocalLeader>t', ":Ctest\n",  "Cargo test", {buf=0} )
map.n('<LocalLeader>u', ":Cupdate\n","Cargo update", {buf=0} )

map.i(';;', ";<Esc>", "Escape at EOL", {buf=0, remap=false})
