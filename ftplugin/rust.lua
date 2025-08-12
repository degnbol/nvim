local map = require "utils/keymap"

map.n('<leader>cc',     ":Crun\n",   "Cargo run", {buffer=true} )
map.n('<LocalLeader>r', ":Crun\n",   "Cargo run", {buffer=true} )
map.n('<LocalLeader>b', ":Cbuild\n", "Cargo build", {buffer=true} )
map.n('<LocalLeader>c', ":Ccheck\n", "Cargo check", {buffer=true} )
map.n('<LocalLeader>C', ":Cclean\n", "Cargo clean", {buffer=true} )
map.n('<LocalLeader>t', ":Ctest\n",  "Cargo test", {buffer=true} )
map.n('<LocalLeader>u', ":Cupdate\n","Cargo update", {buffer=true} )

map.i(';;', ";<Esc>", "Escape at EOL", {buffer=true, remap=false})
