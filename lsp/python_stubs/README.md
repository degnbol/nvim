Getting completion for some packages may take more effort if they are compiled 
C modules.
I checked and the completion was also missing from VS code so it's not a thing 
that gets supported out-of-the-box.
Python stubs are sometimes generated from python source code. Some packages are 
written in other languages, e.g. C modules, so this isn't an option.
Instead, we can run python, import the package and inspect it to get its 
members. I tried the `basedpyright --createstub` but it didn't generate 
anything except the basic high level wrapper for gemmi.
I also tried monkeytype which didn't work either.
The `stubgen` tool from `mypy` seems to be the best option. See `gemmi.sh`.
These stubs are used with the LSP, currently basedpyright, which knows about 
these stubs because the stubs path setting points here. See 
`lsp/basedpyright.lua`

