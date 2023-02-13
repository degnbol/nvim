#!/usr/bin/env bash
rm -f "$XDG_CONFIG_HOME/nvim/plugin/packer_compiled.lua"
nvim +'PackerSync'
# if it hangs forever you will need to call PackerCompile explicitly.
