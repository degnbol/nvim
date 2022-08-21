#!/usr/bin/env bash
rm -f "$XDG_CONFIG_HOME/nvim/plugin/packer_compiled.lua"
nvim +'PackerSync'
