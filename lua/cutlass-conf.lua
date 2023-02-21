#!/usr/bin/env lua
-- https://github.com/gbprod/cutlass.nvim
require("cutlass").setup {
    cut_key='x',
    override_del = true, -- true -> del key will put in blackhole register
}
