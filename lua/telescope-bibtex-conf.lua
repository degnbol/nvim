#!/usr/bin/env lua
-- https://github.com/nvim-telescope/telescope-bibtex.nvim
telescope = require "telescope"
telescope.load_extension("bibtex")
-- the following is used to detect *.bib file used
-- by looking for \bibliography and \addbibresource.
telescope.setup {
      context = true,
}
