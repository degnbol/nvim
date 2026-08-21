-- A plugin rather than an ftplugin: a ```smiles fence in a markdown buffer takes
-- the same @chem.* groups but has no filetype of its own to hang them on.
require("chem/highlight").setup()
require("chem/tsv").setup()
