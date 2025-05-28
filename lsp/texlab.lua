#!/usr/bin/env lua
-- Don't autoshow completion in cmdline
-- https://github.com/neovim/nvim-lspconfig/blob/master/doc/configs.md#texlab
return {
    -- https://github.com/latex-lsp/texlab/wiki/Configuration
    settings = { texlab = {
        experimental = {
            -- custom citation function completion.
            -- In some projects I use \citea{...} as \citeauthor{...}~\cite{...} and \citePDBlit etc.
            citationCommands = {"citea", "citePDB", "citePDBlit"},
            -- \subref from the subpcation package.
            -- \footref from the footmisc package to refer to footnote a second (or more) time.
            -- In some projects I use \see[...]{...} as (see~\cref{...}...)
            labelReferenceCommands = {"footref", "subref", "see", "seename", "seefull"},
        },
    }}
}
