
-- NOTE: Icons should be consistent with ~/dotfiles/exa/icons.sh
-- Colors:
-- Binary data files. color = "#89e051", cterm_color = "113",
-- Text files. color = "#519aba", cterm_color = "67",
return {
    {
        "nvim-tree/nvim-web-devicons",
        opts = {
            override = {
                blend = {
                    icon = "󰂫",
                    name = "Blender",
                },
                blend1 = {
                    icon = "󰂫",
                    name = "BlenderBackup",
                },
                bed = {
                    icon = "",
                    color = "#89e051",
                    cterm_color = "113",
                    name = "BrowserExtensibleData"
                },
                -- .bg2 is bed graph 2D which is just a tsv file with columns chrom1,start1,end1,chrom2,start2,end2,count
                bg2 = {
                    icon = "",
                    color = "#89e051",
                    cterm_color = "113",
                    name = "BEDGraph2D"
                },
                csv = {
                    icon = "",
                    color = "#89e051",
                    cterm_color = "113",
                    name = "CommaSeparatedValues"
                },
                tab = {
                    icon = "",
                    color = "#89e051",
                    cterm_color = "113",
                    name = "Table"
                },
                tsv = {
                    icon = "",
                    color = "#89e051",
                    cterm_color = "113",
                    name = "TabSeparatedValues"
                },
                ssv = {
                    icon = "",
                    color = "#89e051",
                    cterm_color = "113",
                    name = "SpaceSeparatedValues"
                },
                bw = {
                    icon = "󰚄",
                    color = "#89e051",
                    cterm_color = "113",
                    name = "BigWig",
                },
                fa = {
                    icon = "󰚄",
                    color = "#89e051",
                    cterm_color = "113",
                    name = "Fasta",
                },
                faa = {
                    icon = "󰚄",
                    color = "#89e051",
                    cterm_color = "113",
                    name = "FastaAminoAcids",
                },
                fasta = {
                    icon = "󰚄",
                    color = "#89e051",
                    cterm_color = "113",
                    name = "Fasta",
                },
                -- an HDF5 file with specific groups and datasets for storing genomic data, specifically Hi-C
                cool = {
                    icon = "",
                    color = "#89e051",
                    cterm_color = "113",
                    name = "Cooler",
                },
                h5 = {
                    icon = "",
                    color = "#89e051",
                    cterm_color = "113",
                    name = "HDF5",
                },
                gz = {
                    icon = "",
                    name = "GZip",
                },
                mat = {
                    icon = "󰘨",
                    color = "#89e051",
                    cterm_color = "113",
                    name = "Matrix"
                },
                npy = {
                    icon = "󰘨",
                    color = "#89e051",
                    cterm_color = "113",
                    name = "Numpy"
                },
                pdf = {
                    icon = "",
                    color = "#b30b00",
                    ctermfg = 124,
                    name = "PortableDocumentFormat",
                },
                pkl = {
                    icon = "",
                    color = "#89e051",
                    cterm_color = "113",
                    name = "Pickle",
                },
                adoc = {
                    icon = "",
                    color = "#519aba",
                    cterm_color = "67",
                    name = "Asciidoc",
                },
                crd = {
                    icon = "󰝱",
                    color = "#519aba",
                    cterm_color = "67",
                    name = "Chordpro",
                },
                md = {
                    icon = "",
                    color = "#519aba",
                    cterm_color = "67",
                    name = "Markdown",
                },
                -- NOTE: We need to set both md and markdown.
                markdown = {
                    icon = "",
                    color = "#519aba",
                    cterm_color = "67",
                    name = "Markdown",
                },
                tex = {
                    icon = "",
                    color = "#519aba",
                    cterm_color = "67",
                    name = "LaTeX",
                },
                typ = {
                    icon = "",
                    color = "#519aba",
                    cterm_color = "67",
                    name = "Typst",
                },
                wgsl = {
                    icon = "󰢮",
                    color = "#5586A6",
                    cterm_color = "67",
                    name = "WGSL",
                },
            }
        }
    }
}
