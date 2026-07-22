vim.filetype.add {
    extension = {
        osascript = "applescript",
        aql = "aql",
        usage = "complgen",
        crd = "asciidoc.chordpro",
        cypher = "cypher", cyp = "cypher", cql = "cypher",
        mlr = "miller",
        pdb = "pdb",
        piano = "piano",
        pml = "pymol",
        cif = "star",
        info = "text",
        tsv = "tsv", tab = "tsv", bed = "tsv",
        ipynb = "python",
        gitconfig = "gitconfig",
        ndjson = "jsonl", -- jsonl is native; ndjson is the same format
    },
    filename = {
        ["karabiner.json"] = "json.karabiner",
    },
    pattern = {
        [".*%.blend%.py"] = "python.blender",
        [".*%.py"] = {
            function(_, bufnr)
                local first = vim.api.nvim_buf_get_lines(bufnr, 0, 1, false)[1] or ""
                if first:match("blender") then
                    return "python.blender"
                end
            end,
            { priority = 10 },
        },
    },
}
