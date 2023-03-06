-- show "|" on indented lines
return {"lukas-reineke/indent-blankline.nvim", enabled=false, config=function ()
    local g = vim.g

    g.indent_blankline_char = "▏"
    g.indent_blankline_filetype_exclude = {"help", "terminal", "dashboard", "tsv"}
    g.indent_blankline_buftype_exclude = {"terminal"}
    g.indent_blankline_show_trailing_blankline_indent = false
    g.indent_blankline_show_first_indent_level = false

    g.indentLine_concealcursor = 'inc'
    g.indentLine_conceallevel = 2
    g.conceallevel = 2
end}
