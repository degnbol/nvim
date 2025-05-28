#!/usr/bin/env lua
return {
    settings = { yaml = { schemas = {
        [vim.opt.runtimepath:get()[1] .. "/lua/completion/asciidoc/asciidoc-theme.json"] = "*.adoc.yml",
    }}}
}
