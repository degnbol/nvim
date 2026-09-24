return {
	settings = {
		yaml = {
			schemas = {
				[vim.fn.stdpath("config") .. "/lua/completion/asciidoc/asciidoc-theme.json"] = "*.adoc.yml",
			},
		},
	},
}
