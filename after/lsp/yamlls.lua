-- Plan files of the pi `plans` extension, in repo stores and the global store.
local plan_schema = vim.fs.dirname(vim.fn.stdpath("config")) .. "/pi/extensions/plans/schema.json"
local data_home = vim.env.XDG_DATA_HOME or (vim.env.HOME .. "/.local/share")

return {
	settings = {
		yaml = {
			schemas = {
				[vim.fn.stdpath("config") .. "/lua/completion/asciidoc/asciidoc-theme.json"] = "*.adoc.yml",
				[plan_schema] = { "**/.agents/plans/*.yaml", data_home .. "/agents/plans/*.yaml" },
			},
		},
	},
}
