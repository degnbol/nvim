return {
	settings = {
		ltex_plus = {
			language = "en-GB",
			disabledRules = {
				["en-GB"] = {
					"COMMA_PARENTHESIS_WHITESPACE",
					"MORFOLOGIK_RULE_EN_GB", -- misspellings
				},
			},
			dictionary = {
				["en-GB"] = { ":" .. vim.fn.stdpath("config") .. "/spell/custom.utf8.add" },
			},
		},
	},
}
