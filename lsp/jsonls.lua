return {
    cmd = { "vscode-json-language-server", "--stdio" },
    filetypes = { "json", "jsonc", "json.karabiner" },
    root_dir = require('utils/init').symlink_root_dir({ '.git' }),
    settings = {
        json = {
            schemas = {
                {
                    -- Complex modifications in assets folder (not main karabiner.json)
                    fileMatch = { "**/complex_modifications/*.json" },
                    url = "https://gitlab.com/notpushkin/karabiner-jsonschema/-/raw/master/schema.json"
                }
            },
            validate = { enable = true },
        }
    }
}
