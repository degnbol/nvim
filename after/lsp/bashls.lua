
return {
    filetypes = { "sh", "bash", "zsh", "sh.zsh" },
    -- SC1071 ("only supports sh/bash/…") is shellcheck's only output for a zsh shebang.
    settings = { bashIde = { shellcheckArguments = "-e SC1071" } },
}
