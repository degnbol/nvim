-- Also sourced for cpp by the runtime's ftplugin/cpp.vim.
local map = require "utils/keymap"

map.n("<localleader>r", function() require("compile_db").regenerate(0) end,
    "Regenerate compile_commands.json and restart LSP", { buffer = true })
