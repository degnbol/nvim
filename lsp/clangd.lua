return {
    -- Used only for files without a compile command (no compile_commands.json or compile_flags.txt).
    -- Required here, not at file load: vim.lsp.enable loads this file on every startup,
    -- and the flags come from pkg-config calls.
    before_init = function(params)
        params.initializationOptions = vim.tbl_extend("force", params.initializationOptions or {}, {
            fallbackFlags = require "c_fallback_flags",
        })
    end,
}
