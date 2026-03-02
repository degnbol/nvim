
-- :MasonInstall r-languageserver
-- Patch: resolve ... formals in wrapper functions (e.g. scale_y_log10 → scale_y_continuous)
local patch = vim.fn.stdpath("config") .. "/lsp_ext/r_lsp_dots.R"
return {
    cmd = { "R", "--no-echo", "-e", "source('" .. patch .. "'); languageserver::run()" },
}
