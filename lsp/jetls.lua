-- See install/jetls.sh
-- This or julials is currently chosen in ftplugin/julia.lua
local JETLS = vim.fs.abspath("~/.local/share/JETLS.jl")
return {
    cmd = {
        "julia",
        "--startup-file=no",
        "--history-file=no",
        "--project=" .. JETLS,
        "--threads=auto",
        JETLS .. "/runserver.jl",
    },
    filetypes = {"julia"},
}
