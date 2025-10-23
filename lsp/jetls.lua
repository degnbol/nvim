-- See install/jetls.sh
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
