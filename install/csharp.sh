if ! command -v dotnet > /dev/null; then
    if [ `uname` = "Darwin" ]; then
        brew install dotnet
    else
        echo "Requires dotnet"
        exit 1
    fi
fi
# preferred install method according to
# https://github.com/neovim/nvim-lspconfig/blob/master/doc/server_configurations.md#csharp_ls
dotnet tool install --global csharp-ls
# otherwise do something like
# nvim --headless -c "MasonInstall csharp-language-server" +q
# instead of using ensure_installed

