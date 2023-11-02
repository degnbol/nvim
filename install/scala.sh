#!/usr/bin/env zsh

if [ `uname` = "Darwin" ]; then
    brew install coursier/formulas/coursier
    cs setup
    echo 'export PATH="$PATH:$HOME/Library/Application\ Support/Coursier/bin"' >> ~/.zshrc
else
    curl -fL "https://github.com/VirtusLab/coursier-m1/releases/latest/download/cs-aarch64-pc-linux.gz" | gzip -d > cs
    chmod +x ./cs
    ./cs setup
    echo 'export PATH="$PATH:$HOME/.local/share/coursier/bin"' >> ~/.zshrc
    rm ./cs
fi

