#!/usr/bin/env zsh

# install pandoc
if ! command -v pandoc > /dev/null; then
    if command -v brew > /dev/null; then
        brew install pandoc
    else
        mamba install pandoc
    fi
fi

if [ `uname` = "Darwin" ]; then
    os="macOS"
else
    os="Linux"
fi

cd $0:h
# install latest converter filter
# https://github.com/marhop/pandoc-unicode-math/releases/latest
release="v3.1.0"
pandocTypes="1.23"
wget https://github.com/marhop/pandoc-unicode-math/releases/download/$release/pandoc-unicode-math_${os}_pandoc-types-$pandocTypes.zip
unzip pandoc-unicode-math*.zip && rm pandoc-unicode-math*.zip
chmod +x pandoc-unicode-math*

