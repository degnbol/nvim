#!/usr/bin/env zsh
cd $0:h

# install pandoc
if ! type pandoc &> /dev/null; then
    if type brew &> /dev/null; then
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

# install latest converter filter
# https://github.com/marhop/pandoc-unicode-math/releases/latest
release="v3.1.0"
pandocTypes="1.23"
wget https://github.com/marhop/pandoc-unicode-math/releases/download/$release/pandoc-unicode-math_${os}_pandoc-types-$pandocTypes.zip
unzip pandoc-unicode-math*.zip && rm pandoc-unicode-math*.zip
chmod +x pandoc-unicode-math*

