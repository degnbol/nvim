#!/usr/bin/env zsh
cd $0:h/..

if ! command -v aspell > /dev/null; then
    if [ `uname` = "Darwin" ]; then
        brew install aspell
    else
        echo "install aspell"
        exit 1
    fi
fi

mkdir -p spell/
cd spell/

aspell -d en dump master | aspell -l en expand > en.dic
aspell -d da dump master | aspell -l da expand > da.dic

for lang in da en; do
    echo $lang
    if [[ -f "$lang.utf-8.spl" ]]; then
        echo "Not overwriting: $lang.utf-8.spl"
    else
        nvim --clean --headless +"mkspell $PWD/$lang{,.dic} | q"
    fi
done

