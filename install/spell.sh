if ! command -v aspell > /dev/null; then
    if [ `uname` = "Darwin" ]; then
        brew install aspell
    else
        echo "requires aspell"
        exit 1
    fi
fi

cd $XDG_CONFIG_HOME/nvim
mkdir -p spell
aspell -d en dump master | aspell -l en expand > spell/en.dic
aspell -d da dump master | aspell -l da expand > spell/da.dic
# if `set spelllang=da` or `set spelllang=en,da` then install danish understood by vim spell with
nvim +'set spell|q'
# then follow prompt and install danish utf8

