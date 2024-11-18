#!/usr/bin/env zsh
cd $0:h
mkdir -p ./pymol_settings_descriptions/

cat ./pymol_settings.txt | while read line; do
    echo $line
    wget -O temp.html "https://pymolwiki.org/index.php/$line"
    grep -A2 'id="Overview"' temp.html | grep -A2 '<p>' | grep -B2 '<\p>' > ./pymol_settings_descriptions/$line.txt
done
