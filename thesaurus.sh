wget https://github.com/vim/vim/files/2634525/thesaurus_pkg.zip
unzip thesaurus_pkg.zip
mkdir -p thesaurus
mv thesaurus_pkg/thesaurus.txt thesaurus/english.txt
rm -rf thesaurus_pkg{,.zip}
