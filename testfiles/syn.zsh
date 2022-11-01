#!/usr/bin/env zsh
for file in */*.tsv; do
    echo $file:r:t
    echo ${file:r:h}
    echo "$file:t"
done

# comment with apostrophe '
if [ -f "filename.pdf" ]; then
    ./script.jl -a 3 --arg davs
elif [[ -n 4 ]]; then
    fisk subcommand ok hej
fi

git clone --recurse-submodules git@github.com:degnbol/dotfiles.git --cwd=$HOME
kitty +kitten ssh username@server.dk
kitty +kitten ssh username@server

sed 's/a /b_/g' file.txt
sed "s/a_/b /g" file.txt
time brew install lars
sudo ls ~/hejsa.txt
ls *nice.pdf
# comment
echo "string"
ROOT=`git root`
mkdir $hejsa/fisk
cd .; cd -; cd ~; cd ~/nvim
sudo mkdir a/
source file.sh
. ~/dir/file.sh
~/dir/file.sh

echo "\ "
$0:h/pdb2tsv.sh | mlr --t2p filter '$atom == "CA"' +\
    cut -f x,y,z,resi,chain then uniq -a + count


