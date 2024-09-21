
shift+cmd+click in skim to go to section in neovim (skim must be started by neovim for this to work).

]] in insert mode to close env.

goto definition works for e.g. references to labels, e.g. with gd or ctrl+].


# NEOVIM HANGS
Sometimes nvim hangs after compilation.
Turns out it is hanging with VimtexView as well.
It is because skim is being messy.
Try to restart skim. See if the pdf was made.


# LSP
texlab is pretty ok at finding the bib file.
If a file e.g. tables/table.tex doesn't have any mention of it it can still get 
completion for bib as long as another file references tables/table.tex that itself is adding bib.
Worst case scenario (but this shouldn't be necessary if the file is \input'ed / \include'ed etc.)
use `\bibliography{../bibliography.bib}`.
It's an old deprecated version of `addbibresource` but it doesn't read the bib twice.

# tex

hyphen: -. Compound-word
en-dash: --. 1--4
em-dash: ---. Sentence --- inserted sentence --- rest of main sentence.
[1](https://www-users.york.ac.uk/~pjh503/LaTeX/misc.html#:~:text=Dashes%20and%20Hyphens,in%20text%20%2D%2D%2D%20like%20that!)
[2](https://tex.stackexchange.com/questions/3819/dashes-vs-vs)
