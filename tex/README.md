A pandoc filter is used as a formatter to convert between normal latex math 
like \alpha and unicode math like α. The latter should also be made into the 
exact same pdf output using the unicode-math package and the right latex 
flavour, e.g. xelatex or lualatex. It is mostly intended as a fail-safe, in 
case I don't want to use unicode in my writing, for whichever reason, e.g. 
collaborating with others, something not working, etc.

Things to note:
-   It discards preamble so should only be run on text inside the document 
    environment.
-   Some subtle maybe changes, e.g. introducing a hypertarget to figure text, 
    should only be used within math environments.
-   Converts inline math in $$ to \(\), however, this doesn't change 
    the output.
-   It compacts some lines, which again shouldn't change the output.
