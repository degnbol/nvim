return {
--

-- easily change to zsh
s({
    trig="#!",
    dscr="zsh shebang",
    snippetType='autosnippet',
    condition=conds.line_begin,
},
{t{"#!/usr/bin/env zsh", ""}}),

}
