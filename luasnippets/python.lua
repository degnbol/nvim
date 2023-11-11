return {
--

s({trig="#!", dscr="Python shebang", snippetType="autosnippet", condition=conds.line_begin},
{t{"#!/usr/bin/env python3", ""}}),

s({
    trig='""""',
    dscr="Tripple quotes",
    snippetType='autosnippet',
},
fmta([["""
<>
"""
]], i(1))),

-- TODO: make these work better, i.e. show up if typing numpy at line start or 
-- in case one types import n then it shouldn't write import twice
s({ trig="import numpy as np", dscr="Import numpy.", wordTrig=false, },
{t{"import numpy as np", ""}}),
s({trig="import pandas as pd", dscr="Import pandas.", wordTrig=false,},
{t{"import pandas as pd", ""}}),
s({trig="import scipy as sp", dscr="Import scipy.", wordTrig=false,},
{t{"import scipy as sp", ""}}),
s({trig="import matplotlib.pyplot as plt", dscr="Import pyplot.", wordTrig=false,},
{t{"import matplotlib.pyplot as plt", ""}}),

s({trig="DF", dscr="Pandas DataFrame.", snippetType="autosnippet"},
{t"pd.DataFrame(", i(1), t")"}),

}
