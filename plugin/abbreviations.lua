-- Insert-mode abbreviations (typos, contractions, short forms).
-- `abc` mirrors vim-abolish's :Abolish — brace expansion + case variants.
-- `ab` is a literal single iabbrev (no case variants).
-- See lua/utils/iabbrev.lua for the expansion algorithm.

local m = require("utils.iabbrev")
local function ab(lhs, rhs) m.iabbrev(lhs, rhs, false) end
local function abc(lhs, rhs) m.iabbrev(lhs, rhs) end

-- ── Lazy apostrophe ──────────────────────────────────────────────────────
abc("{ca,is,are,do,does,did,has,have,had,was,were,would,should,could,wo}nt", "{}n't")
abc("{let,that,there,here,who,what}s", "{}'s")
abc("{they,you}re", "{}'re")
abc("yall", "y'all")
abc("itøs", "it's")
ab("THeres", "There's")
-- I-pronoun contractions: lowercase `i` is too common in code/comments to
-- expand globally — only the explicitly-cased forms below.
-- ab("im", "I'm")  -- triggered in e.g. julia code
-- ab("ive", "I've")  -- annoys when writing "five"
-- ab("id", "I'd")  -- id is a word on its own
ab("Im", "I'm")
ab("IM", "I'm")
ab("Ive", "I've")
ab("Ill", "I'll")
ab("Id", "I'd")
ab("youll", "you'll")

-- ── Typos ────────────────────────────────────────────────────────────────
-- extra capitalization (Shift held too long)
abc("TH{ere,en,e,is}", "Th{}")
-- flipped letters
ab("fucntion", "function")
abc("liek", "like")
abc("ahve", "have")
abc("waht", "what")
abc("sohw", "show")
abc("blaance", "balance")
abc("sohuld", "should")
abc("tihnk", "think")
abc("shoudl", "should")
abc("udnerstand", "understand")
abc("palce", "place")
abc("unqiue", "unique")
abc("simplicies", "simplices")
abc("{despa,sepe}rat{e,es,ed,ing,ely,ion,ions,or}", "{despe,sepa}rat{}")
-- misspellings
abc("flourescent{,ly}", "fluorescent{}")
abc("eucledian", "Euclidean")
abc("lifes", "lives")
abc("pertruding", "protruding")
abc("effecient", "efficient")
abc("persuit", "pursuit")
abc("oc{,c}uring", "occurring")
abc("feasab{ility,le}", "feasib{ility,le}")
abc("preceed{,ed,ing}", "preced{,ed,ing}")
abc("embarass{,ing,ingly}", "embarrass{,ing,ingly}")
abc("corespond{,s,ing}", "correspond{,s,ing}")
abc("discernable", "discernible")
abc("rudamentary", "rudimentary")
abc("occurence", "occurrence")
abc("occur{,r}ance", "occurrence")
abc("hi{,e}ra{,r}ch{y,ical}", "hi{e}{r}arch{}")

-- ── Shortcuts ────────────────────────────────────────────────────────────
abc("aa{,s}", "amino acid{,s}")
abc("paren{,s}", "parenthes{i,e}s")
abc("algo{,s}", "algorithm{}")
abc("combo{,s}", "combination{}")
abc("tho", "though")
abc("altho", "although")
abc("eventho", "even though")
abc("inspite", "in spite")
abc("defacto", "de facto")
abc("thru", "through")
abc("passthru", "passthrough")
abc("probs", "probably")
abc("ppl", "people")
abc("dialog", "dialogue")
abc("avail", "available")
abc("bc", "because")
abc("melb", "Melbourne")
abc("prio", "priority")
abc("prios", "priorities")
abc("noone", "no one")
abc("ie", "i.e.")  -- i..e isn't valid as a keyword; snippet in luasnippets/all.lua handles the dotted form
abc("wildtype", "wild type")
-- TODO: only do some abbrevs in regular text
-- abc("intro", "introduction")

-- ── Auto-capitalization ──────────────────────────────────────────────────
ab("english", "English")
ab("danish", "Danish")

-- ── Special characters ───────────────────────────────────────────────────
ab("oC", "°C")

-- Language-specific abbreviations toggled by `iminsert` (Danish vs English).
-- See lua/utils/keymap for the iminsert toggle binding.
local function toggle_dansk_abbrev()
    local buf = vim.api.nvim_get_current_buf()
    local function silent_cmd(cmd) pcall(function() vim.cmd(cmd) end) end
    if vim.bo.iminsert ~= 0 then
        vim.cmd("abbrev feks f.eks.")
        silent_cmd("unabbrev eg")
        silent_cmd("unabbrev Eg")
        silent_cmd("iunabbrev <buffer> ti")
        silent_cmd("iunabbrev <buffer> i")
    else
        silent_cmd("iunabbrev <buffer> feks")
        vim.cmd("iabbrev eg e.g.")
        vim.cmd("iabbrev Eg E.g.")
        if vim.bo[buf].filetype == "asciidoc" then
            -- For regular text where we wouldn't be talking about a variable i
            -- or in Danish where i is a word.
            vim.cmd("iabbrev <buffer> i I")
            vim.cmd("iabbrev <buffer> ti it")
        end
    end
end
_G.ToggleDanskAbbrev = toggle_dansk_abbrev  -- kept global for :call from autocmds elsewhere

toggle_dansk_abbrev()
local grp = vim.api.nvim_create_augroup("ToggleDanskAbbrev", { clear = true })
vim.api.nvim_create_autocmd("User", {
    group = grp,
    pattern = "ToggleDansk",
    callback = toggle_dansk_abbrev,
})
vim.api.nvim_create_autocmd("FileType", {
    group = grp,
    callback = toggle_dansk_abbrev,
})
