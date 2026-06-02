-- Insert-mode abbreviations (typos, contractions, short forms).
-- `abc` mirrors vim-abolish's :Abolish — brace expansion + case variants.
-- `ab` is a literal single iabbrev (no case variants).
-- `pab`/`pabc` are the same but register buffer-local; used in the prose-only
-- section below for abbreviations that would collide with code identifiers
-- (algo, dialog, prio, …) or English spellings that are valid in code (im, id).
-- See lua/utils/iabbrev.lua for the expansion algorithm.

local m = require("utils.iabbrev")
local function ab(lhs, rhs) m.iabbrev(lhs, rhs, false) end
local function abc(lhs, rhs) m.iabbrev(lhs, rhs) end

-- pab/pabc register buffer-local iabbrevs. An optional `pred` argument
-- delegates expansion to a runtime predicate (via `<expr>`), so the abbrev
-- only fires when the predicate returns true. Used for typst/tex where
-- prose is interleaved with math/code regions.
local function pab(lhs, rhs, pred) m.iabbrev(lhs, rhs, false, true, pred) end
local function pabc(lhs, rhs, pred) m.iabbrev(lhs, rhs, true, true, pred) end

-- TeX prose detection, via vimtex (lervag/vimtex).
-- Default is "not prose"; exceptions are explicit:
--   * Math (any form) is detected by in_mathzone — already not prose.
--   * Must be inside \begin{document} … \end{document} (rules out preamble).
--   * If inside a nested environment, that env must be in `tex_prose_envs`.
--   * If inside a \command{}, that command must be in `tex_prose_cmds`.
-- Starred forms (`section*`, `figure*`) are normalised before lookup.
local tex_prose_envs = {
    abstract = true,
    quote = true, quotation = true,
    itemize = true, enumerate = true, description = true,
    center = true,
    figure = true, table = true,  -- caption args carry the prose
    frame = true,                 -- beamer
}
local tex_prose_cmds = {
    -- sectioning
    section = true, subsection = true, subsubsection = true,
    chapter = true, part = true,
    paragraph = true, subparagraph = true,
    -- titles / metadata
    title = true, subtitle = true, author = true, date = true,
    -- captions / footnotes
    caption = true, footnote = true, footnotetext = true, marginpar = true,
    -- inline text formatting
    textbf = true, textit = true, emph = true,
    textrm = true, textsf = true, textsl = true, textsc = true,
    underline = true,
    -- quotations
    enquote = true,
}

local function tex_prose()
    if vim.fn["vimtex#syntax#in_mathzone"]() == 1 then return false end
    -- Preamble (outside \begin{document}) is never prose.
    local doc_pos = vim.fn["vimtex#env#is_inside"]("document")
    if not doc_pos or doc_pos[1] == 0 then return false end
    -- If nested in another env, it must be in the prose-env allowlist.
    local env = vim.fn["vimtex#env#get_inner"]()
    if type(env) == "table" and env.name and env.name ~= "" then
        local ename = (env.name):gsub("%*$", "")
        if not tex_prose_envs[ename] then return false end
    end
    -- If inside a command, its name must be in the prose-cmd allowlist.
    local cmd = vim.fn["vimtex#cmd#get_current"]()
    if type(cmd) == "table" and cmd.name and cmd.name ~= "" then
        local cname = (cmd.name):gsub("^\\", ""):gsub("%*$", "")
        return tex_prose_cmds[cname] == true
    end
    return true
end

-- Typst is bimodal: markup (prose) vs code. Math and raw are distinct modes.
-- A `content` block (`[…]`) inside code re-enters markup, so it counts as
-- prose. First mode-changer ancestor wins.
local function typst_prose()
    local ok, node = pcall(vim.treesitter.get_node)
    if not ok or not node then return true end
    while node do
        local t = node:type()
        if t == "math" or t == "formula"
            or t == "code"
            or t == "raw_blck" or t == "raw_span" then
            return false
        elseif t == "content" then
            return true
        end
        node = node:parent()
    end
    return true  -- top-level markup
end

-- ── Lazy apostrophe ──────────────────────────────────────────────────────
abc("{ca,is,are,do,does,did,has,have,had,was,were,would,should,could,wo}nt", "{}n't")
abc("{let,that,there,here,who,what}s", "{}'s")
abc("{they,you}re", "{}'re")
abc("yall", "y'all")
abc("itøs", "it's")
ab("THeres", "There's")
-- I-pronoun contractions: `im`, `ive`, `id` collide with code (julia vars,
-- "five" suffix, HTML/DB id) so they live in the prose-only section below.
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
-- ;-for-' (adjacent keys on QWERTY). Generic suffix abbrevs: don;t → don't,
-- it;s → it's, I;m → I'm, etc. `;X` is end-id (non-kwd then kwd), so it's a
-- valid iabbrev lhs; `it;s` itself isn't. Limited to single-letter suffixes —
-- `;ll`/`;ve`/`;re` can't be expressed (kwd before final kwd).
abc(";s", "'s")
abc(";t", "'t")
abc(";m", "'m")
abc(";d", "'d")

-- ── Shortcuts ────────────────────────────────────────────────────────────
abc("aa{,s}", "amino acid{,s}")
abc("paren{,s}", "parenthes{i,e}s")
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
abc("melb", "Melbourne")
abc("noone", "no one")
abc("ie", "i.e.")  -- i..e isn't valid as a keyword; snippet in luasnippets/all.lua handles the dotted form
abc("wildtype", "wild type")

-- ── Auto-capitalization ──────────────────────────────────────────────────
ab("english", "English")
ab("danish", "Danish")

-- ── Special characters ───────────────────────────────────────────────────
ab("oC", "°C")

-- ── Prose only ───────────────────────────────────────────────────────────
-- Registered buffer-local for prose filetypes via FileType autocmd. For
-- tex/typst (mixed prose+math+code), a treesitter predicate gates expansion
-- so abbrevs only fire in prose regions.
local function setup_prose_abbrevs(pred)
    pab("im", "I'm", pred)
    pab("ive", "I've", pred)
    pab("id", "I'd", pred)
    pabc("intro", "introduction", pred)
    pabc("algo{,s}", "algorithm{}", pred)
    pabc("dialog", "dialogue", pred)
    pabc("avail", "available", pred)
    pabc("bc", "because", pred)
    pabc("prio", "priority", pred)
    pabc("prios", "priorities", pred)
end

-- Explicit prose filetypes. Most don't parse as markdown (tex, rst, …) so can't
-- be inferred from the treesitter language. `quarto` is here too: nothing
-- registers it to the markdown parser in this config, so get_lang("quarto")
-- stays "quarto".
local prose_fts = {
    markdown = true, asciidoc = true, tex = true, typst = true, text = true,
    gitcommit = true, mail = true, rst = true, org = true, neorg = true,
    quarto = true,
}
local ft_predicate = { tex = tex_prose, typst = typst_prose }

--- A filetype is prose if it's in the allowlist or parses as markdown. The
--- get_lang fallback catches custom filetypes registered to the markdown parser
--- (e.g. agentic.nvim's `AgenticInput`) with no per-filetype maintenance here.
--- @param ft string single (non-compound) filetype
--- @return boolean
local function is_prose_ft(ft)
    return prose_fts[ft] or vim.treesitter.language.get_lang(ft) == "markdown"
end

-- Decide prose-ness per dot-component rather than by exact name. Compound
-- filetypes that inherit from markdown (`markdown.mdx`, `vimwiki.markdown`, …)
-- load markdown's syntax and ftplugin via neovim's dotted runtime loading, but a
-- bare-name autocmd pattern never matches them, because pattern matching treats
-- the filetype as one opaque string. A buffer is prose when any of its
-- components is. The math-gating predicate comes from whichever component
-- defines one (tex/typst).
vim.api.nvim_create_autocmd("FileType", {
    group = vim.api.nvim_create_augroup("ProseAbbrevs", { clear = true }),
    callback = function(ev)
        local is_prose, pred = false, nil
        for _, part in ipairs(vim.split(ev.match, ".", { plain = true })) do
            if is_prose_ft(part) then is_prose = true end
            if ft_predicate[part] then pred = ft_predicate[part] end
        end
        if is_prose then setup_prose_abbrevs(pred) end
    end,
})

-- TODO: make this a bit more convenient, for autocorrecting language-specific
-- spelling errors in prose only (not in code/comments).
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
