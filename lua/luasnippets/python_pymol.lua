#!/usr/bin/env lua
local ls = require("luasnip")
local s = ls.snippet
local sn = ls.snippet_node
local isn = ls.indent_snippet_node
local t = ls.text_node
local i = ls.insert_node
local f = ls.function_node
local c = ls.choice_node
local d = ls.dynamic_node
local r = ls.restore_node
local events = require("luasnip.util.events")
local ai = require("luasnip.nodes.absolute_indexer")
local extras = require("luasnip.extras")
local l = extras.lambda
local rep = extras.rep
local p = extras.partial
local m = extras.match
local n = extras.nonempty
local dl = extras.dynamic_lambda
local fmt = require("luasnip.extras.fmt").fmt
local fmta = require("luasnip.extras.fmt").fmta
local conds = require("luasnip.extras.expand_conditions")
local postfix = require("luasnip.extras.postfix").postfix
local types = require("luasnip.util.types")
local parse = require("luasnip.util.parser").parse_snippet
local ms = ls.multi_snippet
local k = require("luasnip.nodes.key_indexer").new_key

local lsu = require "utils/luasnip"
local re = lsu.re

local function in_string()
    -- very simple way of checking if we are in string with treesitter
    local captures = vim.treesitter.get_captures_at_cursor()
    return vim.tbl_contains(captures, "string")
end
local cond = {show_condition=in_string}

return {
    ms({
        {trig="[^%w]/", trigEngine="pattern"},
        {trig="///", snippetType="autosnippet"},
        common={dscr="Selection macro https://pymolwiki.org/index.php/Selection_Macros"},
    },
    fmta([[/<>/<>/<>/<>/<>]], {
        i(1, "OBJ"),
        i(2, "SEGI"),
        i(3, "CHAIN"),
        i(4, "1-4"),
        i(5, "CA"),
    }), cond),

    s({trig="all", dscr="All atoms currently loaded into PyMOL"}, {t"all"}, cond),
    s({trig="*", dscr="All atoms currently loaded into PyMOL"}, {t"*"}, cond),
    s({trig="none", dscr="Empty selection"}, {t"none"}, cond),
    s({trig="enabled", dscr="Atoms from enabled objects"}, {t"enabled"}, cond),
    s({trig="not", dscr="Inverts selection"}, {t"not ", i(1, "S1")}, cond),
    s({trig="!", dscr="Inverts selection"}, {t"! ", i(1, "S1")}, cond),
    s({trig="and", dscr="Atoms included in both S1 and S2"}, {t"and ", i(1,"S2")}, cond),
    s({trig="&", dscr="Atoms included in both S1 and S2"}, {t"& ", i(1,"S2")}, cond),
    s({trig="or", dscr="Atoms included in either S1 and S2"}, {t"or ", i(1,"S2")}, cond),
    s({trig="|", dscr="Atoms included in either S1 and S2"}, {t"| ", i(1,"S2")}, cond),
    s({trig="first", dscr="First atom in S1 (single atom selection)"}, {t"first ", i(1,"S1")}, cond),
    s({trig="last", dscr="Last atom in S1 (single atom selection)"}, {t"last ", i(1,"S1")}, cond),
    s({trig="model", dscr='Atoms from object'}, {t"model ", i(1,"OBJ")}, cond),
    s({trig="m.", dscr='Atoms from object'}, {t"m. ", i(1,"OBJ")}, cond),
    s({trig="chain", dscr='Atoms from chain'}, {t"chain ", i(1,"CHAIN")}, cond),
    s({trig="c.", dscr='Atoms from chain'}, {t"c. ", i(1,"CHAIN")}, cond),
    s({trig="segi", dscr='Atoms from segment'}, {t"segi ", i(1,"SEG")}, cond),
    s({trig="s.",   dscr='Atoms from segment ("label_asym_id" in mmCIF)'}, {t"s. ", i(1,"SEG")}, cond),
    s({trig="resn",   dscr='Residue name'}, {t"resn ", i(1,"ALA")}, cond),
    s({trig="r.",   dscr='Residue name'}, {t"r. ", i(1,"ALA")}, cond),
    s({trig="resi",   dscr='Residue number or range'}, {t"resi ", i(1,"100-200")}, cond),
    s({trig="i.",   dscr='Residue number or range'}, {t"i. ", i(1,"100-200")}, cond),
    s({trig="name",   dscr='Atom name'}, {t"name ", i(1,"CA")}, cond),
    s({trig="n.",   dscr='Atom name'}, {t"n. ", i(1,"CA")}, cond),
    s({trig="alt",   dscr='Alternate location'}, {t"alt ", i(1,"A")}, cond),
    s({trig="index", dscr='Internal per-object atom index (changes with sorting)'}, {t"index ", i(1,"123")}, cond),
    s({trig="idx.", dscr='Internal per-object atom index (changes with sorting)'}, {t"idx. ", i(1,"123")}, cond),
    s({trig="id", dscr='ID column from PDB file'}, {t"id ", i(1,"123")}, cond),
    s({trig="rank", dscr='Per-object atom index at load time (see also retain_order)'}, {t"rank ", i(1,"123")}, cond),
    s({trig="pepseq", dscr='Protein residue sequence with given one-letter (see also FindSeq)'}, {t"pepseq ", i(1,"ACDEF")}, cond),
    s({trig="ps.", dscr='Protein residue sequence with given one-letter (see also FindSeq)'}, {t"ps. ", i(1,"ACDEF")}, cond),
    s({trig="label", dscr='Atoms with given label'}, {t"label ", i(1,"LABEL")}, cond),
    s({trig="in", dscr='Atoms in S1 whose identifiers name, resi, resn, chain and segi all match atoms in S2'}, {t"in ", i(1,"S2")}, cond),
    s({trig="like", dscr='Atoms in S1 whose identifiers name and resi match atoms in S2'}, {t"like ", i(1,"S2")}, cond),
    s({trig="byobject", desc='Expands S1 to complete objects'}, {t"byobject"}, cond),
    s({trig="bysegi", desc='Expands S1 to complete segments'}, {t"bysegi ", i(1,"S1")}, cond),
    s({trig="bs.", desc='Expands S1 to complete segments'}, {t"bs. ", i(1,"S1")}, cond),
    s({trig="bychain", desc='Expands S1 to complete chains'}, {t"bychain ", i(1,"S1")}, cond),
    s({trig="bc.", desc='Expands S1 to complete chains'}, {t"bc. ", i(1,"S1")}, cond),
    s({trig="byres", desc='Expands S1 to complete residues'}, {t"byres ", i(1,"S1")}, cond),
    s({trig="br.", desc='Expands S1 to complete residues'}, {t"byres ", i(1,"S1")}, cond),
    s({trig="bycalpha", desc='A atoms of residues with at least one atom in S1'}, {t"bycalpha ", i(1,"S1")}, cond),
    s({trig="bca.", desc='A atoms of residues with at least one atom in S1'}, {t"bca. ", i(1,"S1")}, cond),
    s({trig="bymolecule", desc='Expands S1 to complete molecules (connected with bonds)'}, {t"bymolecule ", i(1,"S1")}, cond),
    s({trig="bm.", desc='Expands S1 to complete molecules (connected with bonds)'}, {t"bm. ", i(1,"S1")}, cond),
    s({trig="byfragment", desc=''}, {t"byfragment ", i(1,"S1")}, cond),
    s({trig="bf.", desc=''}, {t"bf. ", i(1,"S1")}, cond),
    s({trig="byring", desc='All rings of size ≤ 7 which have at least one atom in S1'}, {t"byring ", i(1,"S1")}, cond),
    s({trig="bycell", desc='Expands selection to unit cell'}, {t"bycell ", i(1,"S1")}, cond),
    s({trig="bound_to", desc='Atoms directly bonded to S1, may include S1'}, {t"bound_to ", i(1,"S1")}, cond),
    s({trig="bto.", desc='Atoms directly bonded to S1, may include S1'}, {t"bto. ", i(1,"S1")}, cond),
    s({trig="neighbor", desc='Atoms directly bonded to S1, excludes S1'}, {t"neighbor ", i(1,"S1")}, cond),
    s({trig="nbr.", desc='Atoms directly bonded to S1, excludes S1'}, {t"nbr. ", i(1,"S1")}, cond),
    s({trig="extend", desc='Expands S1 by N bonds connected to atoms in S1'}, {i(1,"S1"), t" extend ", i(2,"3")}, cond),
    s({trig="xt.", desc='Expands S1 by N bonds connected to atoms in S1'}, {i(1,"S1"), t" xt. ", i(2,"3")}, cond),
    s({trig="within", desc='Atoms in S1 that are within N Angstroms of any atom in S2'}, {i(1,"S1"), t" within ", i(2,"12.3"), t" of ", i(3,"S2")}, cond),
    s({trig="w.", desc='Atoms in S1 that are within N Angstroms of any atom in S2'}, {i(1,"S1"), t" w. ", i(2,"12.3"), t" of ", i(3,"S2")}, cond),
    s({trig="around", desc='Atoms with centers within N Angstroms of the center of any atom in S1'}, {i(1,"S1"), t" around ", i(2,"12.3")}, cond),
    s({trig="a.", desc='Atoms with centers within N Angstroms of the center of any atom in S1'}, {i(1,"S1"), t" a. ", i(2,"12.3")}, cond),
    s({trig="expand", desc='Expands S1 by atoms within 12.3 Angstroms of the center of any atom in S1'}, {i(1,"S1"), t" expand ", i(2,"12.3")}, cond),
    s({trig="x.", desc='Expands S1 by atoms within 12.3 Angstroms of the center of any atom in S1'}, {i(1,"S1"), t" x. ", i(2,"12.3")}, cond),
    s({trig="gap", desc='Atoms whose VDW radii are separated from the VDW radii of S1 by a minimum of N Angstroms.'}, {i(1,"S1"), t" gap ", i(2,"1.2")}, cond),
    s({trig="near_to", desc='Same as within, but excludes S2 from the selection (and thus is identical to S1 and S2 around 12.3)'}, {i(1,"S1"), t" near_to ", i(2,"12.3"), t" of ", i(3,"S2")}, cond),
    s({trig="nto.", desc='Same as within, but excludes S2 from the selection (and thus is identical to S1 and S2 around 12.3)'}, {i(1,"S1"), t" nto. ", i(2,"12.3"), t" of ", i(3,"S2")}, cond),
    s({trig="beyond", desc='Atoms in S1 that are at least 12.3 Anstroms away from S2'}, {i(1,"S1"), t" beyond ", i(2,"12.3"), t" of ", i(3,"S2")}, cond),
    s({trig="be.", desc='Atoms in S1 that are at least 12.3 Anstroms away from S2'}, {i(1,"S1"), t" beyond ", i(2,"12.3"), t" of ", i(3,"S2")}, cond),
    s({trig="partial_charge", desc=''}, {t"partial_charge ", i(1,"< 1.2")}, cond),
    s({trig="pc.", desc=''}, {t"pc. ", i(1,"< 1.2")}, cond),
    s({trig="formal_charge", desc=''}, {t"formal_charge ", i(1,"= 1")}, cond),
    s({trig="fc.", desc=''}, {t"fc. ", i(1,"= 1")}, cond),
    s({trig="b", desc='B-factor less than 100.0'}, {t"b ", i(1,"< 100.0")}, cond),
    s({trig="q", desc='Occupancy less than 1.0'}, {t"q ", i(1,"< 1.0")}, cond),
    s({trig="ss", desc='Atoms with secondary structure H (helix) or S (sheet)'}, {t"ss ", i(1,"H+S")}, cond),
    s({trig="elem", desc='Atoms of element, e.g. C (carbon)'}, {t"elem ", i(1,"C")}, cond),
    s({trig="e.", desc='Atoms of element, e.g. C (carbon)'}, {t"e. ", i(1,"C")}, cond),
    s({trig="bonded", desc='Atoms which have at least one bond'}, {t"bonded"}, cond),
    s({trig="protected", desc='Atoms protected from editing with the protect command'}, {t"protected"}, cond),
    s({trig="fixed", desc='Fixed Atoms (no movement allowed)'}, {t"fixed"}, cond),
    s({trig="fxd.", desc='Fixed Atoms (no movement allowed)'}, {t"fxd."}, cond),
    s({trig="restrained", desc='Restrained Atoms (typically harmonically constrained)'}, {t"restrained"}, cond),
    s({trig="rst.", desc='Restrained Atoms (typically harmonically constrained)'}, {t"rst."}, cond),
    s({trig="masked", dec='Atoms made unselectable with mouse by the mask command'}, {t"masked"}, cond),
    s({trig="msk.", dec='Atoms made unselectable with mouse by the mask command'}, {t"msk."}, cond),
    s({trig="flag", desc='Atoms with flag N. 0=focus, 1=free, 2=restrain, 3=fix, 4=exclude, 5=study, ... see https://pymolwiki.org/index.php/Flag'}, {t"flag ", i(1,"25")}, cond),
    s({trig="f.", desc='Atoms with flag N. 0=focus, 1=free, 2=restrain, 3=fix, 4=exclude, 5=study, ... see https://pymolwiki.org/index.php/Flag'}, {t"f. ", i(1,"25")}, cond),
    s({trig="organic", desc='Non-polymer organic compounds (e.g. ligands, buffers)'}, {t"organic"}, cond),
    s({trig="org.", desc='Non-polymer organic compounds (e.g. ligands, buffers)'}, {t"org."}, cond),
    s({trig="inorganic", desc='Non-polymer inorganic atoms/ions'}, {t"inorganic"}, cond),
    s({trig="ino.", desc='Non-polymer inorganic atoms/ions'}, {t"ino."}, cond),
    s({trig="solvent", desc='Water molecules'}, {t"solvent"}, cond),
    s({trig="sol.", desc='Water molecules'}, {t"sol."}, cond),
    s({trig="polymer", desc='Protein or Nucleic Acid'}, {t"polymer"}, cond),
    s({trig="pol.", desc='Protein or Nucleic Acid'}, {t"pol."}, cond),
    s({trig="polymer.protein", desc='Protein'}, {t"polymer.protein"}, cond),
    s({trig="polymer.nucleic", desc='Nucleic Acid'}, {t"polymer.nucleic"}, cond),
    s({trig="guide", desc="Protein CA and nucleic acid C4*/C4'"}, {t"guide"}, cond),
    s({trig="hetatm", desc='Atoms loaded from PDB HETATM records'}, {t"hetatm"}, cond),
    s({trig="hydrogens", desc='Hydrogen atoms'}, {t"hydrogens"}, cond),
    s({trig="h.", desc='Hydrogen atoms'}, {t"h."}, cond),
    s({trig="backbone",	desc='Polymer backbone atoms'}, {t"backbone"}, cond),
    s({trig="bb.",	desc='Polymer backbone atoms'}, {t"bb."}, cond),
    s({trig="sidechain", desc='Polymer non-backbone atoms'}, {t"sidechain"}, cond),
    s({trig="sc.", desc='Polymer non-backbone atoms'}, {t"sc."}, cond),
    s({trig="metals", desc='Metal atoms'}, {t"metals"}, cond),
    s({trig="donors", desc='Hydrogen bond donor atoms'}, {t"donors"}, cond),
    s({trig="don.", desc='Hydrogen bond donor atoms'}, {t"don."}, cond),
    s({trig="acceptors", desc='Hydrogen bond acceptor atoms'}, {t"acceptors"}, cond),
    s({trig="acc.", desc='Hydrogen bond acceptor atoms'}, {t"acc."}, cond),
    s({trig="visible", desc='Atoms in enabled objects with at least one visible representation'}, {t"visible"}, cond),
    s({trig="v.", desc='Atoms in enabled objects with at least one visible representation'}, {t"v."}, cond),
    s({trig="rep", desc='Atoms with given representation'}, {t"rep ", i(1,"cartoon")}, cond),
    s({trig="color", desc='Atoms with given atom-color (by color index)'}, {t"color ", i(1,"blue")}, cond),
    s({trig="cartoon_color", desc='Atoms with given atom-level cartoon_color (by color index)'}, {t"cartoon_color ", i(1,"blue")}, cond),
    s({trig="ribbon_color", desc='Atoms with given atom-level ribbon_color (by color index)'}, {t"ribbon_color", i(1,"blue")}, cond),
    s({trig="center", desc='Pseudo-atom at the center of the scene'}, {t"center"}, cond),
    s({trig="origin", desc='Pseudo-atom at the origin of rotation'}, {t"origin"}, cond),
    s({trig="state", desc='Atoms with coordinates in state N'}, {t"state ", i(1,"123")}, cond),
    s({trig="present", desc='Atoms with coordinates in the current state'}, {t"present"}, cond),
    s({trig="pr.", desc='Atoms with coordinates in the current state'}, {t"pr."}, cond),
    s({trig="x", desc='Atoms with model-space x coordinate less than N'}, {t"x ", i(1,"< 12.3")}, cond),
    s({trig="y", desc='Atoms with model-space y coordinate less than N'}, {t"y ", i(1,"< 12.3")}, cond),
    s({trig="z", desc='Atoms with model-space z coordinate less than N'}, {t"z ", i(1,"< 12.3")}, cond),
    s({trig="numeric_type", desc=''}, {t"numeric_type ", i(1, "123")}, cond),
    s({trig="nt.", desc=''}, {t"nt. ", i(1, "123")}, cond),
}

