#!/usr/bin/env lua

-- show_conditions
local isAttrLine = function (line_to_cursor)
    -- avoid completing after second : in case attribute has already been written
    return vim.startswith(line_to_cursor, ':') and not line_to_cursor:match(':.*:')
end
local isOptionLine = function (line_to_cursor)
    return vim.startswith(line_to_cursor, '[')
end
local isLineBegin = function (line_to_cursor)
    local r, c = unpack(vim.api.nvim_win_get_cursor(0))
    return line_to_cursor:sub(1,c+1):match('^%s*$')
end

local opta = {show_condition=isAttrLine}
local opto = {show_condition=isOptionLine}
local optb = {show_condition=isLineBegin}

return {
-- syntax snippets
-- block comment. [comment] can also precede a paragraph or open block (-- ... --).
s({
    trig="////",
    dscr='Block comment.',
    condition=conds.line_begin, snippetType="autosnippet",
},
{t{"////", ""}, i(1), t{"", "////"}}),
-- https://docs.asciidoctor.org/asciidoc/latest/blocks/open-blocks/
s({
    trig="--",
    dscr='Open block. Precede by [source] to make it source code or .Title to add the title "Title".',
    condition=conds.line_begin, snippetType="autosnippet",
},
{t{"--", ""}, i(1), t{"", "--"}}),
s({
    trig="++++",
    dscr='A delimited passthrough block. Excludes the block’s content from all substitutions unless the subs attribute is set.',
    condition=conds.line_begin, snippetType="autosnippet",
},
{t{"++++", ""}, i(1), t{"", "++++"}}),
-- ref
s({trig="<<", dscr="Cross reference shorthand.", snippetType="autosnippet"},
{t"<<", i(1, "ID"), c(2, {t"", {t",", i(1, "LABEL")}}), t">>"}),
s({trig="xref", dscr="Cross reference.", snippetType="autosnippet"},
{t"xref:", i(1, "FILE.adoc#ID"), t"[", i(2, "LABEL"), t"]"}),
-- include
s({
    trig="inc",
    dscr="Include file.",
    snippetType='autosnippet',
},
{t"include::", i(1, "FILENAME.adoc"), t"[]"}),
-- conditional
s({
    trig="ifdef",
    dscr="If defined",
},
{t"ifdef::", i(1, "ATTR"), t{"[]", ""}, i(2), t{"", "endif::[]"}},
optb),
s({
    trig="ifndef",
    dscr="If not defined",
},
{t"ifndef::", i(1, "ATTR"), t{"[]", ""}, i(2), t{"", "endif::[]"}},
optb),
s({
    trig="ifeval",
    dscr="If evaluates to true",
},
{t"ifeval::[", i(1, "{sectnumlevels} == 3"), t{"]", ""}, i(2), t{"", "endif::[]"}},
optb),
-- doc attributes
-- https://docs.asciidoctor.org/asciidoc/latest/attributes/document-attributes-ref/
s(
    {
        name="backend",
        trig=":backend:",
        dscr="The backend used to select and activate the converter that creates the output file. Usually named according to the output format (e.g., html5).",
    },
    {t":backend: ", i(1, "html5")},
    opta
),
s(
    {
        name="docdate",
        trig=":docdate:",
        dscr="Last modified date of the source document.",
    },
    {t":docdate: ", i(1, "2019-01-04")},
    opta
),
s(
    {
        trig="docdatetime",
        trig=":docdatetime:",
        dscr="Last modified date and time of the source document.",
    },
    {t":docdatetime: ", i(1, "2019-01-04 19:26:06 UTC")},
    opta
),
s(
    {
        name="doctime",
        trig=":doctime:",
        dscr="Last modified time of the source document.",
    },
    {t":doctime: ", i(1, "19:26:06 UTC")},
    opta
),
s(
    {
        name="docyear",
        trig=":docyear:",
        dscr="Year that the document was last modified.",
    },
    {t":docyear: ", i(1, "2023")},
    opta
),
s(
    {
        name="localdate",
        trig=":localdate:",
        dscr="Date when the document was converted.",
    },
    {t":localdate: ", i(1, "2019-02-17")},
    opta
),
s(
    {
        name="localdatetime",
        trig=":localdatetime:",
        dscr="Date and time when the document was converted.",
    },
    {t":localdatetime: ", i(1, "2019-02-17 19:31:05 UTC")},
    opta
),
s(
    {
        name="localtime",
        trig=":localtime:",
        dscr="Time when the document was converted.",
    },
    {t":localtime: ", i(1, "19:31:05 UTC")},
    opta
),
s(
    {
        name="localyear",
        trig=":localyear:",
        dscr="Year when the document was converted.",
    },
    {t":localyear: ", i(1, "2023")},
    opta
),
s(
    {
        name="outfilesuffix",
        trig=":outfilesuffix:",
        dscr="File extension of the output file (starting with a period) as determined by the backend (.html for html, .xml for docbook, etc.).",
    },
    {t":outfilesuffix: ", i(1, ".pdf")},
    opta
),
s(
    {
        name="attribute-missing",
        trig=":attribute-missing:",
        dscr="Controls how missing attribute references are handled. Default=skip.",
        wordTrig=false,
    },
    {t":attribute-missing: ", c(1, {t"drop", t"drop-line", t"skip", t"warn"})},
    opta
),
s(
    {
        name="attribute-undefined",
        trig=":attribute-undefined:",
        dscr="Controls how attribute unassignments are handled. Default=drop-line.",
        wordTrig=false,
    },
    {t":attribute-undefined: ", c(1, {t"drop", t"drop-line"})},
    opta
),
s(
    {
        name="appendix-caption",
        trig=":appendix-caption:",
        dscr="Label added before an appendix title. Default=Appendix.",
        wordTrig=false,
    },
    {t":appendix-caption: ", i(1, "Appendix")},
    opta
),
s(
    {
        name="appendix-number",
        trig=":appendix-number:",
        dscr="Sets the seed value for the appendix number sequence.",
        wordTrig=false,
    },
    {t":appendix-number: ", i(1, "@")},
    opta
),
s(
    {
        name="appendix-refsig",
        trig=":appendix-refsig:",
        dscr="Signifier added to Appendix title cross references.",
        wordTrig=false,
    },
    {t":appendix-refsig: ", i(1, "Appendix")},
    opta
),
s(
    {
        name="caution-caption",
        trig=":caution-caption:",
        dscr="Text used to label CAUTION admonitions when icons aren’t enabled.",
        wordTrig=false,
    },
    {t":caution-caption: ", i(1, "Caution")},
    opta
),
s(
    {
        name="chapter-number",
        trig=":chapter-number:",
        dscr="Sets the seed value for the chapter number sequence. Book doctype only.",
        wordTrig=false,
    },
    {t":chapter-number: ", i(1, "0")},
    opta
),
s(
    {
        name="chapter-refsig",
        trig=":chapter-refsig:",
        dscr="Signifier added to Chapter titles in cross references. Book doctype only.",
        wordTrig=false,
    },
    {t":chapter-refsig: ", i(1, "Chapter")},
    opta
),
s(
    {
        name="chapter-signifier",
        trig=":chapter-signifier:",
        dscr="Label added to level 1 section titles (chapters). Book doctype only.",
        wordTrig=false,
    },
    {t":chapter-signifier: ", i(1, "")},
    opta
),
s(
    {
        name="example-caption",
        trig=":example-caption:",
        dscr="Text used to label example blocks.",
        wordTrig=false,
    },
    {t":example-caption: ", i(1, "Example")},
    opta
),
s(
    {
        name="example-number",
        trig=":example-number:",
        dscr="Sets the seed value for the example number sequence.",
        wordTrig=false,
    },
    {t":example-number: ", i(1, "0")},
    opta
),
s(
    {
        name="figure-caption",
        trig=":figure-caption:",
        dscr="Text used to label images and figures.",
        wordTrig=false,
    },
    {t":figure-caption: ", i(1, "Figure")},
    opta
),
s(
    {
        name="figure-number",
        trig=":figure-number:",
        dscr="Sets the seed value for the figure number sequence.",
        wordTrig=false,
    },
    {t":figure-number: ", i(1, "0")},
    opta
),
s(
    {
        name="footnote-number",
        trig=":footnote-number:",
        dscr="Sets the seed value for the footnote number sequence.",
        wordTrig=false,
    },
    {t":footnote-number: ", i(1, "0")},
    opta
),
s(
    {
        name="important-caption",
        trig=":important-caption:",
        dscr="Text used to label IMPORTANT admonitions when icons are not enabled.",
        wordTrig=false,
    },
    {t":important-caption: ", i(1, "Important")},
    opta
),
s(
    {
        name="lang",
        trig=":lang:",
        dscr="Language tag specified on document element of the output document. Refer to the lang and xml:lang attributes section of the HTML specification to learn about the acceptable values for this attribute. Header only.",
        wordTrig=false,
    },
    {t":lang: ", i(1, "en")},
    opta
),
s(
    {
        name="last-update-label",
        trig=":last-update-label:",
        dscr="Text used for “Last updated” label in footer. Header only.",
        wordTrig=false,
    },
    {t":last-update-label: ", i(1, "Last updated")},
    opta
),
s(
    {
        name="listing-caption",
        trig=":listing-caption:",
        dscr="Text used to label listing blocks.",
        wordTrig=false,
    },
    {t":listing-caption: ", i(1, "")},
    opta
),
s(
    {
        name="listing-number",
        trig=":listing-number:",
        dscr="Sets the seed value for the listing number sequence.",
        wordTrig=false,
    },
    {t":listing-number: ", i(1, "0")},
    opta
),
s(
    {
        name="manname-title",
        trig=":manname-title:",
        dscr="Label for program name section in the man page. Header only.",
        wordTrig=false,
    },
    {t":manname-title: ", i(1, "Name")},
    opta
),
s(
    {
        name="note-caption",
        trig=":note-caption:",
        dscr="Text used to label NOTE admonitions when icons aren’t enabled.",
        wordTrig=false,
    },
    {t":note-caption: ", i(1, "Note")},
    opta
),
s(
    {
        name="part-refsig",
        trig=":part-refsig:",
        dscr="Signifier added to Part titles in cross references. Book doctype only.",
        wordTrig=false,
    },
    {t":part-refsig: ", i(1, "Part")},
    opta
),
s(
    {
        name="part-signifier",
        trig=":part-signifier:",
        dscr="Label added to level 0 section titles (parts). Book doctype only.",
        wordTrig=false,
    },
    {t":part-signifier: ", i(1, "")},
    opta
),
s(
    {
        name="preface-title",
        trig=":preface-title:",
        dscr="Title text for an anonymous preface when doctype is book.",
        wordTrig=false,
    },
    {t":preface-title: ", i(1, "")},
    opta
),
s(
    {
        name="section-refsig",
        trig=":section-refsig:",
        dscr="Signifier added to title of numbered sections in cross reference text.",
        wordTrig=false,
    },
    {t":section-refsig: ", i(1, "Section")},
    opta
),
s(
    {
        name="table-caption",
        trig=":table-caption:",
        dscr="Text of label prefixed to table titles.",
        wordTrig=false,
    },
    {t":table-caption: ", i(1, "Table")},
    opta
),
s(
    {
        name="table-number",
        trig=":table-number:",
        dscr="Sets the seed value for the table number sequence.",
        wordTrig=false,
    },
    {t":table-number: ", i(1, "0")},
    opta
),
s(
    {
        name="tip-caption",
        trig=":tip-caption:",
        dscr="Text used to label TIP admonitions when icons aren’t enabled.",
        wordTrig=false,
    },
    {t":tip-caption: ", i(1, "Tip")},
    opta
),
s(
    {
        name="toc-title",
        trig=":toc-title:",
        dscr="Title for table of contents. Header only.",
        wordTrig=false,
    },
    {t":toc-title: ", i(1, "Table of Contents")},
    opta
),
s(
    {
        name="untitled-label",
        trig=":untitled-label:",
        dscr="Default document title if document doesn’t have a document title. Header only.",
        wordTrig=false,
    },
    {t":untitled-label: ", i(1, "Untitled")},
    opta
),
s(
    {
        name="version-label",
        trig=":version-label:",
        dscr="See Version Label Attribute. Header only.",
        wordTrig=false,
    },
    {t":version-label: ", i(1, "Version")},
    opta
),
s(
    {
        name="warning-caption",
        trig=":warning-caption:",
        dscr="Text used to label WARNING admonitions when icons aren’t enabled.",
        wordTrig=false,
    },
    {t":warning-caption: ", i(1, "Warning")},
    opta
),
s(
    {
        name="app-name",
        trig=":app-name:",
        dscr="Adds application-name meta element for mobile devices inside HTML document head. Header only.",
        wordTrig=false,
    },
    {t":app-name: ", i(1, "")},
    opta
),
s(
    {
        name="author",
        trig=":author:",
        dscr="Can be set automatically via the author info line or explicitly. See Author Information. Header only.",
    },
    {t":author: ", i(1, "Christian Degnbol Madsen")},
    opta
),
s(
    {
        name="authorinitials",
        trig=":authorinitials:",
        dscr="Derived from the author attribute by default. See Author Information. Header only.",
    },
    {t":authorinitials: ", i(1, "CDM")},
    opta
),
s(
    {
        name="authors",
        trig=":authors:",
        dscr="Can be set automatically via the author info line or explicitly as a comma-separated value list. See Author Information. Header only.",
    },
    {t":authors: ", i(1, "")},
    opta
),
s(
    {
        name="copyright",
        trig=":copyright:",
        dscr="Adds copyright meta element in HTML document head. Header only.",
    },
    {t":copyright: ", i(1, "")},
    opta
),
s(
    {
        name="doctitle",
        trig=":doctitle:",
        dscr="Can be set automatically from doctitle attribute. Header only.",
    },
    {t":doctitle: ", i(1, "")},
    opta
),
s(
    {
        name="description",
        trig=":description:",
        dscr="Adds description meta element in HTML document head. Header only.",
    },
    {t":description: ", i(1, "")},
    opta
),
s(
    {
        name="email",
        trig=":email:",
        dscr="Can be any inline macro, such as a URL. Can be set automatically from author info line. Header only.",
    },
    {t":email: ", i(1, "cdmadsen@student.unimelb.edu.au")},
    opta
),
s(
    {
        name="firstname",
        trig=":firstname:",
        dscr="Can be set automatically from author info line. Header only.",
    },
    {t":firstname: ", i(1, "Christian")},
    opta
),
s(
    {
        name="front-matter",
        trig=":front-matter:",
        dscr="If skip-front-matter is set via the API or CLI, any YAML-style frontmatter skimmed from the top of the document is stored in this attribute.",
    },
    {t":front-matter: ", i(1, "")},
    opta
),
s(
    {
        name="keywords",
        trig=":keywords:",
        dscr="Adds keywords meta element in HTML document head. Header only.",
    },
    {t":keywords: ", i(1, "")},
    opta
),
s(
    {
        name="lastname",
        trig=":lastname:",
        dscr="Can be set automatically from author info line. Header only.",
    },
    {t":lastname: ", i(1, "Madsen")},
    opta
),
s(
    {
        name="middlename",
        trig=":middlename:",
        dscr="Can be set automatically from author info line. Header only.",
    },
    {t":middlename: ", i(1, "Degnbol")},
    opta
),
s(
    {
        name="orgname",
        trig=":orgname:",
        dscr="Adds <orgname> element value to DocBook info element. Header only.",
    },
    {t":orgname: ", i(1, "University of Melbourne")},
    opta
),
s(
    {
        name="revdate",
        trig=":revdate:",
        dscr="Can be set automatically from revision info line. Header only.",
    },
    {t":revdate: ", i(1, "")},
    opta
),
s(
    {
        name="revremark",
        trig=":revremark:",
        dscr="Can be set automatically from revision info line. Header only.",
    },
    {t":revremark: ", i(1, "")},
    opta
),
s(
    {
        name="revnumber",
        trig=":revnumber:",
        dscr="Can be set automatically from revision info line. Header only.",
    },
    {t":revnumber: ", i(1, "")},
    opta
),
s(
    {
        name="title",
        trig=":title:",
        dscr="Value of <title> element in HTML <head> or main DocBook <info> of output document. Used as a fallback when the document title is not specified. See title attribute.",
    },
    {t":title: ", i(1, "")},
    opta
),
s(
    {
        name="idprefix",
        trig=":idprefix:",
        dscr="Prefix of auto-generated section IDs. See Change the ID prefix.",
    },
    {t":idprefix: ", i(1, "valid XML ID start character _")},
    opta
),
s(
    {
        name="idseparator",
        trig=":idseparator:",
        dscr="Word separator used in auto-generated section IDs. See Change the ID word separator.",
    },
    {t":idseparator: ", i(1, "valid XML ID character _")},
    opta
),
s(
    {
        name="leveloffset",
        trig=":leveloffset:",
        dscr="Increases or decreases level of headings below assignment. A leading + or - makes the value relative.",
    },
    {t":leveloffset: ", i(1, "[+-]0–5")},
    opta
),
s(
    {
        name="sectnumlevels",
        trig=":sectnumlevels:",
        dscr="Controls depth of section numbering. Default=3.",
    },
    {t":sectnumlevels: ", i(1, "0–5")},
    opta
),
s(
    {
        name="title-separator",
        trig=":title-separator:",
        dscr="Character used to separate document title and subtitle. Header only.",
        wordTrig=false,
    },
    {t":title-separator: ", i(1, "")},
    opta
),
s(
    {
        name="toc",
        trig=":toc:",
        dscr="Turns on table of contents and specifies its location. Default=auto. Header only.",
    },
    {t":toc: ", c(1, {t"", t"auto", t"left", t"right", t"macro", t"preamble"})},
    opta
),
s(
    {
        name="toclevels",
        trig=":toclevels:",
        dscr="Maximum section level to display. Default=2. Header only.",
    },
    {t":toclevels: ", i(1, "1–5")},
    opta
),
s(
    {
        name="asset-uri-scheme",
        trig=":asset-uri-scheme:",
        dscr="Controls protocol used for assets hosted on a CDN. Header only.",
        wordTrig=false,
    },
    {t":asset-uri-scheme: ", c(1, {t"http", t"https"})},
    opta
),
s(
    {
        name="docinfo",
        trig=":docinfo:",
        dscr="Read input from one or more DocBook info files. No arg=private. Header only.",
    },
    {t":docinfo: ", c(1, {t"", t"shared", t"private", t"shared-head", t"private-head", t"shared-footer", t"private-footer"})},
    opta
),
s(
    {
        name="docinfodir",
        trig=":docinfodir:",
        dscr="Location of docinfo files. Defaults to directory of source file if not specified. Header only.",
    },
    {t":docinfodir: ", i(1, "./")},
    opta
),
s(
    {
        name="docinfosubs",
        trig=":docinfosubs:",
        dscr="AsciiDoc substitutions that are applied to docinfo content. Header only.",
    },
    {t":docinfosubs: ", i(1, "comma-separated substitution names (attributes)")},
    opta
),
s(
    {
        name="doctype",
        trig=":doctype:",
        dscr="Output document type. Default=article. Header only.",
    },
    {t":doctype: ", c(1, {t"article", t"book", t"inline", t"manpage"})},
    opta
),
s(
    {
        name="eqnums",
        trig=":eqnums:",
        dscr=[[Controls automatic equation numbering on LaTeX blocks in HTML output (MathJax). If the value is AMS, only LaTeX content enclosed in an \begin{equation, condition=conda}...\end{equation} container will be numbered. If the value is all, then all LaTeX blocks will be numbered. See equation numbering in MathJax. Default=AMS. Header only.]],
    },
    {t":eqnums: ", c(1, {t"", t"AMS", t"all", t"none"})},
    opta
),
s(
    {
        name="media",
        trig=":media:",
        dscr="Specifies media type of output and enables behavior specific to that media type. PDF converter only. Header only.",
    },
    {t":media: ", c(1, {t"prepress", t"print", t"screen"})},
    opta
),
s(
    {
        name="outfilesuffix",
        trig=":outfilesuffix:",
        dscr="File extension of output file, including dot (.), such as .html. Header only.",
    },
    {t":outfilesuffix: ", i(1, ".html")},
    opta
),
s(
    {
        name="pagewidth",
        trig=":pagewidth:",
        dscr="Page width used to calculate the absolute width of tables in the DocBook output. Header only.",
    },
    {t":pagewidth: ", i(1, "425")},
    opta
),
s(
    {
        name="relfilesuffix",
        trig=":relfilesuffix:",
        dscr="The path suffix (e.g., file extension) to add to relative xrefs. Defaults to the value of the outfilesuffix attribute. (Preferred over modifying outfilesuffix).",
    },
    {t":relfilesuffix: ", i(1, ".html")},
    opta
),
s(
    {
        name="stem",
        trig=":stem:",
        dscr="Enables mathematics processing and interpreter. No arg=asciimath. Header only.",
    },
    {t":stem: ", c(1, {t"", t"asciimath", t"latexmath"})},
    opta
),
s(
    {
        name="table-frame",
        trig=":table-frame:",
        dscr="Controls default value for frame attribute on tables.",
        wordTrig=false,
    },
    {t":table-frame: ", c(1, {t"all", t"ends", t"sides", t"none"})},
    opta
),
s(
    {
        name="table-grid",
        trig=":table-grid:",
        dscr="Controls default value for grid attribute on tables.",
        wordTrig=false,
    },
    {t":table-grid: ", c(1, {t"all", t"cols", t"rows", t"none"})},
    opta
),
s(
    {
        name="table-stripes",
        trig=":table-stripes:",
        dscr="Controls default value for stripes attribute on tables.",
        wordTrig=false,
    },
    {t":table-stripes: ", c(1, {t"none", t"even", t"odd", t"hover", t"all"})},
    opta
),
s(
    {
        name="tabsize",
        trig=":tabsize:",
        dscr="Converts tabs to spaces in verbatim content blocks (e.g., listing, literal).",
    },
    {t":tabsize: ", i(1, "≥0")},
    opta
),
s(
    {
        name="xrefstyle",
        trig=":xrefstyle:",
        dscr="Formatting style to apply to cross reference text.",
    },
    {t":xrefstyle: ", c(1, {t"full", t"short", t"basic"})},
    opta
),
s(
    {
        name="iconfont-cdn",
        trig=":iconfont-cdn:",
        dscr="If not specified, uses the cdnjs.com service. Overrides CDN used to link to the Font Awesome stylesheet. Header only.",
        wordTrig=false,
    },
    {t":iconfont-cdn: ", i(1, "URL")},
    opta
),
s(
    {
        name="iconfont-name",
        trig=":iconfont-name:",
        dscr="Overrides the name of the icon font stylesheet. Header only.",
        wordTrig=false,
    },
    {t":iconfont-name: ", i(1, "font-awesome")},
    opta
),
s(
    {
        name="icons",
        trig=":icons:",
        dscr="Chooses images or font icons instead of text for admonitions. Any other value is assumed to be an icontype and sets the value to empty (image-based icons). Header only.",
    },
    {t":icons: ", c(1, {t"", t"image", t"font"})},
    opta
),
s(
    {
        name="iconsdir",
        trig=":iconsdir:",
        dscr="Location of non-font-based image icons. Defaults to the icons folder under imagesdir if imagesdir is specified and iconsdir is not specified.",
    },
    {t":iconsdir: ", i(1, "./images/icons")},
    opta
),
s(
    {
        name="icontype",
        trig=":icontype:",
        dscr="File type for image icons. Only relevant when using image-based icons.",
    },
    {t":icontype: ", c(1, {t"jpg", t"png", t"gif", t"svg"})},
    opta
),
s(
    {
        name="icon-set",
        trig=":icon-set:",
        dscr="Choose font-based icon set from font awesome solid, brands, regular or foundations icons.",
        wordTrig=false,
    },
    {t":icon-set: ", c(1, {t"fas", t"fab", t"far", t"fi"})},
    opta
),
s(
    {
        name="coderay-css",
        trig=":coderay-css:",
        dscr="Controls whether CodeRay uses CSS classes or inline styles. Header only.",
        wordTrig=false,
    },
    {t":coderay-css: ", c(1, {t"class", t"style"})},
    opta
),
s(
    {
        name="coderay-linenums-mode",
        trig=":coderay-linenums-mode:",
        dscr="Sets how CodeRay inserts line numbers into source listings.",
        wordTrig=false,
    },
    {t":coderay-linenums-mode: ", c(1, {t"inline", t"table"})},
    opta
),
s(
    {
        name="highlightjs-theme",
        trig=":highlightjs-theme:",
        dscr="Name of theme used by highlight.js. Header only.",
        wordTrig=false,
    },
    {t":highlightjs-theme: ", i(1, "github")},
    opta
),
s(
    {
        name="prettify-theme",
        trig=":prettify-theme:",
        dscr="Name of theme used by prettify. Header only.",
        wordTrig=false,
    },
    {t":prettify-theme: ", i(1, "prettify")},
    opta
),
s(
    {
        name="pygments-css",
        trig=":pygments-css:",
        dscr="Controls whether Pygments uses CSS classes or inline styles. Header only.",
        wordTrig=false,
    },
    {t":pygments-css: ", c(1, {t"class", t"style"})},
    opta
),
s(
    {
        name="pygments-linenums-mode",
        trig=":pygments-linenums-mode:",
        dscr="Sets how Pygments inserts line numbers into source listings.",
        wordTrig=false,
    },
    {t":pygments-linenums-mode: ", c(1, {t"table", t"inline"})},
    opta
),
s(
    {
        name="pygments-style",
        trig=":pygments-style:",
        dscr="Name of style used by Pygments. Header only.",
        wordTrig=false,
    },
    {t":pygments-style: ", i(1, "default")},
    opta
),
s(
    {
        name="rouge-css",
        trig=":rouge-css:",
        dscr="Controls whether Rouge uses CSS classes or inline styles. Header only.",
        wordTrig=false,
    },
    {t":rouge-css: ", i(1, {t"class", t"style"})},
    opta
),
s(
    {
        name="rouge-linenums-mode",
        trig=":rouge-linenums-mode:",
        dscr="Sets how Rouge inserts line numbers into source listings. `inline` not yet supported by Asciidoctor. See asciidoctor#3641.",
        wordTrig=false,
    },
    {t":rouge-linenums-mode: ", c(1, {t"inline", t"table"})},
    opta
),
s(
    {
        name="rouge-style",
        trig=":rouge-style:",
        dscr="Name of style used by Rouge. Header only.",
        wordTrig=false,
    },
    {t":rouge-style: ", i(1, "github")},
    opta
),
s(
    {
        name="source-highlighter",
        trig=":source-highlighter:",
        dscr="Specifies source code highlighter. Any other value is permitted, but must be supported by a custom syntax highlighter adapter. Header only.",
        wordTrig=false,
    },
    {t":source-highlighter: ", c(1, {t"coderay", t"highlight.js", t"pygments", t"rouge"})},
    opta
),
s(
    {
        name="source-indent",
        trig=":source-indent:",
        dscr="Normalize block indentation in source code listings.",
        wordTrig=false,
    },
    {t":source-indent: ", i(1, "4")},
    opta
),
s(
    {
        name="source-language",
        trig=":source-language:",
        dscr="Default language for source code blocks.",
        wordTrig=false,
    },
    {t":source-language: ", c(1, {t"julia", t"python", t"R", t"zsh"})},
    opta
),
s(
    {
        name="css-signature",
        trig=":css-signature:",
        dscr="Assign value to id attribute of HTML <body> element. Preferred approach is to assign an ID to document title. Header only.",
        wordTrig=false,
    },
    {t":css-signature: ", i(1, "XML ID")},
    opta
),
s(
    {
        name="max-width",
        trig=":max-width:",
        dscr="Constrains maximum width of document body. Not recommended. Use CSS stylesheet instead. Header only.",
        wordTrig=false,
    },
    {t":max-width: ", i(1, "55em")},
    opta
),
s(
    {
        name="man-linkstyle",
        trig=":man-linkstyle:",
        dscr="Link style in man page output. Header only.",
        wordTrig=false,
    },
    {t":man-linkstyle: ", c(1, {t"blue", t"R", t"<>"})},
    opta
),
s(
    {
        name="page-background-image",
        trig=":page-background-image:",
        dscr="Page background image.",
        wordTrig=false,
    },
    {t":page-background-image: image:", i(1, "FILENAME"), t"[", i(2), t"]"},
    opta
),
s(
    {
        name="page-background-image-recto",
        trig=":page-background-image-recto:",
        dscr="Recto page background image.",
        wordTrig=false,
    },
    {t":page-background-image-recto: image:", i(1, "FILENAME"), t"[", i(2), t"]"},
    opta
),
s(
    {
        name="page-background-image-verso",
        trig=":page-background-image-verso:",
        dscr="Verso page background image.",
        wordTrig=false,
    },
    {t":page-background-image-verso: image:", i(1, "FILENAME"), t"[", i(2), t"]"},
    opta
),
s(
    {
        name="title-page-background-image",
        trig=":title-page-background-image:",
        dscr="Title page background image.",
        wordTrig=false,
    },
    {t":title-page-background-image: ", c(1, {t"none", {t"image:", i(1, "FILENAME"), t"[", i(2), t"]"}})},
    opta
),
-- options
s({trig="float", dscr="Float alignment, left or right."},
{t"float=", c(1, {t"left", t"right"})},
opto),
s({trig="align", dscr="Text alignment, left, center, or right."},
{t"align=", c(1, {t"left", t"center", t"right"})},
opto),
-- macros
s({trig="pass", dscr="Inline passthrough macro."},
{t"pass:[", i(1), t"]"}),
s({trig="image", dscr="Image."},
{
    t"image::", i(1, "FILENAME"),
    t"[", c(2, {
        t"",
         i(1, "TITLE"),
        {i(1, "TITLE"), t",", i(2, "640"), t",", i(3, "480")},
        {t"width=", i(2, "640"), t",height=", i(3, "480")},
    }), t"]",
}, {condition=conds.line_begin}),
-- Math. Use the :stem: attribute
s(
    {
        trig="mm",
        dscr="Inline math, either asciimath or latexmath depending on :stem: attribute.",
        snippetType="autosnippet",
    },
    {t"stem:[", i(1), t"]"}
),
s(
    {
        trig="dm",
        dscr="Block math, either asciimath or latexmath depending on :stem: attribute.",
        snippetType="autosnippet",
        condition=conds.line_begin,
    },
    {t{"[stem]", "++++", ""}, i(1), t{"", "++++"}}
),
-- tables
s(
    {
        trig="tab",
        dscr="Table.",
        snippetType="autosnippet",
        condition=conds.line_begin,
    },
    {
        t"[", c(1, {
            t'%autowidth',
            {t'cols="', i(1, "<1e,^1,>2"), t'"'},
        }), t{']',
        '|===', "|"},
        i(2, "H1|H2|H3"), t{"", "", "|"},
        i(3, {"c1|c2|c3","","|c4","|c5","|c6"}), t{"", "|==="}
    }
),

}
