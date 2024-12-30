local M = {}

local registered = false

M.setup = function()
    if registered then return end
    registered = true

    local has_cmp, cmp = pcall(require, 'cmp')
    if not has_cmp then return end

    -- the rest is modified from :help cmp-develop

    local source = {}

    ---Return whether this source is available in the current context or not (optional).
    ---@return boolean
    function source:is_available()
        return true
    end

    ---Return the keyword pattern for triggering completion (optional).
    ---If this is ommited, nvim-cmp will use a default keyword pattern. See |cmp-config.completion.keyword_pattern|.
    ---@return string
    function source:get_keyword_pattern()
        -- add - to pattern so we allow dashes within attribute names
        return [[\%(\k\|\-\)\+]]
    end

    ---Return trigger characters for triggering completion (optional).
    function source:get_trigger_characters()
        -- NOTE default is ., we tried : ! and % which we consider trigger chars, but are not a part of the completed attribute or option word.
        -- However there were issues with snippets not being listed until regular completions were depleted, e.g. typing 
        -- :ema would show :experimental and :emai would finally update the listing and give you :email. 
        -- A half fix was to have keyword_length=2 for asciidoc and 1 for luasnip. However setting it to 2 means %header will become annoying.
        -- The correct solution was to remove : and ! from trigger chars as to not give regular completion items that "advantage".
        return { '%' }
    end

    -- See other recognized fields at nvim-cmp/lua/cmp/types/lsp.lua
    -- https://github.com/hrsh7th/nvim-cmp/blob/5dce1b778b85c717f6614e3f4da45e9f19f54435/lua/cmp/types/lsp.lua
    -- E.g. "detail" which is code from the filetype, in this case, asciidoc and "documentation" which is following lines of plaintext.
    local attributes = {
        {label="embedded",               insertText="embedded:",               documentation="Only set if content is being converted to an embedded document (i.e., body of document only)."},
        {label="safe-mode-unsafe",       insertText="safe-mode-unsafe:",       documentation="Set if the safe mode is UNSAFE. Read-only."},
        {label="safe-mode-safe",         insertText="safe-mode-safe:",         documentation="Set if the safe mode is SAFE. Read-only."},
        {label="safe-mode-server",       insertText="safe-mode-server:",       documentation="Set if the safe mode is SERVER. Read-only."},
        {label="safe-mode-secure",       insertText="safe-mode-secure:",       documentation="Set if the safe mode is SECURE. Read-only."},
        {label="safe-mode-level",        insertText="safe-mode-level:",        documentation="Numeric value of the safe mode setting. Read-only. (0=UNSAFE, 1=SAFE, 10=SERVER, 20=SECURE)."},
        {label="safe-mode-name",         insertText="safe-mode-name:",         documentation="Textual value of the safe mode setting. Read-only. (UNSAFE, SAFE, SERVER, SECURE)"},
        {label="compat-mode",            insertText="compat-mode:",            documentation="Controls when the legacy parsing mode is used to parse the document."},
        {label="experimental",           insertText="experimental:",           documentation="Enables Button and Menu UI Macros and the Keyboard Macro. Header only."},
        {label="reproducible",           insertText="reproducible:",           documentation="Prevents last-updated date from being added to HTML footer or DocBook info element. Useful for storing the output in a source code control system as it prevents spurious changes every time you convert the document. Alternately, you can use the SOURCE_DATE_EPOCH environment variable, which sets the epoch of all source documents and the local datetime to a fixed value. Header only."},
        {label="skip-front-matter",      insertText="skip-front-matter:",      documentation="Consume YAML-style frontmatter at top of document and store it in front-matter attribute. Header only."},
        {label="nolang",                 insertText="nolang:",                 documentation="Prevents lang attribute from being added to root element of the output document. Header only."},
        {label="partnums",               insertText="partnums:",               documentation="Enables numbering of parts. See Number book parts. Book doctype only."},
        {label="sectanchors",            insertText="sectanchors:",            documentation="Adds anchor in front of section title on mouse cursor hover."},
        {label="sectids",                insertText="sectids:",                documentation="Generates and assigns an ID to any section that does not have an ID. Set by default."},
        {label="sectlinks",              insertText="sectlinks:",              documentation="Turns section titles into self-referencing links."},
        {label="sectnums",               insertText="sectnums:",               documentation="Enable sections numbering to depth specified by sectnumlevels."},
        {label="fragment",               insertText="fragment:",               documentation="Informs parser that document is a fragment and that it shouldn’t enforce proper section nesting. Header only."},
        {label="cache-uri",              insertText="cache-uri:",              documentation="Cache content read from URIs. Header only."},
        {label="data-uri",               insertText="data-uri:",               documentation="Embed graphics as data-uri elements in HTML elements so file is completely self-contained. Header only."},
        {label="hardbreaks-option",      insertText="hardbreaks-option:",      documentation="Preserve hard line breaks."},
        {label="hide-uri-scheme",        insertText="hide-uri-scheme:",        documentation="Hides URI scheme for raw links."},
        {label="nofooter",               insertText="nofooter:",               documentation="Turns off footer. Header only."},
        {label="nofootnotes",            insertText="nofootnotes:",            documentation="Turns off footnotes. Header only."},
        {label="noheader",               insertText="noheader:",               documentation="Turns off header. Header only."},
        {label="notitle",                insertText="notitle:",                documentation="Hides the doctitle in an embedded document. Mutually exclusive with the showtitle attribute. Header only."},
        {label="relfileprefix",          insertText="relfileprefix:",          documentation="The path prefix to add to relative xrefs."},
        {label="show-link-uri",          insertText="show-link-uri:",          documentation="Prints the URI of a link after the link text. PDF converter only."},
        {label="showtitle",              insertText="showtitle:",              documentation="Displays the doctitle in an embedded document. Mutually exclusive with the notitle attribute. Header only."},
        {label="webfonts",               insertText="webfonts:",               documentation="Control whether webfonts are loaded when using the default stylesheet. When set to empty, uses the default font collection from Google Fonts. A non-empty value replaces the family query string parameter in the Google Fonts URL. Header only."},
        {label="iconfont-remote",        insertText="iconfont-remote:",        documentation="Allows use of a CDN for resolving the icon font. Only relevant used when value of icons attribute is font. Header only. Set by default."},
        {label="imagesdir ",             insertText="imagesdir: ",             documentation="Location of image files."},
        {label="coderay-unavailable",    insertText="coderay-unavailable:",    documentation="Instructs processor not to load CodeRay. Also set if processor fails to load CodeRay. Header only."},
        {label="prewrap",                insertText="prewrap:",                documentation="Wrap wide code listings. Set by default."},
        {label="pygments-unavailable",   insertText="pygments-unavailable:",   documentation="Instructs processor not to load Pygments. Also set if processor fails to load Pygments. Header only."},
        {label="rouge-unavailable",      insertText="rouge-unavailable:",      documentation="Instructs processor not to load Rouge. Also set if processor fails to load Rouge. Header only."},
        {label="source-linenums-option", insertText="source-linenums-option:", documentation="Turns on line numbers for source code listings."},
        {label="copycss",                insertText="copycss:",                documentation="Copy CSS files to output. Only relevant when the linkcss attribute is set. Header only. Set by default."},
        {label="linkcss",                insertText="linkcss:",                documentation="Links to stylesheet instead of embedding it. Can’t be unset in SECURE mode. Header only."},
        {label="stylesheet",             insertText="stylesheet: ",            documentation="CSS stylesheet file name. An empty value tells the converter to use the default stylesheet. Header only."},
        {label="basebackend",	         insertText="basebackend:",	           documentation="The generic backend on which the backend is based. Typically derived from the backend value minus trailing numbers (e.g., the basebackend for docbook5 is docbook). May also indicate the internal backend employed by the converter (e.g., the basebackend for pdf is html)."},
        {label="user-home",	             insertText="user-home:",	           documentation="Get the full path of the home directory for the current user. Masked as . if the safe mode is SERVER or SECURE."},
        {label="mantitle",               insertText="mantitle:",               documentation="Metadata for man page output. Header only."},
        {label="manvolnum",              insertText="manvolnum:",              documentation="Metadata for man page output. Header only."},
        {label="manname",                insertText="manname:",                documentation="Metadata for man page output. Header only."},
        {label="manpurpose",             insertText="manpurpose:",             documentation="Metadata for man page output. Header only."},
        {label="mansource",              insertText="mansource:",              documentation="Source (e.g., application and version) the man page describes. Header only."},
        {label="manmanual",              insertText="manmanual:",              documentation="Manual name displayed in the man page footer. Header only."},
        {label="outfile",                insertText="outfile:",                documentation="Full path of the output file. (Cannot be referenced in the content. Only available to the API once the document is converted). Read-only."},
        {label="outdir",                 insertText="outdir:",                 documentation="Full path of the output directory. (Cannot be referenced in the content. Only available to the API once the document is converted)."},
        {label="htmlsyntax",             insertText="htmlsyntax:",             documentation="Syntax used when generating the HTML output. Controlled by and derived from the backend name (html=html or xhtml=html)."},
        {label="relfileprefix",          insertText="relfileprefix: ",         documentation="The path prefix to add to relative xrefs."},
        {label="highlightjsdir",         insertText="highlightjsdir: ",        documentation="Location of the highlight.js source code highlighter library. Default=CDN URL. Header only."},
        {label="prettifydir",            insertText="prettifydir: ",           documentation="Location of non-CDN prettify source code highlighter library. Default=CDN URL. Header only."},
        {label="stylesdir",              insertText="stylesdir: ",             documentation="Location of CSS stylesheets. Header only."},
        {label="toc-class",              insertText="toc-class: ",             documentation="CSS class on the table of contents container. HTML only. Header only."},
        {label="figure-caption",         insertText="figure-caption:",         documentation="Enable figure captions. Enabled by default."},
        {label="table-caption",          insertText="table-chaption:",         documentation='Enable table captions. Enabled by default. If a title is given (.Title) then it will still show but there will be no label, i.e. "Table A.".'},
        -- https://docs.asciidoctor.org/pdf-converter/latest/theme/apply-theme/
        {label="pdf-theme",              insertText="pdf-theme: ",             documentation='The name or file path of the theme to load. Can be set using the --theme CLI option. When using JRuby, the file path may begin with uri:classloader: to reference a location on the classpath.'},
        {label="pdf-themesdir",          insertText="pdf-themesdir: ",         documentation='The directory path where the theme file is located. When using JRuby, the directory path may begin with uri:classloader: to reference a location on the classpath.'},
        {label="pdf-fontsdir",           insertText="pdf-fontsdir: ",          documentation='The directory path or paths where the fonts used by your theme, if any, are located. When using JRuby, each path may begin with uri:classloader: to reference a location on the classpath.'},
    }
    local options = {
        {label="discrete",                                                     documentation='A discrete heading is declared and styled in a manner similar to that of a section title, but it’s not part of the section hierarchy, it cannot have any child blocks, and it’s not included in the table of contents.'},
        {label="options",                insertText="options=",                documentation='Comma-separated list of options in double quotes, e.g. "unbreakable" to disallow table span across page break.'},
        {label="source",                                                       documentation='Set block style as source code (verbatim).'},
        {label="comment",                                                      documentation='Comment out the block.'},
        {label="quote",                                                        documentation='Set block as quote. Optionally follow by arguments attribution, then citation title and information.'},
        {label="pass",                                                         documentation='The pass style and delimited passthrough block exclude the block’s content from all substitutions unless the subs attribute is set.'},
        {label="python",                                                       documentation='Mark source code as Python.'},
        {label="julia",                                                        documentation='Mark source code as Julia.'},
        {label="subs",                   insertText="subs=",                   documentation='You can use the subs attribute to specify a comma-separated list of substitutions. These substitutions will be applied to the content prior to it being reintroduced to the output document.'},
        {label="indent",                 insertText="indent=",                 documentation='Indent text in block.'},
        -- https://docs.asciidoctor.org/asciidoc/latest/macros/image-ref/
        {label="alt",                    insertText="alt=",                    documentation='Alternative text.'},
        {label="title",                  insertText="title=",                  documentation='Title text, e.g. the entire caption for a table or figure.'},
        {label="caption",                insertText="caption=",                documentation='Caption text, i.e. a table or figure caption excluding the label and numbering.'},
        {label="width",                  insertText="width=",                  documentation='Image width.'},
        {label="pdfwidth",               insertText="pdfwidth=",               documentation='The preferred width of the image in the PDF when converting using Asciidoctor PDF.'},
        {label="scaledwidth",            insertText="scaledwidth=",            documentation='The preferred width of the image when converting to PDF using the DocBook toolchain. (Mutually exclusive with scale).'},
        {label="scale",                  insertText="scale=",                  documentation='Scales the original image size by this amount when converting to PDF using the DocBook toolchain. (Mutually exclusive with scaledwidth).'},
        {label="height",                 insertText="height=",                 documentation='Image height.'},
        {label="float",                  insertText="float=",                  documentation='Float "left" or "right".'},
        {label="align",                  insertText="align=",                  documentation='Align text "left", "right", or "center".'},
        {label="format",                 insertText="format=",                 documentation='Image format, e.g. svg.'},
        -- decorations
        {label="underline",                                                    documentation='Underline text with [underline]#text#. Can be combined with e.g. bold or italic by using * or _ instead of #.'},
        {label="line-through",                                                 documentation='Line through text with [line-through]#text#. Can be combined with e.g. bold or italic by using * or _ instead of #.'},
        {label="none",                                                         documentation='Clears an inherited value and no decoration is applied to the text.'},
        {label="big",                                                          documentation='Big text size.'},
        {label="separator",              insertText="separator=",              documentation='Set a custom separator character for a table. Default=| for top-level tables and ! for nested. Also changed by setting format option, e.g. comma for CSV.'},
    }
    local shorthands = {
        {label="header",                                                       documentation='Set first row of table as a header row. Implicitly set when first row is written on one line. Shorthand syntax for options="header".'},
        {label="footer",                                                       documentation='Set last row of table as a footer row.'},
        {label="collapsible",                                                  documentation='Make example collapsible (in html).'},
        {label="linenums",                                                     documentation='Show line numbers, e.g. for code.'},
        {label="autowidth",                                                    documentation='Set column widths automatically by content in a table.'},
        {label="autowidth.stretch",                                            documentation='Set column widths automatically by content in a table and stretch the table to fill the available page width.'},
        {label="rotate",                                                       documentation='Display table in landscape orientation by rotating it 90 degrees counterclockwise.'},
    }

    -- main purpose of this is to separate these completion items from regular text items, 
    -- so regular text suggestions comes after due to my custom sorting in completion.lua
    local types = require('cmp.types')
    for _, item in ipairs(attributes) do
        item["kind"] = types.lsp.CompletionItemKind.Property
    end
    for _, item in ipairs(options) do
        item["kind"] = types.lsp.CompletionItemKind.Variable
    end
    for _, item in ipairs(shorthands) do
        item["kind"] = types.lsp.CompletionItemKind.Property
    end

    ---Invoke completion (required).
    ---@param params cmp.SourceCompletionApiParams
    ---@param callback fun(response: lsp.CompletionResponse|nil)
    function source:complete(params, callback)
        local line = params.context.cursor_before_line
        local prefix = line:sub(1, params.offset-1)
        if prefix == ":" or prefix == ":!" then
            callback {items=attributes, isIncomplete=true}
        elseif vim.startswith(line, "[") then
            if vim.endswith(prefix, "%") then
                if not prefix:match("[0-9]%%$") then
                    callback {items=shorthands, isIncomplete=false}
                end
            else
                callback {items=options, isIncomplete=true}
            end
        else
            -- isIncomplete allows adding other completion items, such as snippets.
            callback {isIncomplete=true}
        end
    end

    ---Resolve completion item (optional). This is called right before the completion is about to be displayed.
    ---Useful for setting the text shown in the documentation window (`completion_item.documentation`).
    ---@param completion_item lsp.CompletionItem
    ---@param callback fun(completion_item: lsp.CompletionItem|nil)
    function source:resolve(completion_item, callback)
        callback(completion_item)
    end

    ---Executed after the item was selected.
    ---@param completion_item lsp.CompletionItem
    ---@param callback fun(completion_item: lsp.CompletionItem|nil)
    function source:execute(completion_item, callback)
        callback(completion_item)
    end

    ---Register your source to nvim-cmp.
    require('cmp').register_source('asciidoc', source)

    cmp.setup.filetype('asciidoc', {
        sources = cmp.config.sources({
            { name = 'luasnip', options = {show_autosnippets=true} },
            { name = 'asciidoc' },
            { name = 'path', option = {trailing_slash=true} },
            -- { name = 'calc' },
            { name = 'buffer', keyword_length=2, group_index=2 },
            { name = 'dictionary', keyword_length=3, max_item_count=10, group_index=2 },
        }),
    })
end

return M
