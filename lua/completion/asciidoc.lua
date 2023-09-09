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

    ---Return the debug name of this source (optional).
    ---@return string
    function source:get_debug_name()
        return 'debug name'
    end

    ---Return LSP's PositionEncodingKind.
    ---@NOTE: If this method is ommited, the default value will be `utf-16`.
    ---@return lsp.PositionEncodingKind
    function source:get_position_encoding_kind()
        return 'utf-16'
    end

    ---Return the keyword pattern for triggering completion (optional).
    ---If this is ommited, nvim-cmp will use a default keyword pattern. See |cmp-config.completion.keyword_pattern|.
    ---@return string
    function source:get_keyword_pattern()
        -- add !, % and - to pattern so we allow e.g. :sectnums!: and shorthand syntax %header
        return [[\%(\k\|\!\|\-\|%\)\+]]
    end

    ---Return trigger characters for triggering completion (optional).
    function source:get_trigger_characters()
        return { ':%' }
    end

    -- See other recognized fields at nvim-cmp/lua/cmp/types/lsp.lua
    -- https://github.com/hrsh7th/nvim-cmp/blob/5dce1b778b85c717f6614e3f4da45e9f19f54435/lua/cmp/types/lsp.lua
    -- E.g. "detail" which is first lines of extra info and "documentation" which is following lines.
    ---Invoke completion (required).
    ---@param params cmp.SourceCompletionApiParams
    ---@param callback fun(response: lsp.CompletionResponse|nil)
    function source:complete(params, callback)
        callback {
            {label=":embedded:",               detail="Only set if content is being converted to an embedded document (i.e., body of document only)."},
            {label=":safe-mode-unsafe:",       detail="Set if the safe mode is UNSAFE. Read-only."},
            {label=":safe-mode-safe:",         detail="Set if the safe mode is SAFE. Read-only."},
            {label=":safe-mode-server:",       detail="Set if the safe mode is SERVER. Read-only."},
            {label=":safe-mode-secure:",       detail="Set if the safe mode is SECURE. Read-only."},
            {label=":safe-mode-level:",        detail="Numeric value of the safe mode setting. Read-only. (0=UNSAFE, 1=SAFE, 10=SERVER, 20=SECURE)."},
            {label=":safe-mode-name:",         detail="Textual value of the safe mode setting. Read-only. (UNSAFE, SAFE, SERVER, SECURE)"},
            {label=":compat-mode:",            detail="Controls when the legacy parsing mode is used to parse the document."},
            {label=":experimental:",           detail="Enables Button and Menu UI Macros and the Keyboard Macro. Header only."},
            {label=":reproducible:",           detail="Prevents last-updated date from being added to HTML footer or DocBook info element. Useful for storing the output in a source code control system as it prevents spurious changes every time you convert the document. Alternately, you can use the SOURCE_DATE_EPOCH environment variable, which sets the epoch of all source documents and the local datetime to a fixed value. Header only."},
            {label=":skip-front-matter:",      detail="Consume YAML-style frontmatter at top of document and store it in front-matter attribute. Header only."},
            {label=":nolang:",                 detail="Prevents lang attribute from being added to root element of the output document. Header only."},
            {label=":partnums:",               detail="Enables numbering of parts. See Number book parts. Book doctype only."},
            {label=":sectanchors:",            detail="Adds anchor in front of section title on mouse cursor hover."},
            {label=":sectids:",                detail="Generates and assigns an ID to any section that does not have an ID. Set by default."},
            {label=":!sectids:",               detail="Disable automatic ID generation. "},
            {label=":sectlinks:",              detail="Turns section titles into self-referencing links."},
            {label=":sectnums:",               detail="Numbers sections to depth specified by sectnumlevels."},
            {label=":sectnums!:",              detail="Disable section numbering from this line onwards."},
            {label=":fragment:",               detail="Informs parser that document is a fragment and that it shouldn’t enforce proper section nesting. Header only."},
            {label=":cache-uri:",              detail="Cache content read from URIs. Header only."},
            {label=":data-uri:",               detail="Embed graphics as data-uri elements in HTML elements so file is completely self-contained. Header only."},
            {label=":hardbreaks-option:",      detail="Preserve hard line breaks."},
            {label=":hide-uri-scheme:",        detail="Hides URI scheme for raw links."},
            {label=":nofooter:",               detail="Turns off footer. Header only."},
            {label=":nofootnotes:",            detail="Turns off footnotes. Header only."},
            {label=":noheader:",               detail="Turns off header. Header only."},
            {label=":notitle:",                detail="Hides the doctitle in an embedded document. Mutually exclusive with the showtitle attribute. Header only."},
            {label=":relfileprefix:",          detail="The path prefix to add to relative xrefs."},
            {label=":show-link-uri:",          detail="Prints the URI of a link after the link text. PDF converter only."},
            {label=":showtitle:",              detail="Displays the doctitle in an embedded document. Mutually exclusive with the notitle attribute. Header only."},
            {label=":webfonts:",               detail="Control whether webfonts are loaded when using the default stylesheet. When set to empty, uses the default font collection from Google Fonts. A non-empty value replaces the family query string parameter in the Google Fonts URL. Header only."},
            {label=":iconfont-remote:",        detail="Allows use of a CDN for resolving the icon font. Only relevant used when value of icons attribute is font. Header only. Set by default."},
            {label=":imagesdir: ",             detail="Location of image files."},
            {label=":coderay-unavailable:",    detail="Instructs processor not to load CodeRay. Also set if processor fails to load CodeRay. Header only."},
            {label=":prewrap:",                detail="Wrap wide code listings. Set by default."},
            {label=":pygments-unavailable:",   detail="Instructs processor not to load Pygments. Also set if processor fails to load Pygments. Header only."},
            {label=":rouge-unavailable:",      detail="Instructs processor not to load Rouge. Also set if processor fails to load Rouge. Header only."},
            {label=":source-linenums-option:", detail="Turns on line numbers for source code listings."},
            {label=":copycss:",                detail="Copy CSS files to output. Only relevant when the linkcss attribute is set. Header only. Set by default."},
            {label=":linkcss:",                detail="Links to stylesheet instead of embedding it. Can’t be unset in SECURE mode. Header only."},
            {label=":stylesheet: ",            detail="CSS stylesheet file name. An empty value tells the converter to use the default stylesheet. Header only."},
            {label=":basebackend:",	           detail="The generic backend on which the backend is based. Typically derived from the backend value minus trailing numbers (e.g., the basebackend for docbook5 is docbook). May also indicate the internal backend employed by the converter (e.g., the basebackend for pdf is html)."},
            {label=":user-home:",	           detail="Get the full path of the home directory for the current user. Masked as . if the safe mode is SERVER or SECURE."},
            {label=":mantitle:",               detail="Metadata for man page output. Header only."},
            {label=":manvolnum:",              detail="Metadata for man page output. Header only."},
            {label=":manname:",                detail="Metadata for man page output. Header only."},
            {label=":manpurpose:",             detail="Metadata for man page output. Header only."},
            {label=":mansource:",              detail="Source (e.g., application and version) the man page describes. Header only."},
            {label=":manmanual:",              detail="Manual name displayed in the man page footer. Header only."},
            {label=":outfile:",                detail="Full path of the output file. (Cannot be referenced in the content. Only available to the API once the document is converted). Read-only."},
            {label=":outdir:",                 detail="Full path of the output directory. (Cannot be referenced in the content. Only available to the API once the document is converted)."},
            {label=":htmlsyntax:",             detail="Syntax used when generating the HTML output. Controlled by and derived from the backend name (html=html or xhtml=html)."},
            {label=":relfileprefix: ",         detail="The path prefix to add to relative xrefs."},
            {label=":highlightjsdir: ",        detail="Location of the highlight.js source code highlighter library. Default=CDN URL. Header only."},
            {label=":prettifydir: ",           detail="Location of non-CDN prettify source code highlighter library. Default=CDN URL. Header only."},
            {label=":stylesdir: ",             detail="Location of CSS stylesheets. Header only."},
            {label=":toc-class: ",             detail="CSS class on the table of contents container. Header only."},
            -- captions
            {label=":figure-caption!:",        detail="Disable figure captions."},
            {label=":table-caption!:",         detail='Disable table captions. If a title is given (.Title) then it will still show but there will be no label, i.e. "Table A.".'},
            -- options
            {label="%header",                  detail='Set first row of table as a header row. Implicitly set when first row is written on one line. Shorthand syntax for options="header".'},
            {label="source",                   detail='Mark open block as source code (verbatim).'},
            {label="python",                   detail='Mark source code as Python.'},
            {label="julia",                    detail='Mark source code as Julia.'},
            {label="indent=",                  detail='Indent text in block.'},
        }
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
            { name = 'asciidoc' },
            { name = 'path', option = {trailing_slash=true} },
            { name = 'nvim_lsp_signature_help' },
            { name = 'luasnip', options = {show_autosnippets=true} },
            { name = 'calc' },
            { name = 'buffer', group_index=2 },
            { name = 'dictionary', keyword_length=3, max_item_count=10, group_index=2 },
        }),
    })
end

return M
