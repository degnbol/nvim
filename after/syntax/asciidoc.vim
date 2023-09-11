" Recognise hard line break
" https://docs.asciidoctor.org/asciidoc/latest/blocks/hard-line-breaks/
syn match linebreak ' +$'

" hide commentleader for simplicity
syn match commentDelimiter '^//'
" Override syntax groups from vim-asciidoctor to contain the commentDelimiter
" Also allow for italic and bold in comments.
syn match asciidoctorComment "^//.*$" contains=@Spell,commentDelimiter,asciidoctorBoldComment,asciidoctorItalicComment,asciidoctorBoldItalicComment
syn region asciidoctorComment start="^////.*$" end="^////.*$" contains=@Spell

" bold and italic in comments
syn region asciidoctorBoldComment matchgroup=Conceal start=/\m\*\*/ end=/\*\*/ contains=@Spell concealends contained
syn region asciidoctorBoldComment matchgroup=Conceal start=/\m\%(^\|[[:punct:][:space:]]\@<=\)\*\ze[^* ].\{-}\S/ end=/\*\%([[:punct:][:space:]]\@=\|$\)/ contains=@Spell concealends oneline contained
syn region asciidoctorItalicComment matchgroup=Conceal start=/\m__/ end=/__/ contains=@Spell concealends contained
syn region asciidoctorItalicComment matchgroup=Conceal start=/\m\%(^\|[[:punct:][:space:]]\@<=\)_\ze[^_ ].\{-}\S/ end=/_\%([[:punct:][:space:]]\@=\|$\)/ contains=@Spell concealends oneline contained
syn region asciidoctorBoldItalicComment matchgroup=Conceal start=/\m\*\*_/ end=/_\*\*/ contains=@Spell concealends contained
syn region asciidoctorBoldItalicComment matchgroup=Conceal start=/\m\%(^\|[[:punct:][:space:]]\@<=\)\*_\ze[^*_ ].\{-}\S/ end=/_\*\%([[:punct:][:space:]]\@=\|$\)/ contains=@Spell concealends oneline contained

syn region asciidoctorBold matchgroup=Conceal start=/\m\%(^\|[[:punct:][:space:]]\@<=\)\*\ze[^* ].\{-}\S/ end=/\*\%([[:punct:][:space:]]\@=\|$\)/ contains=@Spell concealends oneline
syn region asciidoctorItalic matchgroup=Conceal start=/\m\%(^\|[[:punct:][:space:]]\@<=\)_\ze[^_ ].\{-}\S/ end=/_\%([[:punct:][:space:]]\@=\|$\)/ contains=@Spell concealends oneline
syn region asciidoctorBoldItalic matchgroup=Conceal start=/\m\%(^\|[[:punct:][:space:]]\@<=\)\*_\ze[^*_ ].\{-}\S/ end=/_\*\%([[:punct:][:space:]]\@=\|$\)/ contains=@Spell concealends oneline

" disable syntax highlight for filenames and URLs inside comments
syn match filenameCommentNoSpell '[A-Za-z0-9-_]\+\.[a-z0-9]\+' contains=@NoSpell containedin=asciidoctorComment contained
syn match UrlCommentNoSpell '\w\+:\/\/[^[:space:]]\+' contains=@NoSpell containedin=asciidoctorComment contained

" highlight
syn region asciidocHighlight matchgroup=Conceal start=/##/ end=/##/ contains=@Spell concealends
syn region asciidocHighlight matchgroup=Conceal start=/\m\%(^\|[[:punct:][:space:]]\@<=\)#\ze[^# ].\{-}\S/ end=/#\%([[:punct:][:space:]]\@=\|$\)/ contains=@Spell concealends oneline

" underline. Placed after highlight so it trumps it even though their patterns 
" are very similar
syn region asciidocUnderline matchgroup=Conceal start=/\[underline\]##/ end=/##/ contains=@Spell concealends
syn region asciidocUnderline matchgroup=Conceal start=/\m\%(^\|[[:punct:][:space:]]\@<=\)\[underline\]#\ze[^# ].\{-}\S/ end=/#\%([[:punct:][:space:]]\@=\|$\)/ contains=@Spell concealends oneline
syn region asciidocBoldUnderline matchgroup=Conceal start=/\[underline\]\*\*/ end=/\*\*/ contains=@Spell concealends
syn region asciidocBoldUnderline matchgroup=Conceal start=/\m\%(^\|[[:punct:][:space:]]\@<=\)\[underline\]\*\ze[^* ].\{-}\S/ end=/\*\%([[:punct:][:space:]]\@=\|$\)/ contains=@Spell concealends oneline
syn region asciidocItalicUnderline matchgroup=Conceal start=/\[underline\]__/ end=/__/ contains=@Spell concealends
syn region asciidocItalicUnderline matchgroup=Conceal start=/\m\%(^\|[[:punct:][:space:]]\@<=\)\[underline\]_\ze[^_ ].\{-}\S/ end=/_\%([[:punct:][:space:]]\@=\|$\)/ contains=@Spell concealends oneline
" strikethrough
syn region asciidocStrikethrough matchgroup=Conceal start=/\[line-through\]##/ end=/##/ contains=@Spell concealends
syn region asciidocStrikethrough matchgroup=Conceal start=/\m\%(^\|[[:punct:][:space:]]\@<=\)\[line-through\]#\ze[^# ].\{-}\S/ end=/#\%([[:punct:][:space:]]\@=\|$\)/ contains=@Spell concealends oneline
syn region asciidocBoldStrikethrough matchgroup=Conceal start=/\[line-through\]\*\*/ end=/\*\*/ contains=@Spell concealends
syn region asciidocBoldStrikethrough matchgroup=Conceal start=/\m\%(^\|[[:punct:][:space:]]\@<=\)\[line-through\]\*\ze[^* ].\{-}\S/ end=/\*\%([[:punct:][:space:]]\@=\|$\)/ contains=@Spell concealends oneline
syn region asciidocItalicStrikethrough matchgroup=Conceal start=/\[line-through\]__/ end=/__/ contains=@Spell concealends
syn region asciidocItalicStrikethrough matchgroup=Conceal start=/\m\%(^\|[[:punct:][:space:]]\@<=\)\[line-through\]_\ze[^_ ].\{-}\S/ end=/_\%([[:punct:][:space:]]\@=\|$\)/ contains=@Spell concealends oneline

" passthrough
" https://docs.asciidoctor.org/asciidoc/latest/pass/pass-macro/
syn region asciidocPassthrough matchgroup=Conceal start=/++/ end=/++/ contains=@Spell concealends oneline
syn region asciidocPassthrough matchgroup=Conceal start=/\m\%(^\|[[:punct:][:space:]]\@<=\)+\ze[^+ ].\{-}\S/ end=/+\%([[:punct:][:space:]]\@=\|$\)/ contains=@Spell concealends oneline
" passthrough even has a triple version
syn region asciidocPassthrough matchgroup=Conceal start=/+++/ end=/+++/ contains=@Spell concealends oneline

" subscript and superscript.
" concealends removed since we don't have a way of showing it
syn region asciidocSubscript matchgroup=Conceal start=/\M~~/ end=/\M~~/ contains=@Spell oneline
syn region asciidocSuperscript matchgroup=Conceal start=/\^\^/ end=/\^\^/ contains=@Spell oneline
syn region asciidocSubscript matchgroup=Conceal start=/\m\%(^\|[[:punct:][:space:]]\@<=\)\~\ze[^\~ ].\{-}\S/ end=/\~\%([[:punct:][:space:]]\@=\|$\)/ contains=@Spell oneline
syn region asciidocSuperscript matchgroup=Conceal start=/\m\%(^\|[[:punct:][:space:]]\@<=\)\^\ze[^\^ ].\{-}\S/ end=/\^\%([[:punct:][:space:]]\@=\|$\)/ contains=@Spell oneline

