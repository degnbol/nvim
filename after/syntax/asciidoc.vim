" hide commentleader for simplicity
syn match commentDelimiter '^//'
" Override syntax groups from vim-asciidoctor to contain the commentDelimiter
" Also allow for italic and bold in comments.
syn match asciidoctorComment "^//.*$" contains=@Spell,commentDelimiter,asciidoctorBoldComment,asciidoctorItalicComment,asciidoctorBoldItalicComment
syn region asciidoctorComment start="^////.*$" end="^////.*$" contains=@Spell

" bold and italic in comments
syn region asciidoctorBoldComment matchgroup=Conceal start=/\m\*\*/ end=/\*\*/ contains=@Spell concealends oneline contained
syn region asciidoctorBoldComment matchgroup=Conceal start=/\m\%(^\|[[:punct:][:space:]]\@<=\)\*\ze[^* ].\{-}\S/ end=/\*\%([[:punct:][:space:]]\@=\|$\)/ contains=@Spell concealends oneline contained
syn region asciidoctorItalicComment matchgroup=Conceal start=/\m__/ end=/__/ contains=@Spell concealends oneline contained
syn region asciidoctorItalicComment matchgroup=Conceal start=/\m\%(^\|[[:punct:][:space:]]\@<=\)_\ze[^_ ].\{-}\S/ end=/_\%([[:punct:][:space:]]\@=\|$\)/ contains=@Spell concealends oneline contained
syn region asciidoctorBoldItalicComment matchgroup=Conceal start=/\m\*\*_/ end=/_\*\*/ contains=@Spell concealends oneline contained
syn region asciidoctorBoldItalicComment matchgroup=Conceal start=/\m\%(^\|[[:punct:][:space:]]\@<=\)\*_\ze[^*_ ].\{-}\S/ end=/_\*\%([[:punct:][:space:]]\@=\|$\)/ contains=@Spell concealends oneline contained
syn match asciidoctorBoldComment /\%(^\|[[:punct:][:space:]]\@<=\)\*[^* ].\{-}\S\*\%([[:punct:][:space:]]\@=\|$\)/ contains=@Spell contained
" single char *b* bold
syn match asciidoctorBoldComment /\%(^\|[[:punct:][:space:]]\@<=\)\*[^* ]\*\%([[:punct:][:space:]]\@=\|$\)/ contains=@Spell contained
syn match asciidoctorBoldComment /\*\*\S.\{-}\*\*/ contains=@Spell contained
syn match asciidoctorItalicComment /\%(^\|[[:punct:][:space:]]\@<=\)_[^_ ].\{-}\S_\%([[:punct:][:space:]]\@=\|$\)/ contains=@Spell contained
" single char _b_ italic
syn match asciidoctorItalicComment /\%(^\|[[:punct:][:space:]]\@<=\)_[^_ ]_\%([[:punct:][:space:]]\@=\|$\)/ contains=@Spell contained
syn match asciidoctorItalicComment /__\S.\{-}__/ contains=@Spell contained
syn match asciidoctorBoldItalicComment /\%(^\|[[:punct:][:space:]]\@<=\)\*_[^*_ ].\{-}\S_\*\%([[:punct:][:space:]]\@=\|$\)/ contains=@Spell contained
" single char *_b_* bold+italic
syn match asciidoctorBoldItalicComment /\%(^\|[[:punct:][:space:]]\@<=\)\*_[^*_ ]_\*\%([[:punct:][:space:]]\@=\|$\)/ contains=@Spell contained
syn match asciidoctorBoldItalicComment /\*\*_\S.\{-}_\*\*/ contains=@Spell contained
" these needed as well for some reason
syn match formatDelim '_' conceal contained containedin=asciidoctorItalicComment,asciidoctorBoldItalicComment
syn match formatDelim '\*' conceal contained containedin=asciidoctorBoldComment,asciidoctorBoldItalicComment

" disable syntax highlight for filenames and URLs inside comments
syn match filenameCommentNoSpell '[A-Za-z0-9-_]\+\.[a-z0-9]\+' contains=@NoSpell containedin=asciidoctorComment contained
syn match UrlCommentNoSpell '\w\+:\/\/[^[:space:]]\+' contains=@NoSpell containedin=asciidoctorComment contained
