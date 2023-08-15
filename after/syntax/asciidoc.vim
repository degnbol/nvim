" hide commentleader for simplicity
syn match commentDelimiter '//'
" Override syntax groups from vim-asciidoctor to contain the commentDelimiter
syn match asciidoctorComment "^//.*$" contains=@Spell,commentDelimiter
syn region asciidoctorComment start="^////.*$" end="^////.*$" contains=@Spell
