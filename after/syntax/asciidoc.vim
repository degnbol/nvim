" Let's hide the comment leader by matching on it alone.
syntax match commentDelimiter '//'
" redefine to contain the delim
syntax match asciidoctorComment '//.*' contains=commentDelimiter

