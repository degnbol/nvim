; Minimal highlights for testing — full highlights in Phase 2

; Keywords (anonymous tokens inside named parent nodes)
["if" "elif" "else" "for" "while" "do" "in"
 "begin" "end" "func" "subr" "call"
 "print" "printn" "eprint" "eprintn"
 "emit" "emit1" "emitf" "emitp"
 "dump" "edump" "tee"] @keyword

; Named statement nodes that are just keywords
(break_statement) @keyword
(continue_statement) @keyword
(return_statement "return" @keyword)
(filter_statement "filter" @keyword)
(unset_statement "unset" @keyword)

; Type keywords
(type_name) @type

; Booleans
(boolean) @boolean

; Literals
(string) @string
(string_content) @string
(escape_sequence) @string.escape
(number) @number
(comment) @comment

; Field references
(field_ref "$" @punctuation.special)
(field_ref (identifier) @variable)
(field_ref "*" @operator)

; OOS variable references
(oosvar_ref "@" @punctuation.special)
(oosvar_ref (identifier) @variable)
(oosvar_ref "*" @operator)

; Special variables
(special_variable) @variable.builtin

; Function calls
(function_call function_name: (identifier) @function.call)

; Function/subroutine definitions
(func_definition name: (identifier) @function)
(subr_definition name: (identifier) @function)

; Operators
["+" "-" "*" "/" "//" "%" "**" "."
 ".+" ".-" ".*" "./"
 "==" "!=" "<" ">" "<=" ">=" "<=>"
 "=~" "!=~"
 "&&" "||" "!" "??" "???" "^^"
 "&" "|" "^" "~" "<<" ">>" ">>>"
 "=" "+=" "-=" "*=" "/=" "//=" "%=" "**=" ".="
 "?" ":"] @operator

; Punctuation
["(" ")" "{" "}" "[" "]"] @punctuation.bracket
["," ";"] @punctuation.delimiter
