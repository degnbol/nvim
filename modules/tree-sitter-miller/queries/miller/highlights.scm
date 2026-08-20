; Miller DSL highlights

; --- Keywords ---

; Anonymous keyword tokens inside named parent nodes
["if" "elif" "else" "for" "while" "do" "in"
 "begin" "end" "func" "subr" "call"
 "print" "printn" "eprint" "eprintn"
 "emit" "emit1" "emitf" "emitp"
 "dump" "edump" "tee"] @keyword

; Named statement nodes that are single keyword tokens (alias nodes)
(break_statement) @keyword
(continue_statement) @keyword
(return_statement "return" @keyword)
(filter_statement "filter" @keyword)
(unset_statement "unset" @keyword)

; Type keywords
(type_name) @type

; --- Literals ---

(boolean) @boolean
(string) @string
(string_content) @string
(escape_sequence) @string.escape
(number) @number
(comment) @comment

; --- References ---

; Field references ($name, $*, $[expr])
(field_ref "$" @punctuation.special)
(field_ref (identifier) @variable)
(field_ref "*" @operator)

; Out-of-stream variable references (@name, @*, @[expr])
(oosvar_ref "@" @punctuation.special)
(oosvar_ref (identifier) @variable)
(oosvar_ref "*" @operator)

; Special variables (NR, NF, FILENAME, etc.)
(special_variable) @variable.builtin

; --- Functions ---

; Builtin function calls — must come before the generic function.call rule
(function_call function_name: (identifier) @function.builtin
  (#any-of? @function.builtin
    "abs" "acos" "acosh" "antimode" "any" "append" "apply" "arrayify"
    "asin" "asinh"
    "asserting_absent" "asserting_array" "asserting_bool" "asserting_boolean"
    "asserting_empty" "asserting_empty_map" "asserting_error" "asserting_float"
    "asserting_int" "asserting_map" "asserting_nonempty_map" "asserting_not_array"
    "asserting_not_empty" "asserting_not_map" "asserting_not_null" "asserting_null"
    "asserting_numeric" "asserting_present" "asserting_string"
    "atan" "atan2" "atanh"
    "bitcount" "boolean"
    "capitalize" "cbrt" "ceil" "clean_whitespace" "collapse_whitespace"
    "concat" "contains" "cos" "cosh" "count"
    "depth" "dhms2fsec" "dhms2sec" "distinct_count"
    "erf" "erfc" "every" "exec" "exp" "expm1"
    "flatten" "float" "floor" "fmtifnum" "fmtnum" "fold" "format"
    "fsec2dhms" "fsec2hms"
    "get_keys" "get_values" "gmt2localtime" "gmt2nsec" "gmt2sec"
    "gssub" "gsub"
    "haskey" "hexfmt" "hms2fsec" "hms2sec" "hostname"
    "index" "int" "invqnorm"
    "is_absent" "is_array" "is_bool" "is_boolean" "is_empty" "is_empty_map"
    "is_error" "is_float" "is_int" "is_map" "is_nan" "is_nonempty_map"
    "is_not_array" "is_not_empty" "is_not_map" "is_not_null" "is_null"
    "is_numeric" "is_present" "is_string"
    "joink" "joinkv" "joinv" "json_parse" "json_stringify"
    "kurtosis"
    "latin1_to_utf8" "leafcount" "leftpad" "length"
    "localtime2gmt" "localtime2nsec" "localtime2sec"
    "log" "log10" "log1p" "logifit" "lstrip"
    "madd" "mapdiff" "mapexcept" "mapselect" "mapsum"
    "max" "maxlen" "md5" "mean" "meaneb" "median"
    "mexp" "min" "minlen" "mmul" "mode" "msub"
    "nsec2gmt" "nsec2gmtdate" "nsec2localdate" "nsec2localtime" "null_count"
    "os"
    "percentile" "percentiles" "pow"
    "qnorm"
    "reduce" "regextract" "regextract_or_else" "rightpad" "round" "roundm" "rstrip"
    "sec2dhms" "sec2gmt" "sec2gmtdate" "sec2hms" "sec2localdate" "sec2localtime"
    "select" "sgn" "sha1" "sha256" "sha512"
    "sin" "sinh" "skewness" "sort" "sort_collection"
    "splita" "splitax" "splitkv" "splitkvx" "splitnv" "splitnvx"
    "sqrt" "ssub" "stat" "stddev"
    "strfntime" "strfntime_local" "strftime" "strftime_local"
    "string" "strip" "strlen" "strmatch" "strmatchx"
    "strpntime" "strpntime_local" "strptime" "strptime_local"
    "sub" "substr" "substr0" "substr1"
    "sum" "sum2" "sum3" "sum4"
    "sysntime" "system" "systime" "systimeint"
    "tan" "tanh" "tolower" "toupper" "truncate" "typeof"
    "unflatten" "unformat" "unformatx"
    "upntime" "uptime" "urand" "urand32" "urandelement" "urandint" "urandrange"
    "utf8_to_latin1"
    "variance" "version"))

; User-defined function calls (not in builtin list)
(function_call function_name: (identifier) @function.call)

; Function/subroutine definitions
(func_definition name: (identifier) @function)
(subr_definition name: (identifier) @function)

; --- Parameters ---

(parameter (identifier) @variable.parameter)

; --- Operators ---

["+" "-" "*" "/" "//" "%" "**" "."
 ".+" ".-" ".*" "./"
 "==" "!=" "<" ">" "<=" ">=" "<=>"
 "=~" "!=~"
 "&&" "||" "!" "??" "???" "^^"
 "&" "|" "^" "~" "<<" ">>" ">>>"
 "?" ":"] @operator

; Assignment operators
["=" "+=" "-=" "*=" "/=" "//=" "%=" "**=" ".="] @operator

; Redirect operators (>, >>, |) in output statements
(redirect [">" ">>" "|"] @keyword.operator)

; --- Punctuation ---

["(" ")" "{" "}" "[" "]"] @punctuation.bracket
["," ";"] @punctuation.delimiter
