; reclassify word "in" to the syntax group @repeat
"in" @repeat

; create new syntax group @flag that is a regex match for --flags etc
((word) @flag
  (#match? @flag "^[-+]{1,2}[a-zA-Z0-9_=-]+$"))

; reclassify "-f" in "[ -f ... ]" to @flag from @string
(test_operator) @flag

; find word OR concatenation nodes that look like paths
([(word) (concatenation)] @path
  (#match? @path "\~\?[a-zA-Z0-9_/:$@.*-]*[/.*][a-zA-Z0-9_/:$@.*-]*")
  ; don't consider a single dot a path, since it messes with sourcing (". ./script.sh")
  ; discovering the "not-eq" name was found by starting at https://tree-sitter.github.io/tree-sitter/using-parsers#an-example-program
  ; following Predicates link "WebAssembly binding" https://github.com/tree-sitter/tree-sitter/tree/master/lib/binding_web
  ; leading to a bindings.js
  (#not-eq? @path ".")
)

; same as above except in string and ~ is not allowed anymore
((string) @path
  (#match? @path "^\"[a-zA-Z0-9_/:$@.*-]*[/.*][a-zA-Z0-9_/:$@.*-]*\"$"))

((word) @operator
  (#lua-match? @operator "^+$"))
((word) @operator
  (#match? @operator "^then$"))

((command_name) @function.builtin
  (#match? @function.builtin "^\.$"))


(command
  name: (command_name) @hassubcmd
  .
  argument: (word) @function.call
  (#any-of? @hassubcmd
   ; builtins taken from zshPrecommand https://github.com/vim/vim/blob/master/runtime/syntax/zsh.vim
    "noglob" "nocorrect" "exec" "command" "builtin" "time" "sudo"
    ; extras
    "git" "pip" "pipx" "brew" "conda")
  ; only a subcommand if it's a word, not e.g. a flag
  (#match? @function.call "^[a-zA-Z][a-zA-Z0-9_-]*$"))

; in the testfile $0:h/script.sh wasn't captured as a function.call
; this was written by using TSPlaygroundToggle and seeing the tree for the uncaptured node
(command
  name: (command_name
          (concatenation
            (word) @function.call)))

