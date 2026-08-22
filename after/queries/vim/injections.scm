; extends

; Highlight the shell command in `:!` filter statements (`:.!shfmt …`) as zsh,
; both in the ui2 cmdline and in .vim files. `:!` always runs via 'shell' (zsh
; here), so target zsh unconditionally — no bash/sh parser is installed anyway.
; include-children is mandatory: the command's args are named children
; (filter_command + command_argument), so without it only inter-arg whitespace
; is injected and the shell tree comes out empty.
((bang_filter_statement (command) @injection.content)
 (#set! injection.language "zsh")
 (#set! injection.include-children))
