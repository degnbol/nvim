syn keyword Key alt shift meta control C S tab esc insert leftalt rightalt leftshift rightshift capslock insert
syn keyword @function.builtin oneshot overload macro
hi def link Key Constant
syn match Delimiter ','
syn match Delimiter '\.'
syn match Delimiter '-'
syn match Operator '*'
syn match Operator '='
syn match Comment '^#.*'
syn region Title matchgroup=Delimiter start='\[' end='\]' contains=Variable
syn region Arguments matchgroup=Delimiter start='(' end=')' contains=ALL
