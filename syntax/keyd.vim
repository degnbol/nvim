syn keyword Key C M A S esc escape plus minus equal backspace tab leftbrace rightbrace enter leftcontrol semicolon apostrophe grave leftshift backslash comman dot slash rightshift leftalt space capslock numlock scrolllock rightcontrol rightalt home up pageup left right end down pagedown insert delete macro mute volumedown volumeup pause leftmeta rightmeta compose
syn keyword @function.builtin oneshot overload macro layer
hi def link Key Constant
syn match Delimiter ','
syn match Delimiter '\.'
syn match Delimiter '-'
syn match Operator '*'
syn match Operator '='
syn match Comment '^#.*'
syn region Title matchgroup=Delimiter start='\[' end='\]' contains=Variable
syn region Arguments matchgroup=Delimiter start='(' end=')' contains=ALL
