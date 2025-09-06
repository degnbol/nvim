syn keyword KeydKey C M A S f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12 esc escape plus minus equal backspace tab leftbrace rightbrace enter control leftcontrol semicolon apostrophe grave leftshift backslash comma dot slash rightshift alt leftalt space capslock numlock scrolllock rightcontrol rightalt home up pageup left right end down pagedown insert delete macro brightnessup brightnessdown scale dashboard mute volumedown volumeup pause playpause previoussong nextsong leftmeta rightmeta meta compose
syn keyword Error ctrl
syn keyword KeydOption overload_tap_timeout
syn keyword @title.builtin ids global main containedin=Title
syn keyword @function.builtin oneshot overload macro layer swap swapm
hi def link KeydKey Constant
hi def link KeydOption Identifier
hi def @title.builtin gui=underdouble,italic
syn match Delimiter ','
syn match Delimiter '\.'
syn match Delimiter '-'
syn match Operator '*'
syn match Operator '='
syn match Comment '^#.*'
syn region Title matchgroup=Delimiter start='\[' end='\]' contains=Variable
syn region Arguments matchgroup=Delimiter start='(' end=')' contains=ALL
