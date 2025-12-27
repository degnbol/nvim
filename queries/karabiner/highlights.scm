; Karabiner-specific property names highlighted as @property.builtin
; Applied on top of base json highlights via ftplugin/karabiner.lua
((pair
  key: (string
    (string_content) @_key) @property.builtin)
  (#any-of? @_key
    "global" "show_in_menu_bar"
    "rules" "description" "manipulators" "type" "from" "to"
    "key_code" "modifiers" "mandatory" "optional" "conditions" "is_keyboard"
    "frontmost_application_if" "frontmost_application_unless"
    "device_if" "device_unless" "keyboard_type_if" "keyboard_type_unless"
    "input_source_if" "input_sources" "input_source_id" "input_source_unless" "variable_if" "variable_unless"
    "event_changed_if" "event_changed_unless"
    "bundle_identifiers" "file_paths" "identifiers" "vendor_id" "product_id"
    "set_variable" "name" "value" "shell_command" "select_input_source"
    "mouse_key" "pointing_button" "consumer_key_code"
    "profiles" "selected" "parameters" "complex_modifications"
    "simple_modifications" "fn_function_keys" "devices"
    "virtual_hid_keyboard" "country_code" "keyboard_type_v2")
  (#set! priority 95))

; Builtin values for "type" key
((pair
  key: (string (string_content) @_key)
  value: (string (string_content) @string.builtin))
  (#eq? @_key "type")
  (#any-of? @string.builtin
    "basic" "mouse_motion_to_scroll"
    ; condition types
    "frontmost_application_if" "frontmost_application_unless"
    "device_if" "device_unless" "keyboard_type_if" "keyboard_type_unless"
    "input_source_if" "input_source_unless" "variable_if" "variable_unless"
    "event_changed_if" "event_changed_unless")
  (#set! priority 95))

; Builtin modifier values in arrays
((pair
  key: (string (string_content) @_key)
  value: (array (string (string_content) @string.builtin)))
  (#any-of? @_key "mandatory" "optional")
  (#any-of? @string.builtin
    "any" "caps_lock" "command" "control" "fn" "left_command" "left_control"
    "left_option" "left_shift" "option" "right_command" "right_control"
    "right_option" "right_shift" "shift")
  (#set! priority 95))

; Builtin key_code values
((pair
  key: (string (string_content) @_key)
  value: (string (string_content) @string.builtin))
  (#eq? @_key "key_code")
  (#any-of? @string.builtin
    ; letters
    "a" "b" "c" "d" "e" "f" "g" "h" "i" "j" "k" "l" "m"
    "n" "o" "p" "q" "r" "s" "t" "u" "v" "w" "x" "y" "z"
    ; numbers
    "0" "1" "2" "3" "4" "5" "6" "7" "8" "9"
    ; function keys
    "f1" "f2" "f3" "f4" "f5" "f6" "f7" "f8" "f9" "f10" "f11" "f12"
    "f13" "f14" "f15" "f16" "f17" "f18" "f19" "f20"
    ; modifiers
    "caps_lock" "left_control" "left_shift" "left_option" "left_command"
    "right_control" "right_shift" "right_option" "right_command" "fn"
    ; special keys
    "return_or_enter" "escape" "delete_or_backspace" "delete_forward"
    "tab" "spacebar" "hyphen" "equal_sign" "open_bracket" "close_bracket"
    "backslash" "non_us_pound" "semicolon" "quote" "grave_accent_and_tilde"
    "comma" "period" "slash" "non_us_backslash"
    ; arrows and navigation
    "up_arrow" "down_arrow" "left_arrow" "right_arrow"
    "page_up" "page_down" "home" "end"
    ; keypad
    "keypad_num_lock" "keypad_slash" "keypad_asterisk" "keypad_hyphen"
    "keypad_plus" "keypad_enter" "keypad_period" "keypad_equal_sign"
    "keypad_0" "keypad_1" "keypad_2" "keypad_3" "keypad_4"
    "keypad_5" "keypad_6" "keypad_7" "keypad_8" "keypad_9"
    ; misc
    "insert" "print_screen" "scroll_lock" "pause"
    "menu" "application" "help" "power" "execute" "select" "stop" "again" "undo"
    "cut" "copy" "paste" "find" "mute" "volume_increment" "volume_decrement"
    "locking_caps_lock" "locking_num_lock" "locking_scroll_lock"
    "vk_none")
  (#set! priority 95))

; Builtin pointing_button values
((pair
  key: (string (string_content) @_key)
  value: (string (string_content) @string.builtin))
  (#eq? @_key "pointing_button")
  (#any-of? @string.builtin
    "button1" "button2" "button3" "button4" "button5"
    "button6" "button7" "button8" "button9" "button10")
  (#set! priority 95))

; Builtin consumer_key_code values
((pair
  key: (string (string_content) @_key)
  value: (string (string_content) @string.builtin))
  (#eq? @_key "consumer_key_code")
  (#any-of? @string.builtin
    "display_brightness_decrement" "display_brightness_increment"
    "dictation" "rewind" "play_or_pause" "fast_forward" "mute"
    "volume_decrement" "volume_increment" "eject"
    "al_terminal_lock_or_screensaver" "al_calculator"
    "ac_search" "ac_home" "ac_back" "ac_forward" "ac_stop" "ac_refresh"
    "ac_bookmarks" "ac_next_keyboard_layout_select" "ac_desktop_show_all_windows"
    "apple_display_brightness_decrement" "apple_display_brightness_increment"
    "apple_top_case_display_brightness_decrement" "apple_top_case_display_brightness_increment"
    "menu" "menu_pick" "menu_up" "menu_down" "menu_left" "menu_right" "menu_escape")
  (#set! priority 95))
