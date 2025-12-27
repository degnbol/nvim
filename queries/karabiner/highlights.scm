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
