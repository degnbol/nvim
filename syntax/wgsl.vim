" Bevy/naga_oil preprocessor directives (not in the wgsl treesitter grammar)
syntax match @keyword.import /^#import\>/
syntax match @keyword.import /^#define_import_path\>/
syntax match @keyword.directive /^#ifdef\>/
syntax match @keyword.directive /^#ifndef\>/
syntax match @keyword.directive /^#else\>/
syntax match @keyword.directive /^#endif\>/
