-- Compile flags for C/C++ files without a compile command: include paths of libraries
-- outside the compiler's default search path, e.g. Homebrew's include/glib-2.0/.
-- Without them, every `gboolean f(...)` is an invalid declaration, dropped from the
-- AST, and goto-definition on f finds nothing.
return require("utils/pkg_config").cflags({ "glib-2.0", "hdf5" })
