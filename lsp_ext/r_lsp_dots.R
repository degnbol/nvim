# Monkey-patch R languageserver to resolve ... in wrapper functions.
# When a function has only `...` as formals (e.g. scale_y_log10),
# trace the body to find the function it delegates to and return
# that function's formals instead. Also fixes parameter documentation
# by pointing the resolve handler at the underlying function.

resolve_dots = function(fn, depth = 3L) {
  if (depth <= 0L || is.primitive(fn) || is.null(body(fn))) return(NULL)
  b = body(fn)
  exprs = if (identical(b[[1]], as.name("{"))) as.list(b)[-1] else list(b)
  dots = as.name("...")
  env = environment(fn)
  for (expr in exprs) {
    if (!is.call(expr)) next
    if (!any(vapply(as.list(expr)[-1], identical, logical(1), dots))) next
    fname = expr[[1]]
    if (!is.name(fname)) next
    target = tryCatch(get(as.character(fname), envir = env), error = function(e) NULL)
    if (!is.function(target)) next
    fmls = formals(target)
    if (!identical(names(fmls), "...")) {
      return(list(formals = fmls, funct = as.character(fname)))
    }
    resolved = resolve_dots(target, depth - 1L)
    if (!is.null(resolved)) return(resolved)
  }
  NULL
}

# Find the resolved function name for a ... wrapper (cached per session)
.dots_cache = new.env(parent = emptyenv())

resolve_dots_cached = function(funct, package) {
  key = paste0(package, "::", funct)
  if (exists(key, envir = .dots_cache)) return(get(key, envir = .dots_cache))
  fn = tryCatch(get(funct, envir = asNamespace(package)), error = function(e) NULL)
  if (!is.function(fn)) { assign(key, NULL, envir = .dots_cache); return(NULL) }
  fmls = formals(if (is.primitive(fn)) args(fn) else fn)
  if (length(fmls) == 1L && identical(names(fmls), "...")) {
    result = resolve_dots(fn)
  } else {
    result = NULL
  }
  assign(key, result, envir = .dots_cache)
  result
}

setHook(packageEvent("languageserver", "onLoad"), function(...) {
  lsp_ns = asNamespace("languageserver")
  PkgNS = get("PackageNamespace", envir = lsp_ns)

  # Patch get_formals to resolve ... wrappers
  PkgNS$set("public", "get_formals", function(funct, exported_only = TRUE) {
    if (!self$exists_funct(funct, exported_only = exported_only)) return(NULL)
    ns = asNamespace(self$package_name)
    fn = get(funct, envir = ns)
    fmls = formals(if (is.primitive(fn)) args(fn) else fn)
    if (length(fmls) == 1L && identical(names(fmls), "...")) {
      resolved = resolve_dots(fn)
      if (!is.null(resolved)) return(resolved$formals)
    }
    fmls
  }, overwrite = TRUE)

  # Wrap arg_completion: add textEdit for " = " insertion and fix
  # data$funct for resolved ... wrappers so docs resolve correctly.
  # textEdit is more reliable than insertText with some completion clients.
  orig_arg_completion = get("arg_completion", envir = lsp_ns)
  patched_arg_completion = function(uri, workspace, point, token, funct,
                                     package = NULL, exported_only = TRUE) {
    items = orig_arg_completion(uri, workspace, point, token, funct,
                                package, exported_only)
    if (length(items)) {
      # Replace insertText with textEdit for reliable " = " insertion
      row = point$row     # 0-indexed
      col = point$col     # 0-indexed, cursor position
      start_col = col - nchar(token)
      for (i in seq_along(items)) {
        items[[i]]$textEdit = list(
          range = list(
            start = list(line = row, character = start_col),
            end   = list(line = row, character = col)
          ),
          newText = items[[i]]$insertText
        )
        items[[i]]$insertText = NULL
      }
      # Fix data$funct for resolved ... wrappers
      if (!is.null(package)) {
        resolved = resolve_dots_cached(funct, package)
        if (!is.null(resolved)) {
          for (i in seq_along(items)) {
            items[[i]]$data$funct = resolved$funct
          }
        }
      }
    }
    items
  }
  assignInNamespace("arg_completion", patched_arg_completion, ns = "languageserver")
})
