# DRAFT PLAN: Julia function-definition form toggle

Toggle the function definition **under the cursor** between Julia's two forms:

```julia
f(x) = x^2                    ⟷    function f(x)
                                       x^2
                                   end
```

Replaces the two register macros in `ftplugin/julia.lua` (`@i` long→short,
`@f` short→long) — brittle keystroke sequences that can't handle
`where`/return-type signatures, multi-statement bodies, or overload sets.

## Grammar (verified against `julia.so` via `get_string_parser`)

- **Long form** `function f(x) … end` → `function_definition` with children
  `signature` (wrapping the callee) and `block` (the body).
- **Short form** `f(x) = …` → `assignment`, then `operator "="`, then the body
  expression. The LHS distinguishes a short function def from a plain variable
  assignment. LHS node shapes:

  | source | `assignment` LHS node |
  |---|---|
  | `f(x) = …` | `call_expression` |
  | `f(x)::Int = …` | `typed_expression` (first child `call_expression`) |
  | `f(x) where T = …` | `where_expression` (first child `call_expression`) |
  | `a = 1` | `identifier` → **must not match** |

- **Macro prefix** `@inline f(x) = x` → `macrocall_expression >
  macro_argument_list > assignment`. The def is nested; the macro is a parent.
- **Docstring** `"""…"""` above a `function` → the `string_literal` is a
  **preceding sibling** of `function_definition`, outside its range.

## Module: `lua/julia_fndef.lua`

### `toggle()` — the primitive (cursor def only)

1. Node under cursor (`vim.treesitter.get_node()`); ascend to the **nearest**
   `function_definition`, or `assignment` whose LHS is a `call_expression` /
   `typed_expression` / `where_expression` (first child `call_expression`).
   Stop at that node — do **not** climb to an enclosing `macrocall_expression`
   (preserves `@inline`); the rewrite range excludes any preceding docstring
   sibling automatically. None found → `vim.notify` "no function definition
   under cursor", return.
2. Branch on node type:
   - **short → long**: `sig` = full LHS text (everything before the top-level
     `operator "="`), `body` = text after it. `where T` / `::Int` live on the
     LHS, so they carry over verbatim and `function f(x)::Int … end` /
     `function f(x) where T … end` are valid. Replace node range with
     `function <sig>\n    <body>\nend`.
   - **long → short**: only if `block` has a **single named child** (one
     expression). `<signature> = <expr>` (the `signature` text already wraps any
     `::Int`/`where T`). Multi-statement body → notify "body is not a single
     expression", return (no inline form exists).
3. `nvim_buf_set_text` over the node's range → single undo step. Indentation is
   computed **deterministically**: the emitted body gets the def line's leading
   whitespace + one `shiftwidth` (honouring `expandtab`), `end` gets the def
   line's leading whitespace. `=`/`indentexpr` is *not* used — the buffer's
   `indentexpr` (`vim.treesitter.indentexpr()`, set in `after/indent/julia.vim`)
   is not a core function and silently no-ops, so a `=` pass is unreliable.
   `vim.lsp.buf.format` is likewise avoided (depends on julials' formatter).

### `conform()` — coerce all same-name methods to the cursor def's current form

Semantics **(A) conform-to-cursor**: the cursor def is the unchanged exemplar;
every same-name method in the buffer is coerced to *its current* form.
Workflow: `toggle()` the too-long overload to long form, then `conform()` so the
sibling overloads follow.

1. Resolve the cursor def (same ascent as `toggle()`). Its current form is the
   target. Extract the **callee text** at the head of its `call_expression` —
   a bare `identifier` (`f`) or qualified `field_expression` (`Base.show`).
2. Treesitter-query the buffer for all defs (both forms) whose signature callee
   text matches exactly (so `foo` ≠ `Base.foo`). Scope = current buffer only.
3. Coerce each non-matching-form def to the target form. Long→short can't inline
   a multi-statement body → collect those, skip them, and `vim.notify` which
   names/lines were skipped (don't fail the whole op).
4. Apply edits **bottom-to-top** (or via extmarks) so ranges stay valid; wrap
   the batch in one undo block.

## Keymaps (`ftplugin/julia.lua`, buffer-local)

- `n <localleader>f` → `require("julia_fndef").toggle()`.
- `n <localleader>F` → `require("julia_fndef").conform()` (uppercase = wider
  scope: the whole overload set).

Remove the `@i`/`@f` `setreg` lines (and the explanatory comment above them)
once these land.

## Out of scope / YAGNI

- **Visual range** — superseded by name-based `conform()`, which keys on *what*
  the methods are, not *where* they sit. No selection edge cases.
- **Cross-file / cross-module conform** — Julia methods are module-keyed;
  buffer-scope is enough. Caveat accepted: two unrelated same-name defs in
  separate `module` blocks in one file would both match — restrict to the
  enclosing module only if it bites.
- Anonymous functions (`x -> …`), `do` blocks — different syntax, separate concern.
- No whole-buffer "convert *all* defs" mode — there is no deterministic canonical
  form; conform is scoped to one name at a time.
