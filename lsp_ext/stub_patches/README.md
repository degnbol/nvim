Patches for `.pyi` stub trees, run where the runtime package is importable.
Kept out of `../python_stubs/`: `stubPath` resolves `.py` files too, ahead of
site-packages, so a helper module there would shadow a project's import of the
same name.

## The sequence

`patch_stubs.py <stub-dir-or-.pyi>` runs, on a staged copy:

1. the docstring pass (`docify_stubs.py`), on every package,
2. the package's repair, if `REPAIRS` names one (rdkit: `patch_pybind_stubs.py`),

then writes the changed files back, the marker file last. The marker is the
head line of the top-level `__init__.pyi` (or of a single-module stub):
`# patch_stubs: <package> <version> <digest>`, the digest covering every
script here and `requirements.txt`. Changing any of them, or the installed
version, re-runs the whole sequence. Exit status and invocation: the
`patch_stubs.py` docstring.

One sequence under one marker, because run independently the two passes undo
each other. The repair drops rdkit's C++-only constructor docstrings, which
docify reads as documentation and restores. docify gives property setters
their getter's docstring, which the repair's setter pattern has to accept.
`patch_stubs.py` alone decides which package a tree describes, whether it
applies, whether it is current and which repair runs.

## The docstring rule

A declaration whose runtime counterpart has documentation of its own gets it
written into the stub, on every overload. "Of its own" is
`docify.get_doc_def`: a module constant does not inherit its class's docstring
(`pi` does not get `float`'s). With `if_needed=True`, anything whose source
`inspect.getsourcefile` finds is left to basedpyright's own source fallback,
so in practice only compiled implementations are documented.

## docify's gaps

[docify](https://github.com/AThePeanut4/docify), pinned in `requirements.txt`,
writes basedpyright's docified typeshed. `docify_stubs.patch_docify` rebinds
two of its module globals to close three gaps on third-party trees:

- `get_qualname` raises on a member of a PEP 695 generic class, whose scope
  chain holds a libcst `AnnotationScope`.
- A class that exists only in the stub has no runtime object. `CLASS_ALIASES`
  in `patch_stubs.py` binds it on its module to the runtime class standing for
  it (numpy's `_ArrayOrScalarCommon` → `ndarray`).
- No visitor for annotated assignments (`divide: _UFunc_Nin2_Nout1[...]`).
  `docify_stubs.Transformer` adds a string literal on the following line.

## Tests

```zsh
./test_stub_tree.py
./test_patch_pybind_stubs.py
uv run --no-project --python 3.13 --with-requirements requirements.txt --with pytest \
    pytest -p no:cacheprovider --assert=plain test_docify_stubs.py test_patch_stubs.py
```

Add `--with rdkit` to the last to run the rdkit sequence test. `--assert=plain`
keeps pytest from rewriting rdkit's own `test*` modules, which the repair
introspects.
