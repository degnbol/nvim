# PLAN: close the runtime-vs-stub gaps in rdkit's type stubs

## Goal

Three classes of defect remain in `lsp_ext/python_stubs/rdkit/` after
`patch_pybind_stubs.py`'s current passes, all of them cases where the stub says
less than the live interpreter knows.

Every count below was measured with `basedpyright --outputjson` at
`typeCheckingMode = "standard"` against rdkit 2026.03.5, with the vendored tree
as the only stub source (no installed `rdkit-stubs` in the probe environment —
its presence shifts the totals by a couple of names, so pin the probe when
re-measuring). Probe environment and scripts: `~/.cache/rdkit_lsp_probe`.

| defect | scale |
| --- | --- |
| module-level names present at runtime, absent from the stub → false `reportAttributeAccessIssue` | `Chem.AllChem` 534/793 attributes, `Chem.Descriptors` 174/420, `Chem.Lipinski` 24/35, `Chem.Crippen` 2/8 |
| pure-Python modules with no `.pyi` → `"unknown import symbol"` | 75, of which 9 are outside the demo/test/`Contrib` exclusions below and 4 are library modules: `Chem.PandasTools`, `Chem.PandasPatcher`, `Chem.Draw.IPythonConsole`, `Chem.Draw.InteractiveRenderer` |
| methods monkey-patched in Python over their Boost original → stub keeps the stale C++ signature | ~25 unique `(class, method)` pairs, of which 8 are non-dunder |

`Chem`, `Chem.rdchem`, `Chem.Draw` and `Chem.rdMolDescriptors` measure 0
unresolved attributes, so the first defect is not tree-wide: it hits modules
whose names arrive by a mechanism pybind11-stubgen does not record.

Passes A and B are new passes in `patch_pybind_stubs.py`, so they reach both the
vendored tree and the `--in-place` rewrite of an environment's own
`rdkit-stubs` (`lua/autocmds/stub_patches.lua`). Pass C belongs in `RUNME.sh`
instead — see there for why. Name them for the operation:
`synthesise_missing_names`, `reexport_missing_names`, `type_boost_overrides`;
the fixer's docstring enumerates three defects today and becomes five.

## Pass A — name top-up from the runtime

`AllChem.AddHs`, `Descriptors.MolWt`, `Crippen.MolLogP` and ~730 siblings are
errors today. Two different origins, and the distinction decides the fix:

- **519 of AllChem's 534** are objects owned by another module whose stub already
  declares them with a real signature — `AllChem.AddHs` *is*
  `rdmolops.AddHs`, and `rdmolops.pyi` types it. AllChem's `.py` gets them by
  `from rdkit.Chem import *`; the stub lists ~180 explicit imports and drops the
  rest.
- **The descriptor modules** build their names as module-level lambdas in loops
  (`Lipinski.NOCount`, `GraphDescriptors.Chi0n`, `Fragments.fr_Al_COO` — 85 of
  them, `MolSurf` 38, `EState_VSA` 10). No module's stub declares these, so
  there is nothing to point at.

### The guard both phases need

Skip any name **bound at the stub's own top level** — a declaration, a plain
`from x import y`, or an `import x`. Star imports are not followed.

This guard is what makes the pass safe, and it is not optional for either
phase. Appending at end of file makes the last binding win, and a name's
runtime `__module__` is sometimes a *worse* home than the declaration already in
the file. Measured with the guard removed (1343 aliases appended across `Chem`,
`AllChem`, `Draw`, `Descriptors`): 694 errors fixed, 0 new — and 4 regressions,
all from `rdkit.Chem.inchi`, the Python wrapper that owns names the stub
already declared better:

```
rdkit.Chem.InchiToInchiKey  was (inchi: str) -> str  now (inchi: Unknown) -> Unknown
rdkit.Chem.MolToInchiKey    was (mol, options: str = '') -> str  now -> Unknown
```

With the guard, over 6 modules: 696 fixed, 0 new, 0 regressions.

### Classifying a name

`getattr(obj, "__module__")` alone misclassifies non-functions, in two measured
ways: 97 names tree-wide are submodule objects with no `__module__` at all, and
an *instance* inherits `__module__` from its class, so 5 names
(`MACCSkeys.GenMACCSKeys`, `Pairs.GetAtomPairFingerprintAsIntVect`, …) name an
origin module that has no such attribute. Aliasing those produces a stub file
carrying its own `unknown import symbol` error while the call site silently
degrades to `Unknown` — the gap masked rather than fixed.

So classify by identity against the claimed origin, in three branches:

1. `isinstance(obj, types.ModuleType)` → **re-export the submodule**:
   `from <parent> import <name> as <name>`. This is the mechanism behind the
   `from rdkit.Chem import PandasTools` failure, so the ones that still do not
   resolve are exactly pass C's targets. Verify on one no-stub submodule before
   building the branch: an alias to a module pyright cannot resolve may just
   relocate the error.
2. `origin` is set and `getattr(sys.modules[origin], name, None) is obj` →
   **re-export**: `from <origin> import <name> as <name>`.
3. otherwise → **synthesise locally** from the runtime object.

Phase ordering: synthesis must run over the whole tree before any re-export is
emitted, since branch 2's target may be a name phase 1 has yet to write. That
means two passes over the tree — read, synthesise, write, then re-read — which
is the one place these passes do not fit `patch_dir`'s single
read-modify-write loop. The ordering does not affect consumers: an alias whose
target is missing degrades to `Unknown` at the call site and errors only inside
the stub file.

### Synthesising a declaration

From `inspect.signature`, with **every default rendered as `= ...`**, never the
runtime repr. `str(inspect.signature(obj))` is unusable verbatim: of 214
candidates, 102 do not parse as Python, e.g.

```
def fr_Al_COO(mol, countUnique=True, pattern=<rdkit.Chem.rdchem.Mol object at 0x1058f1ee0>): ...
```

which is both a syntax error and address-dependent, so the vendored tree would
differ run to run. `= ...` accepts any argument type while keeping arity
checking (measured: `Lipinski.NOCount()` still reports "Argument missing for
parameter x"). Append the runtime docstring for hover.

Parameter names survive for most candidates (`NOCount(x)`,
`fr_Al_COO(mol, countUnique = ..., pattern = ...)`) but not all: `MolWt`,
`MolLogP` and `MolMR` are `lambda *x, **y`, so they yield no parameter
information at all. Return types are not recoverable without calling, so they
stay unannotated.

Synthesis cannot shadow a better type *today* — a probe touching all 190
guard-surviving candidates resolves 0 of them — but that is a property of this
rdkit version, not a consequence of the guard. Assert the 0 in the end-to-end
test rather than reasoning from the rule.

### Importing the tree is the hazard

The current fixer imports only modules carrying untyped properties, on purpose.
Pass A must import a module to enumerate its names, and rdkit's tree contains
demos that compute and print on import. Exclude by **path component** —
`Contrib`, `sping`, `tests`, `TestRunner`, `conftest`, `demo*` — not by
substring on the dotted name, which would silently drop `rdkit.DataStructs` and
`rdkit.DataManip` (11 stubs; `DataStructs` alone measures 63 phase-2 gains). A
`Data` exclusion is pointless either way: `rdkit/Data/` has no `__init__.py`, so
it is never walked.

Under those exclusions only `rdkit.ML.Cluster.murtagh_test` prints on import, so
redirecting `sys.stdout` to devnull across the import (keeping the script's own
report on a saved handle) suffices. Seven modules raise on import; route them
through the existing `skipped` reporting path.

## Pass B — methods patched in Python

`rdkit/Chem/__init__.py` runs `rdchem.Mol.GetAtoms = lambda self:
_GetAtomsIterator(self)` at import. The stub keeps the C++-derived
`-> typing.Iterable[Atom]`, which is wrong in both directions: it loses `len()`
and `[i]` (both hard errors today) and names a type the call never returns.

### Detection: drive from the stub

For each method a stub file declares, rewrite it when the runtime attribute on
that class is a `types.FunctionType` **and** the class itself is defined in a C
extension module (`sys.modules[cls.__module__].__file__` is not a `.py`). The
two halves solve different problems:

- **Driving from the stub deduplicates.** The same `Mol` is reachable as
  `Chem.Mol`, `AllChem.Mol`, `PandasPatcher.Mol` and `rdchem.Mol`; a walk of
  runtime classes yields ~46 hits for ~25 unique pairs, and would rewrite the
  class in every re-exporting module's stub. A stub-driven pass visits each
  declaration once.
- **The extension-module test excludes pure-Python subclasses of Boost
  classes.** `FeatMaps.FeatMapPoint` subclasses a Boost class, so *every* method
  is a `FunctionType` and overrides nothing. Nothing about the stub side
  distinguishes these; only where the class is defined does.

Do **not** filter on whether a Boost base carries the attribute — that rejects
the primary case. `Mol` *is* the Boost class, not a subclass of one:

```
Mol mro: ['rdkit.Chem.rdchem.Mol', 'Boost.Python.instance', 'builtins.object']
boost bases have GetAtoms?: [('instance', False)]
```

The metaclass does not discriminate either (`type(Mol)` and
`type(FeatMapPoint)` are both `Boost.Python.class`).

Scope: `__iter__`, `__next__` and `__getitem__` are in — 14 `rdBase` vector
classes and `MatchTypeVect` patch `__iter__`, and typing it gives `for x in
vect` a real element type. Lifecycle dunders (`__init__`, `__reduce__`,
`__getstate__`, `__setstate__`, `__ge__`) are out: a stale signature there is
harmless.

### Typing the rewrite

For a patched method taking nothing but `self`, the return type is readable off
a live call on a no-arg instance — the same trick `property_types` already uses,
and `Mol()` is constructible. Methods needing arguments keep their annotations
(see Ceiling).

`_GetAtomsIterator` is itself a stubgen stub of a pure-Python class, so every
method is unannotated and `atoms[0]` yields `Unknown`. Its generic shape cannot
be introspected, so it needs a literal override: a table from
`(module, class name)` to replacement class text **plus a module-level
preamble**, applied as a class-block substitution and asserted to have matched,
so an rdkit rename fails loudly instead of silently dropping the fix.

```python
class _GetRDKitObjIterator(typing.Generic[_RDKitObj]):
    def __init__(self, mol: Mol) -> None: ...
    def __iter__(self) -> _GetRDKitObjIterator[_RDKitObj]: ...
    def __next__(self) -> _RDKitObj: ...
    def __getitem__(self, i: int) -> _RDKitObj: ...
    def __len__(self) -> int: ...
class _GetAtomsIterator(_GetRDKitObjIterator[Atom]): ...
class _GetBondsIterator(_GetRDKitObjIterator[Bond]): ...
```

The preamble is `import typing` and `_RDKitObj = typing.TypeVar('_RDKitObj')`;
`Chem/__init__.pyi` has neither today. It is the one part of pass B that is not
naturally idempotent, so insert each line only when absent — a second
`--in-place` run must leave the tree byte-identical.

`Atom`/`Bond` are in scope via the existing `from .rdchem import *`, and
`rdchem.pyi` referring back to `rdkit.Chem._GetAtomsIterator` needs no new
import (`Chem/rdchem.pyi:6` already has `import rdkit.Chem`); pyright resolves
stub import cycles. Verified:

```
Type of "atoms" is "_GetAtomsIterator"
Type of "atoms[0]" is "Atom"
Type of "len(atoms)" is "int"
Type of "[a.GetSymbol() for a in atoms]" is "list[str]"
```

A one-line `-> typing.Sequence[Atom]` is the cheaper alternative: same
`len`/index/iterate, no new classes, no preamble — but it hides the real type
name on hover and over-promises `.count`/`.index`, which the runtime object
lacks.

## Pass C — generate the missing module stubs, in `RUNME.sh`

`from rdkit.Chem import PandasTools` errors because rdkit's stub build
environment has no pandas, and pyright never falls back to source for a module
missing from a stub package. Two non-options, both tested: deleting a `.pyi`
makes the module *unresolvable* rather than source-inferred, and a `partial\n`
`py.typed` marker has no effect on a `stubPath` tree (PEP 561 partial stubs
apply to `<pkg>-stubs` in site-packages — see Ceiling for what it does do
there).

So the gap has to be filled, with `stubgen --inspect-mode` — the generator the
other packages in `RUNME.sh` already use. This belongs in `RUNME.sh`, not the
fixer: the fixer's contract is to *fix* a package's bundled stubs and
`stub_patches.lua:10` promises nothing in it is package-specific, whereas
generation needs a build dependency (mypy) and only ever runs for the vendored
tree. As a loop in `RUNME.sh` feeding stub-less modules to `stubgen`, with the
output merged into the tree before `patch_pybind_stubs.py --out` runs, the
`--out`-only gate and the report wording disappear.

`RUNME.sh:44` currently bypasses the `stub()` helper
(`uv run --no-project --with rdkit python3 patch_pybind_stubs.py rdkit --out rdkit`),
so that environment has **no mypy, no pandas, no IPython**. All three matter:

- `PandasTools` imports fine without pandas but defines nothing — `LoadSDF` sits
  inside the `else:` of a `try: import pandas` (PandasTools.py:235) — so stubgen
  would emit an *empty* stub for the headline target. Worse than none, because
  it looks fixed.
- Of the 9 non-excluded stub-less modules, 6 cannot import in a bare
  environment: `Draw.IPythonConsole`, `Draw.InteractiveRenderer` (IPython),
  `DSViewer` (win32com), `MolDb.Loader_sa` (sqlalchemy), `utils.chemdraw` (the
  ChemDraw application), `Pharm2D.LazyGenerator` (raises `NotImplementedError`).

Add `--with mypy --with pandas --with ipython`, and let the loop tolerate and
report per-module import failure rather than aborting.

## Re-patching already-fixed environments

`stub_patches.lua` treats any tree whose head line starts with
`# fix_pybind_stubs:` as done, so environments patched by today's script would
never pick up these passes. Extend the marker with a fingerprint of the fixer
itself — `# fix_pybind_stubs: rdkit 2026.3.5 <sha256[:8] of patch_pybind_stubs.py>`
(the version is `importlib.metadata.version`'s PEP 440 form, not
`rdkit.__version__`'s `2026.03.5`) — and have the Lua side compare against
`vim.fn.sha256` of the script it already holds a path to. Verified equal to
`hashlib.sha256(path.read_bytes()).hexdigest()[:8]`; `util.readtext` reads the
whole file in binary, so there is no newline mismatch.

Two implementation details: the check at `stub_patches.lua:120` becomes a prefix
test *plus* a parse of the head line's last field, because the Lua side cannot
reconstruct the environment's rdkit version; and the marker format is
documented in three places — the fixer docstring, `README.md`, and the
`StubFix.marker` annotation — all of which move together.

The cost of fingerprinting the script is that any edit to it re-prompts every
environment. The fingerprint also covers only the script, so a tree patched in
an environment that later gains pandas stays "done".

## Tests

`test_patch_pybind_stubs.py` covers the text transforms on inline stub snippets;
passes A and B are runtime-driven, so they need a module to introspect. Build a
throwaway module in `tmp_path` — a Boost class cannot be faked, but both the
`FunctionType`-on-a-class detection and the signature synthesis work on a plain
class — and assert:

- a lambda-built name is synthesised at its origin, once, with its parameters
  and `= ...` for every default;
- a name the origin module does not actually hold (`__module__` lying) is not
  aliased;
- a name bound at the stub's top level by a plain `from x import y` is left
  alone — the `Chem.inchi` regression above, as a regression test;
- a submodule object is re-exported, not synthesised;
- the class-block override raises when its target class is absent, and its
  preamble is inserted exactly once across two runs;
- two runs leave the tree byte-identical (verified for guarded phases 1+2 over
  `Lipinski`/`AllChem`/`Descriptors`: 720 appends, then 0). This is the property
  the whole `--in-place` re-patch story rests on.

Add a slow end-to-end check that runs basedpyright over a probe file and asserts
the revealed types quoted above, plus the 0-shadowing count. It is the only test
that catches a pyright behaviour change, and the harness these measurements came
from.

## Ceiling

- **417 methods return `typing.Any`** and stay that way. Only 32 still carry a
  C++ signature line, and those are unmappable
  `boost::python::objects::iterator_range<…>` templates; the other 385 were
  opaque at the C++ level too (`_object*`, `boost::python::object`) —
  `GetSubstructMatches`, `Conformer.GetPositions`, `RingInfo.AtomRings`. Nothing
  static can recover them; only a curated table would.
- **Patched methods taking arguments** (`FilterMatcher.*`, `PropertyMol.SetProp`)
  keep their stale annotations. Typing them means calling them, which needs
  fixture arguments per method — a curated table again.
- **Descriptor return types** stay unannotated. They are recoverable by calling
  each entry of `Descriptors.descList` on a probe molecule, which is a bigger
  step than pass A's introspection: it executes library code for ~200 names.
  Worth revisiting separately.
- **`--in-place` stays weaker than the vendored tree** by exactly pass C, so a
  project owning a `pyrightconfig.json` keeps the `PandasTools` import error.
  Writing a `partial\n` `py.typed` into the *installed* `rdkit-stubs` would fix
  that with no mypy and no stubgen — measured: `PandasTools.LoadSDF` resolves
  with its full parameter list while the Boost stubs still win where they exist.
  Rejected, because a partial-marked installed package then overrides
  `stubPath` for every project sharing that environment, and the vendored tree's
  property fixes go with it (`SmilesParserParams().allowCXSMILES` measured
  `bool` → `Unknown`). Trading those across an environment for one import is the
  wrong swap, but the option is a one-line reversal if that judgement changes.
