"""Write runtime docstrings into a stub tree with docify.

docify (https://github.com/AThePeanut4/docify) produces basedpyright's docified
typeshed. On a third-party tree it has three gaps: it raises on a member of a
PEP 695 generic class (closed by ``patch_docify``), it cannot reflect a class
that exists only in the stub (closed by ``document``'s aliases), and it has no
visitor for annotated assignments (closed by this module's ``Transformer``).
"""
from __future__ import annotations

import importlib
import importlib.metadata
import inspect
import pathlib
import re
import tokenize
import types
import warnings
from collections.abc import Callable, Mapping, Sequence

import docify
import libcst as cst
import libcst.matchers as m
import libcst.metadata as meta

from stub_tree import comment_keyword_targets

_REQUIREMENTS = pathlib.Path(__file__).with_name("requirements.txt")
_ANNOTATED_NAME = m.SimpleStatementLine([m.AnnAssign(target=m.Name())])
_DOCSTRING = m.SimpleStatementLine([m.Expr(m.SimpleString()), m.ZeroOrMore()])


def patch_docify() -> None:
    """Make docify's ``Transformer`` use this module's ``get_qualname``.

    Idempotent. It rebinds ``docify.get_qualname``, a docify internal, so it is
    only trusted on the version ``requirements.txt`` pins.

    Raises:
        SystemExit: ``requirements.txt`` pins no docify, or the installed one is
            another version.
    """
    pin = re.search(r"(?m)^docify==(\S+)", _REQUIREMENTS.read_text())
    if not pin:
        raise SystemExit(f"no docify pin in {_REQUIREMENTS}")
    pinned = pin.group(1)
    installed = importlib.metadata.version("docify")
    if installed != pinned:
        raise SystemExit(f"docify {installed} installed, {pinned} pinned")
    docify.get_qualname = get_qualname


def get_qualname(scope: meta.Scope, name: str) -> str:
    """The ``__qualname__`` Python gives a declaration: libcst's
    ``Scope.get_qualified_names_for``, where docify's raises on any scope but a
    global or class one.

    Args:
        scope: the scope a declaration's name is bound in.
        name: the declared name.

    Returns:
        The dotted name of the declaration below its module.
    """
    # docify.get_obj returns None on the AttributeError a `<locals>` part
    # raises, so a function-local name gets no docstring.
    (qualname,) = {q.name for q in scope.get_qualified_names_for(name)
                   if q.source is meta.QualifiedNameSource.LOCAL}
    return qualname


class Transformer(docify.Transformer):
    """docify's transformer, also documenting annotated assignments.

    A string literal goes on the line after an assignment to a bare name whose
    runtime value has documentation of its own, by ``docify.get_doc_def``'s rule.
    """

    def _assignment_doc(self, statement: cst.SimpleStatementLine) -> str | None:
        """Cleaned documentation of an annotated assignment's runtime value.

        Args:
            statement: an original-tree line holding one annotated assignment to a name.

        Returns:
            The documentation, or None where there is none of its own, the
            declaration is unreachable at this Python version, or
            ``inspect.getsourcefile`` finds the value's source.
        """
        if self.get_metadata(docify.UnreachableProvider, statement, False):
            return None
        target = statement.body[0].target
        scope = self.get_metadata(meta.ScopeProvider, target, None)
        if scope is None:
            return None
        qualname = get_qualname(scope, target.value)
        found = docify.get_obj(self.mod, qualname)
        if found is None:
            return None
        scope_obj, obj = found
        if not self.check_if_needed(obj):
            return None
        doc = docify.get_doc_def(scope_obj, obj, qualname, target.value)
        return inspect.cleandoc(doc) or None if isinstance(doc, str) else None

    def _block_indent(self, block: cst.CSTNode) -> str:
        """Indentation of a block's statements, summed over its enclosing blocks.

        Args:
            block: an original-tree module or indented block.

        Returns:
            The whitespace prefixing each of its statements.
        """
        indent = ""
        node = block
        while node is not None:
            if isinstance(node, cst.IndentedBlock):
                indent += node.indent if node.indent is not None else self.module.default_indent
            node = self.get_metadata(meta.ParentNodeProvider, node, None)
        return indent

    def _document_assignments(
        self,
        original: cst.Module | cst.IndentedBlock,
        updated: cst.Module | cst.IndentedBlock,
    ) -> list[cst.BaseStatement]:
        """A block's statements with a docstring after each documented annotated
        assignment that lacks one.

        Args:
            original: the block as parsed, for metadata lookups.
            updated: the same block with its children already transformed, which
                keeps the statement count of the original.

        Returns:
            The updated block's statements, docstrings inserted.
        """
        indent = self._block_indent(original)
        statements = []
        for i, (before, after) in enumerate(zip(original.body, updated.body)):
            statements.append(after)
            if not m.matches(before, _ANNOTATED_NAME):
                continue
            if i + 1 < len(original.body) and m.matches(original.body[i + 1], _DOCSTRING):
                continue
            doc = self._assignment_doc(before)
            if doc:
                literal = cst.SimpleString(docify.docquote_str(doc, indent))
                statements.append(cst.SimpleStatementLine([cst.Expr(literal)]))
        return statements

    # At block level, since whether a string literal already follows an
    # assignment is a sibling, invisible from the statement itself.
    def leave_IndentedBlock(self, original_node, updated_node):
        return updated_node.with_changes(
            body=self._document_assignments(original_node, updated_node))

    def leave_Module(self, original_node, updated_node):
        # Before docify's own, which may prepend a module docstring and so break
        # the pairing of original and updated statements.
        updated_node = updated_node.with_changes(
            body=self._document_assignments(original_node, updated_node))
        return super().leave_Module(original_node, updated_node)


def document(
    modules: Sequence[tuple[str, pathlib.Path]],
    safe: Callable[[Sequence[str]], bool],
    aliases: Mapping[tuple[str, str], str],
) -> dict[str, BaseException]:
    """Document stub files in place from their imported runtime modules.

    Every stub that tokenizes has its keyword-named assignment targets commented
    out, including the stubs of modules not safe to import. Calls
    ``patch_docify``, binds each alias as an attribute of its imported module,
    and ignores the warnings importing and reading the runtime raise.

    Args:
        modules: (dotted module name, stub path) pairs, parents first.
        safe: whether a module is safe to import, given the components of its
            dotted name below the package root.
        aliases: (module, stub-only class) -> name of the runtime class on that
            module standing for it.

    Returns:
        {module: error} for each module skipped. A stub that does not tokenize
        is left untouched. One whose module fails to import, or that libcst
        cannot parse, gets its keyword targets commented and nothing else.
    """
    patch_docify()
    skipped = {}
    # Importing and reading the runtime warns (deprecations, mostly), about
    # packages nothing here can act on.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for name, path in modules:
            try:
                # libcst rejects a keyword target (`None: T`).
                text = comment_keyword_targets(path.read_text(encoding="utf-8"))
            except (SyntaxError, tokenize.TokenError) as e:
                skipped[name] = e
                continue
            path.write_text(text, encoding="utf-8")
            if not safe(name.split(".")[1:]):
                continue
            try:
                module = importlib.import_module(name)
            except (Exception, SystemExit) as e:
                skipped[name] = e
                continue
            for (owner, stub_class), runtime_class in aliases.items():
                if owner == name:
                    setattr(module, stub_class, getattr(module, runtime_class))
            try:
                path.write_text(with_docstrings(text, name, module), encoding="utf-8")
            except cst.ParserSyntaxError as e:
                skipped[name] = e
    return skipped


def with_docstrings(text: str, name: str, module: types.ModuleType) -> str:
    """A stub's source with docstrings written in from its runtime module.

    Only what ``inspect.getsourcefile`` cannot find the source of is documented
    (docify's ``if_needed``).

    Args:
        text: stub source.
        name: dotted name of the module the stub describes.
        module: that module, imported.

    Returns:
        The documented source.

    Raises:
        libcst.ParserSyntaxError: libcst cannot parse text.
    """
    tree = cst.parse_module(text)
    return cst.MetadataWrapper(tree).visit(Transformer(name, module, True)).code
