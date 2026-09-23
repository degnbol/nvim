#!/usr/bin/env python3
"""Tests for the stub-tree helpers. Run: ./test_stub_tree.py"""
import stub_tree as t


def test_comments_keyword_enum_member():
    src = "    None: typing.ClassVar[E]  # value = E.None\n"
    assert t.comment_keyword_targets(src).startswith("# ")


def test_keyword_commenting_leaves_docstring_prose_alone():
    # `from: <url>` inside a docstring reads as an annotated assignment.
    src = 'def f():\n    """\n    Plane of best fit\n    from: https://example.org\n    """\n'
    assert t.comment_keyword_targets(src) == src
    assert t.comment_keyword_targets("    None: int\n") == "#     None: int\n"


def test_side_effect_exclusion_matches_whole_path_components():
    assert t.side_effect_free(("DataStructs", "cDataStructs"))
    assert t.side_effect_free(("DataManip", "Metric"))
    assert not t.side_effect_free(("Chem", "Subshape", "demoCombined"))
    assert not t.side_effect_free(("sping", "PDF"))
    assert not t.side_effect_free(("__main__",))


if __name__ == "__main__":
    tests = [f for name, f in sorted(globals().items()) if name.startswith("test_")]
    for test in tests:
        test()
        print(f"ok {test.__name__}")
    print(f"{len(tests)} passed")
