; Each emphasis_delimiter node is one `*` or `_`.
(emphasis) @emphasis.outer

((emphasis) @emphasis.inner
  (#offset! @emphasis.inner 0 1 0 -1))

(strong_emphasis) @emphasis.outer

((strong_emphasis) @emphasis.inner
  (#offset! @emphasis.inner 0 2 0 -2))
