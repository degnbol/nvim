; extends
; The key of `["b"] = 2` is `["b"]`, as the lhs of `t["b"] = 2` is `t["b"]`.
(field
  "["? @assignment.lhs
  name: (_) @assignment.lhs
  "]"? @assignment.lhs
  value: (_) @assignment.rhs)
