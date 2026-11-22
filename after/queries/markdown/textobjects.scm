; extends

; A fence without a closing delimiter runs to the end of the document.
(fenced_code_block
  (code_fence_content) @fence.inner
  (fenced_code_block_delimiter) .) @fence.outer
