= Glossaries/acronyms/abbreviations

There's several packages for this.
- glossarium:0.5.9
 - featured package.
 - has `@Ref` to capitalize `@ref` which is nice and I haven't seen it in another pkg.
 - LSP completes things like `@ref:pl` which might be nice.
 - Seems to not have the features it promises on the github yet, i.e.:
  - `@ref:longplural` or similar
  - `@ref:description`
 - Also doesn't work:
  - `@ref:long:pl`
 - No support for articles, i.e. "a" vs "an"
- glossy:0.9.0
 - Has emphasis on `:` modifiers, very nice.
 - Works with chaining modifiers.
 - Downside: no `@Ref`, has to do `@ref:cap`
 - Support for articles, i.e. "a" vs "an"
- Many others that don't have the `@` notation...

Current best option: glossy

