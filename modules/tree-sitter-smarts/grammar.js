/// <reference types="tree-sitter-cli/dsl" />
// @ts-check

const symbols = require('./symbols.js');

/** @param {string[]} chars @returns {string} escaped body of a regex class */
function classBody(chars) {
  // `[` too: to a POSIX-class-aware engine `[[` opens a nested class.
  return [...new Set(chars)].join('').replace(/([\^\]\\\-\[])/g, '\\$1');
}

// Every character legal outside a bracket atom.
const bareChars = [
  ...symbols.bareAliphatic.join(''),
  ...symbols.bareAromatic.join(''),
  ...'Aa*.%()[]-=#$:/\\~@&,;!+<>',
  ...'0123456789',
];

module.exports = grammar({
  name: 'smarts',

  // The record layer is whitespace-significant, so every separator and newline
  // is threaded explicitly.
  extras: () => [],

  conflicts: $ => [
    // After an atom, `-` may open a ring closure (`C-1CCCCC-1`) or a chain bond
    // (`C-C`). LR(1) cannot decide.
    [$.branched_atom],
  ],

  rules: {
    source_file: $ => seq(
      optional($.header),
      repeat($._line),
      // A final line with no trailing newline.
      optional($._line_body),
    ),

    _line: $ => seq(optional($._line_body), $._eol),

    _line_body: $ => choice(
      seq(optional($._ws), $.comment),
      $.record,
      $._ws,
    ),

    // `C#C` is ethyne, so `#` is a triple bond mid-line; at line start no
    // structure may begin with a bond, which makes the comment unambiguous.
    comment: _ => token(seq('#', /[^\r\n]*/)),

    // The separator before the CX layer is optional because `ReactionFromSmarts`
    // reads `CC>>CC|$a;b;c;d$|` and applies the labels, while `MolFromSmiles`
    // rejects the same spelling without the space. Which of the two reads a
    // record is not something this rule can see, and accepting a line RDKit
    // rejects is the better error of the two. A name still needs its separator:
    // `|` cannot occur in a structure, but a name is free text, so an optional
    // separator would read every unparseable tail as one.
    record: $ => seq(
      $.structure,
      optional(seq(optional($._ws), $.cx_layer)),
      optional(seq($._ws, $.name)),
      optional($._ws),
    ),

    // The extended layer: a `|` pair around a comma-separated list of fields,
    // each read by RDKit's dispatch on its first character (CXSmilesOps.cpp,
    // `parse_it`). That dispatch skips whatever it does not recognise, so a
    // stray comma is accepted here too, but an unknown field is an error rather
    // than an annotation silently dropped. Nothing inside is whitespace
    // significant and text fields may hold spaces, so the `|` pair is the only
    // thing bounding the layer.
    cx_layer: $ => seq('|', repeat(choice(
      $.cx_coordinates,
      $.cx_atom_labels,
      $.cx_atom_values,
      $.cx_atom_props,
      $.cx_data_sgroup,
      $.cx_polymer_sgroup,
      $.cx_field,
      ',',
    )), '|'),

    // One point per atom in writing order. A component may be left empty — the
    // writer omits a zero z, as in `|(-1.29904,-0.25,;0,0.5,;1.29904,-0.25,)|`
    // — and an empty point puts its atom at the origin, so each separator
    // stands on its own rather than joining two numbers.
    cx_coordinates: $ => seq('(', repeat(choice($.coordinate, ',', ';')), ')'),

    coordinate: _ => /-?[0-9]+\.?[0-9]*|-?\.[0-9]+/,

    // Atom labels, one per atom, and the same block spelled `$_AV:` carries the
    // molfile values instead. Only `;` and `$` end an entry, so it may hold
    // anything else — spaces and `>` included.
    cx_atom_labels: $ => seq('$', optional($._cx_entries), '$'),
    cx_atom_values: $ => seq('$_AV:', optional($._cx_entries), '$'),

    // An entry may be empty, so what has to be there is either a label or a
    // separator; a rule that could match neither would match the empty string,
    // which cannot be generated.
    _cx_entries: $ => choice(
      $.label,
      seq(optional($.label), repeat1(seq(';', optional($.label)))),
    ),

    // A delimiter is written literally as a character reference — `&#44;` for a
    // comma — which read_text_to() resolves in every text field of the layer.
    // Only a label has to know: a reference ends in the `;` that separates two
    // labels, so a label that read one as plain text would split in two and every
    // atom after it would take the wrong one. Elsewhere a reference is text, as
    // no field is delimited by a character one can contain. Plain runs stop
    // before an `&`, or longest match would swallow the reference; a lone `&` is
    // itself a plain character.
    label: $ => repeat1(choice($._label_chars, $.escape, '&')),

    _label_chars: _ => /[^|;$&\r\n]+/,

    escape: _ => /&#[0-9]+;/,

    // `atomProp:index.name.value`, repeated with `:` between. A name ends at the
    // `.` before its value and the value at the `,` or `|` after it, so a value
    // may hold a `.` and neither may hold a `:`.
    cx_atom_props: $ => seq(
      'atomProp', ':', $.atom_prop, repeat(seq(':', $.atom_prop)),
    ),

    atom_prop: $ => seq($.index, '.', $.prop_name, '.', optional($.prop_value)),

    prop_name: _ => /[^|.:,\r\n]+/,
    prop_value: _ => /[^|:,\r\n]+/,

    // `SgD:[atoms]:[field name]:[data]:[query operator]:[unit]:[tag]:[coords]`,
    // the text fields commonly all empty (`|SgD:2,1:FIELD:info::::|`). The query
    // operator is where a bare `>` occurs outside a reaction.
    cx_data_sgroup: $ => seq(
      'SgD', ':', $._cx_index_list,
      repeat(seq(':', optional(choice($.sgroup_field, $.cx_coordinates)))),
    ),

    // `Sg:[type]:[atoms]:[subscript]:[superscript]:[head bonds]:[tail bonds]:`,
    // where a trailing field may be missing entirely rather than left empty.
    cx_polymer_sgroup: $ => seq(
      'Sg', ':', $.sgroup_type, ':', $._cx_index_list,
      optional(seq(':', optional($.sgroup_field),
        optional(seq(':', optional($.sgroup_field),
          optional(seq(':', optional($._cx_index_list),
            optional(seq(':', optional($._cx_index_list))))))))),
    ),

    // The keys of RDKit's sgroupTypemap; anything else is rejected there.
    sgroup_type: _ => choice(
      'n', 'mon', 'mer', 'co', 'xl', 'mod', 'mix', 'f', 'any', 'gen', 'c',
      'grf', 'alt', 'ran', 'blk',
    ),

    // RDKit tells the trailing coordinates of a data sgroup from a text field by
    // position, having read exactly five of the latter first. Here it is the
    // opening `(` that does it, which only rules out a text field starting with
    // one; a field may still hold parens anywhere else.
    //
    // The precedence keeps a one-letter text field a text field: `t` is also the
    // marker of the trans field, and a string literal beats a regex of the same
    // length. read_text_to() is greedy in the same way, so only a comma ends an
    // sgroup — `|SgD:0:x:y:::,SgH:0:1|` needs its separator.
    sgroup_field: _ => token(prec(1, /[^|:,(\r\n][^|:,\r\n]*/)),

    // Every other field names atoms or bonds by index: `,` separates entries,
    // `.` groups the indices of one, and `:` introduces the value an entry sets,
    // where `rb` and `s` may write `*` for "any". Binding rightward keeps a
    // trailing separator with the list — read_int_list() accepts one and the
    // writer emits it, as in `Sg:n:6,1,2,4::hh&#44;f:6,0,:4,2,`.
    cx_field: $ => seq($.cx_name, ':', $._cx_index_groups),

    _cx_index_groups: $ => prec.right(seq(
      $._cx_index_group,
      repeat(seq(',', $._cx_index_group)),
      optional(','),
    )),

    _cx_index_group: $ => seq(
      $._cx_index, repeat(seq(choice('.', ':'), $._cx_index)),
    ),

    _cx_index_list: $ => prec.right(seq(
      $.index, repeat(seq(',', $.index)), optional(','),
    )),

    _cx_index: $ => choice($.index, alias('*', $.wildcard)),

    index: _ => /[0-9]+/,

    // Dative and hydrogen bonds, zero-order bonds, unsaturation, substitution,
    // variable attachment, wedges, ring-bond count, link nodes, sgroup
    // hierarchy, double-bond stereo, radicals of 1-7 electrons, and the three
    // enhanced-stereo groups. A stereo group's number is part of its name
    // because it selects which group is named, as the electron count does.
    cx_name: _ => choice(
      'C', 'H', 'Z', 'u', 's', 'm', 'w', 'wU', 'wD', 'rb', 'LN', 'SgH',
      'ctu', 'c', 't', 'a', /o[0-9]+/, /&[0-9]+/, /\^[1-7]/,
    ),

    // May not begin with `|`: that opens the CX layer, and a `name` allowed to
    // start with one would outlex the marker by longest match, absorbing the
    // whole layer silently with no ERROR node.
    name: _ => /[^| \t\r\n][^\r\n]*/,

    // `SmilesMolSupplier`'s signature is `titleLine=True`, so the first line is
    // expected to be column names. Header and record are distinguished in the
    // lexer, by longest match, which is the only layer that can decide it: the
    // choice of token at the first character precedes any parse state, so a
    // dynamic precedence between the two rules has nothing to attach to
    // (verified — a broad `/[^\s]+/` header field swallows `CCO ethanol`
    // whole). The first field therefore has to be a token no structure can
    // match: its run up to the first `[` must hold a character that is not
    // legal outside a bracket atom, where the alphabet is only the bare element
    // symbols plus punctuation. `SMILES` qualifies on `M`; `CCO` does not.
    header: $ => seq(
      alias($._header_head, $.header_field),
      repeat(seq($._ws, alias($._header_word, $.header_field))),
      $._eol,
    ),

    _header_head: _ => token(new RegExp(
      `[^ \\t\\r\\n\\[]*[^ \\t\\r\\n\\[${classBody(bareChars)}][^ \\t\\r\\n]*`)),

    _header_word: _ => /[^ \t\r\n]+/,

    _eol: _ => /\r?\n/,

    // Space and tab are interchangeable as the field separator.
    _ws: _ => /[ \t]+/,

    // `MolFromSmarts` takes a lone structure, `ReactionFromSmarts` requires the
    // arrow; the grammar spans both. A reaction has exactly two `>`: one is "a
    // reaction requires at least one >>" and three is "multi-step reactions not
    // supported", both rejected. So `C>>C` states an empty agent side rather
    // than leaving it out, and every side may be empty — expressed here rather
    // than as an empty-matching `side` rule, which cannot be generated.
    structure: $ => choice(
      $.chain,
      seq(optional($.chain), '>', optional($.chain), '>', optional($.chain)),
    ),

    // A single optional recursive tail: `branched_atom (bond? chain)*` splits
    // `CCC` two ways and does not generate.
    chain: $ => seq(
      $.branched_atom,
      optional(seq(optional(choice($.bond, '.')), $.chain)),
    ),

    // `(` always follows an atom, so a component group leading a reaction side
    // (`(CC.O)>>CC`) is an ERROR node. RDKit reads it as one reactant and writes
    // it back out in the same form; fixing it means moving `.` out of `chain`.
    branched_atom: $ => seq($._atom, repeat($.ring_bond), repeat($.branch)),

    branch: $ => seq('(', optional($.bond), $.chain, ')'),

    _atom: $ => choice(
      $.organic_symbol,
      $.aromatic_symbol,
      $.wildcard,
      $.bracket_atom,
    ),

    // The marker only announces how many digits follow, so it is not part of the
    // number: one digit bare, exactly two after `%`, any number inside `%(...)`.
    // Three tokens under one name, differing only in length.
    ring_bond: $ => seq(optional($.bond), choice(
      alias($._ring_digit, $.ring_number),
      seq('%', alias($._ring_digits, $.ring_number)),
      seq('%', '(', alias($._ring_digits_long, $.ring_number), ')'),
    )),

    _ring_digit: _ => /[0-9]/,
    _ring_digits: _ => /[0-9][0-9]/,
    _ring_digits_long: _ => /[0-9]+/,

    bond: $ => seq(
      $.bond_prim,
      repeat(seq(optional(choice('&', ',', ';')), $.bond_prim)),
    ),

    // `!` is a prefix on any primitive, not a special `!@` token: `C!-!=C`
    // reads as `C!-&!=C`.
    //
    // `->` and `<-` are the dative bonds, one lexer token each in both dialects
    // (smiles.ll, smarts.ll) and primitives like any other, so they nest in a
    // bond expression (`[C]-,->[O]`) and close a ring (`N->1CC1`). A doubled
    // backslash is the same bond as a single one: the lexers read `[\\]{1,2}`.
    bond_prim: $ => seq(
      optional('!'),
      choice('-', '=', '#', '$', ':', '/', '\\', '\\\\', '~', '@', '->', '<-'),
    ),

    // Outside brackets the alphabet is a much smaller set than inside, and the
    // two must stay disjoint: `Sn` is sulfur + nitrogen, `[Sn]` is tin.
    organic_symbol: _ => choice(...symbols.bareAliphatic),
    aromatic_symbol: _ => choice(...symbols.bareAromatic),
    wildcard: _ => choice('*', 'A', 'a'),

    bracket_atom: $ => seq('[', $._atom_expr, ']'),

    // Named for their operator: SMARTS precedence is `&` > `,` > `;`, so the
    // outside-in nesting is low-AND → OR → high-AND. `&` is implicit when
    // primitives are juxtaposed.
    _atom_expr: $ => choice($._prim, $.and_low, $.or_expr, $.and_high),
    and_low: $ => prec.left(1, seq($._atom_expr, ';', $._atom_expr)),
    or_expr: $ => prec.left(2, seq($._atom_expr, ',', $._atom_expr)),
    and_high: $ => prec.left(3, seq($._atom_expr, optional('&'), $._atom_expr)),

    _prim: $ => choice(
      $.negation,
      $.isotope,
      $.element_symbol,
      alias($._bracket_aromatic, $.aromatic_symbol),
      $.hydrogen,
      $.wildcard,
      $.chirality,
      $.charge,
      $.atomic_number,
      $.property_primitive,
      $.recursive,
      $.atom_class,
    ),

    negation: $ => prec(4, seq('!', $._prim)),

    // Each of these is a token: otherwise the implicit `&` between primitives
    // makes `[--]` (charge −2 vs two negations), `[#6]` and `[C:1]` ambiguous.
    isotope: _ => /[0-9]+/,
    element_symbol: _ => choice(...symbols.bracketElement),
    _bracket_aromatic: _ => choice(...symbols.bracketAromatic),

    // The letter names the property and the digits count it, so the two are
    // separate. Both need the digits bound rightward, or they are equally well
    // read as an isotope — `[C13]` really is carbon of mass 13, so the choice is
    // the lexer's and RDKit binds them to the primitive too. `H` also has to
    // outrank the element symbol of the same spelling; `[H]` is element hydrogen
    // to RDKit but nodes here as an H-count with no digit.
    hydrogen: $ => prec.right(1, seq(
      alias('H', $.property_name), optional($.count))),

    chirality: _ => token(seq('@', optional(choice(
      '@',
      seq(choice('TH', 'AL', 'SP', 'TB', 'OH'), /[0-9]*/),
    )))),

    charge: _ => token(choice(/\+\+*/, /\+[0-9]+/, /--*/, /-[0-9]+/)),

    // Its own node, separate from `property_primitive`: `[#7]` is nitrogen, so
    // the element-colouring layer has to reach it.
    atomic_number: _ => token(seq('#', /[0-9]+/)),

    // `^` is the one letter whose digit is mandatory, and the one whose digit
    // names a hybridisation rather than counting anything.
    property_primitive: $ => prec.right(choice(
      seq(
        alias(new RegExp(`[${classBody(symbols.propertyPrimitives)}]`), $.property_name),
        optional($.count),
      ),
      seq(
        alias(new RegExp(`[${classBody(symbols.propertyPrimitivesWithDigit)}]`), $.property_name),
        $.count,
      ),
    )),

    count: _ => /[0-9]+/,

    // Nests a chain, not a structure: `[$(C>>C)]` is rejected by RDKit. `$` is
    // the recursion marker and the parens are ordinary grouping, so they are
    // three tokens. Safe despite `$` also being the quadruple bond, because
    // RDKit rejects a bond before a branch: `C$(CC)C` never occurs.
    recursive: $ => seq('$', '(', $.chain, ')'),

    // `:` prefixes the map number rather than being part of it, so the two are
    // separate nodes. Splitting is safe here although the other primitives must
    // each be one token: `:` starts nothing else inside a bracket atom, and the
    // digits are only reachable through it.
    atom_class: $ => seq(':', $.class_number),
    class_number: _ => /[0-9]+/,
  },
});
