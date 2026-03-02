/// <reference types="tree-sitter-cli/dsl" />
// @ts-check

// PyMOL atom selection algebra grammar
// Reference: https://pymolwiki.org/index.php/Selection_Algebra

const SELECTORS = [
  'name', 'resn', 'resi', 'chain', 'segi', 'alt', 'model',
  'index', 'id', 'rank', 'pepseq', 'label', 'elem', 'flag',
  'ss', 'rep', 'color', 'cartoon_color', 'ribbon_color',
  'numeric_type', 'formal_charge', 'partial_charge', 'state',
  'b', 'q', 'x', 'y', 'z',
];

const SELECTOR_SHORTHANDS = [
  'n.', 'r.', 'i.', 'c.', 's.', 'm.', 'idx.', 'ps.',
  'e.', 'f.', 'nt.', 'pc.', 'fc.', 'v.', 'pr.',
];

const BUILTINS = [
  'all', 'none', 'enabled', 'visible', 'bonded', 'protected',
  'fixed', 'restrained', 'masked', 'organic', 'inorganic',
  'solvent', 'polymer', 'guide', 'hetatm', 'hydrogens',
  'backbone', 'sidechain', 'metals', 'donors', 'acceptors',
  'present', 'center', 'origin',
];

const BUILTIN_SHORTHANDS = [
  'org.', 'ino.', 'sol.', 'pol.', 'h.', 'bb.', 'sc.',
  'don.', 'acc.', 'fxd.', 'rst.', 'msk.',
];

const EXPANSION_KEYWORDS = [
  'byres', 'bychain', 'bymolecule', 'byfragment', 'bysegi',
  'byobject', 'bycalpha', 'byring', 'bycell', 'bound_to',
  'neighbor', 'extend', 'first', 'last', 'in', 'like',
];

const EXPANSION_SHORTHANDS = [
  'br.', 'bc.', 'bm.', 'bf.', 'bs.', 'bca.', 'bto.',
  'nbr.', 'xt.',
];

const PROXIMITY_KEYWORDS = [
  'within', 'around', 'expand', 'gap', 'near_to', 'beyond',
];

const PROXIMITY_SHORTHANDS = [
  'w.', 'a.', 'x.', 'nto.', 'be.',
];

const REPRESENTATIONS = [
  'everything', 'lines', 'sticks', 'spheres', 'dots', 'surface', 'mesh',
  'cartoon', 'ribbon', 'labels', 'nb_spheres', 'nonbonded', 'volume',
  'slice', 'licorice', 'cell', 'extent',
];

// Case-insensitive keyword regex
function kw(word) {
  const escaped = word.replace(/\./g, '\\.');
  const pattern = escaped.split('').map(ch => {
    if (/[a-z]/i.test(ch)) return '[' + ch.toLowerCase() + ch.toUpperCase() + ']';
    return ch;
  }).join('');
  return new RegExp(pattern);
}

function kwChoice(words) {
  return choice(...words.map(w => kw(w)));
}

module.exports = grammar({
  name: 'pymol_select',

  extras: $ => [/\s/],

  word: $ => $.identifier,

  rules: {
    // Top-level allows bare operators for f-string fragments where
    // interpolation splits the string_content, e.g. f"{x} and polymer"
    // produces a fragment starting with "and".
    source: $ => repeat(choice(
      $.expression,
      $.logical_operator,
      $.not_operator,
    )),

    expression: $ => choice(
      $.paren_expr,
      $.not_expr,
      $.binary_expr,
      $.proximity_expr,
      $.expansion_expr,
      $.selector_expr,
      $.builtin,
      $.representation,
      $.macro,
      $.object_ref,
      // No bare identifier — tightens the grammar so non-pymol strings
      // (like "hello world") produce ERROR instead of parsing as atoms.
    ),

    paren_expr: $ => seq(
      '(',
      $.expression,
      ')',
    ),

    binary_expr: $ => prec.left(1, seq(
      $.expression,
      $.logical_operator,
      $.expression,
    )),

    not_expr: $ => prec(2, seq(
      $.not_operator,
      $.expression,
    )),

    not_operator: _ => choice(/[nN][oO][tT]/, '!'),

    logical_operator: _ => choice(
      /[aA][nN][dD]/,
      /[oO][rR]/,
      '&', '|',
    ),

    proximity_expr: $ => prec(2, seq(
      $.proximity_keyword,
      $.number,
      optional(seq($.of, $.expression)),
    )),

    proximity_keyword: _ => choice(
      kwChoice(PROXIMITY_KEYWORDS),
      kwChoice(PROXIMITY_SHORTHANDS),
    ),

    expansion_expr: $ => prec(2, seq(
      $.expansion_keyword,
      $.expression,
    )),

    expansion_keyword: _ => choice(
      kwChoice(EXPANSION_KEYWORDS),
      kwChoice(EXPANSION_SHORTHANDS),
    ),

    of: _ => /[oO][fF]/,

    selector_expr: $ => prec(3, seq(
      $.selector,
      $.value,
    )),

    selector: _ => choice(
      kwChoice(SELECTORS),
      kwChoice(SELECTOR_SHORTHANDS),
    ),

    value: $ => prec.right(repeat1(choice(
      $.comparison_operator,
      // Multi-part values like "A+B+C" or "1-100" as a single token.
      // Avoids shift/reduce conflicts between value continuation and
      // new expression when identifiers overlap with selector keywords.
      $.multi_value,
      $.identifier,
      $.number,
      $.wildcard,
    ))),

    // Matches values joined by + or - without spaces: A+B+C, 1-100, C+N+O
    // prec(2) ensures this wins over _macro_segment at the same length.
    multi_value: _ => token(prec(2, /[A-Za-z0-9_*]+[+\-][A-Za-z0-9_*]+([+\-][A-Za-z0-9_*]+)*/)),

    builtin: _ => choice(
      kwChoice(BUILTINS),
      kwChoice(BUILTIN_SHORTHANDS),
      /[pP][oO][lL][yY][mM][eE][rR]\.[pP][rR][oO][tT][eE][iI][nN]/,
      /[pP][oO][lL][yY][mM][eE][rR]\.[nN][uU][cC][lL][eE][iI][cC]/,
    ),

    representation: _ => kwChoice(REPRESENTATIONS),

    // Macro-style selection: /object/segment/chain/resi/atom
    // A single token captures the entire macro path to avoid parser-level
    // conflicts between _macro_segment and identifier.
    // Three regex forms reject file paths while accepting pymol macros:
    macro: _ => token(prec(4, choice(
      // 1a. Leading //: //A/100/, ///A/
      /\/\/[A-Za-z0-9_*+\-]*(?:\/[A-Za-z0-9_*+\-]*)*/,
      // 1b. Internal //: /obj//A/100/CA — has filled segments before the //
      /\/[A-Za-z0-9_*+\-]+(?:\/[A-Za-z0-9_*+\-]+)*\/\/[A-Za-z0-9_*+\-]*(?:\/[A-Za-z0-9_*+\-]*)*/,
      // 2. Leading / with all 5 segments filled: /pept/seg/A/100/CA
      //    Exactly 4 internal slashes. 5-segment absolute paths are vanishingly rare.
      /\/[A-Za-z0-9_*+\-]+\/[A-Za-z0-9_*+\-]+\/[A-Za-z0-9_*+\-]+\/[A-Za-z0-9_*+\-]+\/[A-Za-z0-9_*+\-]+/,
      // 3. Trailing form (no leading /): A/100/CA, 100/CA
      //    Right-to-left pymol macro. Limited to 1-4 slashes.
      /[A-Za-z0-9_*+\-]+(?:\/[A-Za-z0-9_*+\-]*){1,4}/,
    ))),

    // Explicit object reference with prefix (%, ?)
    object_ref: $ => seq($.object_prefix, $.identifier),

    identifier: _ => /[A-Za-z_][A-Za-z0-9_]*/,
    number: _ => token(prec(1, /\d+(\.\d+)?/)),
    wildcard: _ => '*',
    comparison_operator: _ => choice('<', '>', '=', '<=', '>='),
    object_prefix: _ => choice('%', '?'),
  },
});
