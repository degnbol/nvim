/// <reference types="tree-sitter-cli/dsl" />
// @ts-check

// Tree-sitter grammar for Miller DSL (the language used in put/filter/tee verbs)
// Reference: https://miller.readthedocs.io/en/latest/reference-dsl/

// Precedence levels (higher = tighter binding)
const PREC = {
  ASSIGN: 1,
  TERNARY: 2,
  NULL_COALESCE: 3,
  OR: 4,
  XOR: 5,
  AND: 6,
  BITOR: 7,
  BITXOR: 8,
  BITAND: 9,
  EQUALITY: 10,
  COMPARE: 11,
  SHIFT: 12,
  ADD: 13,
  MUL: 14,
  UNARY: 15,
  POWER: 16,
  POSTFIX: 17,
  CALL: 18,
};

module.exports = grammar({
  name: 'miller',

  extras: $ => [/\s/, $.comment],

  word: $ => $.identifier,

  conflicts: $ => [],

  rules: {
    source: $ => repeat(choice(
      $._statement,
      // Catch-all for injection contexts (shell quotes, etc.)
      $._fallback,
    )),

    _statement: $ => seq(
      choice(
        $.begin_block,
        $.end_block,
        $.func_definition,
        $.subr_definition,
        $.assignment,
        $.if_statement,
        $.for_in_statement,
        $.for_c_statement,
        $.while_statement,
        $.do_while_statement,
        $.return_statement,
        $.break_statement,
        $.continue_statement,
        $.call_statement,
        $.filter_statement,
        $.unset_statement,
        $.print_statement,
        $.emit_statement,
        $.tee_statement,
        $.dump_statement,
        // Bare expression as statement (e.g. function call)
        $._expression,
        // Bare block (per-record)
        $.block,
        // Guarded block: (condition) { ... }
        $.guarded_block,
      ),
      // Optional statement terminator
      optional(';'),
    ),

    // --- Blocks ---

    block: $ => prec(1, seq('{', repeat($._statement), '}')),

    begin_block: $ => seq('begin', $.block),
    end_block: $ => seq('end', $.block),

    guarded_block: $ => prec(1, seq(
      '(', $._expression, ')',
      $.block,
    )),

    // --- Functions and subroutines ---

    func_definition: $ => seq(
      'func',
      field('name', $.identifier),
      $.parameter_list,
      optional(seq(':', $.type_name)),
      $.block,
    ),

    subr_definition: $ => seq(
      'subr',
      field('name', $.identifier),
      $.parameter_list,
      $.block,
    ),

    parameter_list: $ => seq(
      '(',
      optional(seq(
        $.parameter,
        repeat(seq(',', $.parameter)),
      )),
      ')',
    ),

    parameter: $ => seq(
      optional($.type_name),
      $.identifier,
    ),

    type_name: _ => choice('var', 'int', 'float', 'num', 'str', 'map', 'bool', 'funct'),

    // --- Control flow ---

    if_statement: $ => seq(
      'if',
      '(', $._expression, ')',
      $.block,
      repeat($.elif_clause),
      optional($.else_clause),
    ),

    elif_clause: $ => seq('elif', '(', $._expression, ')', $.block),
    else_clause: $ => seq('else', $.block),

    for_in_statement: $ => prec(1, seq(
      'for',
      '(',
      choice(
        // Simple: for (k in ...)
        seq($.identifier, 'in', $._expression),
        // Key-value: for (k, v in ...)
        seq($.identifier, ',', $.identifier, 'in', $._expression),
        // Destructured: for ((k1, k2), v in ...)
        seq(
          '(', $.identifier, repeat(seq(',', $.identifier)), ')',
          ',', $.identifier, 'in', $._expression,
        ),
      ),
      ')',
      $.block,
    )),

    for_c_statement: $ => seq(
      'for',
      '(',
      // init
      optional($._for_clause_list),
      ';',
      // condition
      optional($._expression),
      ';',
      // update
      optional($._for_clause_list),
      ')',
      $.block,
    ),

    _for_clause_list: $ => seq(
      $._for_clause_item,
      repeat(seq(',', $._for_clause_item)),
    ),

    _for_clause_item: $ => choice(
      $.assignment,
      $._expression,
    ),

    while_statement: $ => seq(
      'while',
      '(', $._expression, ')',
      $.block,
    ),

    do_while_statement: $ => seq(
      'do',
      $.block,
      'while',
      '(', $._expression, ')',
    ),

    return_statement: $ => prec.right(seq('return', optional($._expression))),
    break_statement: _ => 'break',
    continue_statement: _ => 'continue',

    call_statement: $ => seq('call', $.function_call),

    // --- Output statements ---

    print_statement: $ => seq(
      field('keyword', choice('print', 'printn', 'eprint', 'eprintn')),
      optional(seq($.redirect, ',')),
      $._expression,
      repeat(seq(',', $._expression)),
    ),

    emit_statement: $ => seq(
      field('keyword', choice('emit', 'emit1', 'emitf', 'emitp')),
      optional(seq($.redirect, ',')),
      $._expression,
      repeat(seq(',', $._expression)),
    ),

    tee_statement: $ => seq(
      'tee',
      $.redirect,
      ',',
      $._expression,
    ),

    dump_statement: $ => seq(
      choice('dump', 'edump'),
      optional($.redirect),
    ),

    filter_statement: $ => seq('filter', $._expression),

    unset_statement: $ => seq(
      'unset',
      $._expression,
      repeat(seq(',', $._expression)),
    ),

    redirect: $ => seq(
      choice('>', '>>', '|'),
      $._expression,
      // Redirect target: filename expression or stderr/stdout keyword
    ),

    // --- Assignment ---

    assignment: $ => prec.right(PREC.ASSIGN, seq(
      optional($.type_name),
      field('left', $._lvalue),
      field('operator', choice(
        '=', '+=', '-=', '*=', '/=', '//=', '%=', '**=', '.=',
      )),
      field('right', $._expression),
    )),

    _lvalue: $ => choice(
      $.field_ref,
      $.oosvar_ref,
      $.identifier,
      $.index_expression,
    ),

    // --- Expressions ---

    _expression: $ => choice(
      $.binary_expression,
      $.unary_expression,
      $.ternary_expression,
      $.paren_expression,
      $.function_call,
      $.index_expression,
      $.field_ref,
      $.oosvar_ref,
      $.special_variable,
      $.identifier,
      $.number,
      $.string,
      $.boolean,
      $.array_literal,
      $.map_literal,
    ),

    binary_expression: $ => choice(
      // String concat
      prec.left(PREC.ADD, seq($._expression, '.', $._expression)),
      // Vectorised operators
      prec.left(PREC.ADD, seq($._expression, '.+', $._expression)),
      prec.left(PREC.ADD, seq($._expression, '.-', $._expression)),
      prec.left(PREC.MUL, seq($._expression, '.*', $._expression)),
      prec.left(PREC.MUL, seq($._expression, './', $._expression)),
      // Arithmetic
      prec.left(PREC.ADD, seq($._expression, '+', $._expression)),
      prec.left(PREC.ADD, seq($._expression, '-', $._expression)),
      prec.left(PREC.MUL, seq($._expression, '*', $._expression)),
      prec.left(PREC.MUL, seq($._expression, '/', $._expression)),
      prec.left(PREC.MUL, seq($._expression, '//', $._expression)),
      prec.left(PREC.MUL, seq($._expression, '%', $._expression)),
      prec.right(PREC.POWER, seq($._expression, '**', $._expression)),
      // Comparison
      prec.left(PREC.COMPARE, seq($._expression, '<', $._expression)),
      prec.left(PREC.COMPARE, seq($._expression, '>', $._expression)),
      prec.left(PREC.COMPARE, seq($._expression, '<=', $._expression)),
      prec.left(PREC.COMPARE, seq($._expression, '>=', $._expression)),
      prec.left(PREC.COMPARE, seq($._expression, '<=>', $._expression)),
      // Equality
      prec.left(PREC.EQUALITY, seq($._expression, '==', $._expression)),
      prec.left(PREC.EQUALITY, seq($._expression, '!=', $._expression)),
      // Regex
      prec.left(PREC.EQUALITY, seq($._expression, '=~', $._expression)),
      prec.left(PREC.EQUALITY, seq($._expression, '!=~', $._expression)),
      // Logical
      prec.left(PREC.AND, seq($._expression, '&&', $._expression)),
      prec.left(PREC.OR, seq($._expression, '||', $._expression)),
      prec.left(PREC.XOR, seq($._expression, '^^', $._expression)),
      // Bitwise
      prec.left(PREC.BITAND, seq($._expression, '&', $._expression)),
      prec.left(PREC.BITOR, seq($._expression, '|', $._expression)),
      prec.left(PREC.BITXOR, seq($._expression, '^', $._expression)),
      prec.left(PREC.SHIFT, seq($._expression, '<<', $._expression)),
      prec.left(PREC.SHIFT, seq($._expression, '>>', $._expression)),
      prec.left(PREC.SHIFT, seq($._expression, '>>>', $._expression)),
      // Null coalescing
      prec.left(PREC.NULL_COALESCE, seq($._expression, '??', $._expression)),
      prec.left(PREC.NULL_COALESCE, seq($._expression, '???', $._expression)),
    ),

    unary_expression: $ => choice(
      prec(PREC.UNARY, seq('!', $._expression)),
      prec(PREC.UNARY, seq('~', $._expression)),
      prec(PREC.UNARY, seq('-', $._expression)),
    ),

    ternary_expression: $ => prec.right(PREC.TERNARY, seq(
      $._expression,
      '?',
      $._expression,
      ':',
      $._expression,
    )),

    paren_expression: $ => seq('(', $._expression, ')'),

    function_call: $ => prec(PREC.CALL, seq(
      field('function_name', $.identifier),
      '(',
      optional(seq(
        $._expression,
        repeat(seq(',', $._expression)),
        optional(','),
      )),
      ')',
    )),

    index_expression: $ => prec(PREC.POSTFIX, seq(
      $._expression,
      '[',
      $._expression,
      ']',
    )),

    // --- Field and OOS variable references ---

    field_ref: $ => choice(
      // $name
      seq('$', $.identifier),
      // $*
      seq('$', '*'),
      // $[expr]
      seq('$', '[', $._expression, ']'),
    ),

    oosvar_ref: $ => choice(
      // @name
      seq('@', $.identifier),
      // @*
      seq('@', '*'),
      // @name[key] — handled via index_expression on oosvar_ref
      // @[expr]
      seq('@', '[', $._expression, ']'),
    ),

    special_variable: _ => choice(
      'NR', 'NF', 'FNR', 'FILENAME', 'FILENUM',
      'IFS', 'OFS', 'IRS', 'ORS', 'IPS', 'OPS',
      'M_PI', 'M_E', 'ENV',
      'stdout', 'stderr',
    ),

    // --- Literals ---

    array_literal: $ => seq(
      '[',
      optional(seq(
        $._expression,
        repeat(seq(',', $._expression)),
        optional(','),
      )),
      ']',
    ),

    map_literal: $ => seq(
      '{',
      optional(seq(
        $.map_entry,
        repeat(seq(',', $.map_entry)),
        optional(','),
      )),
      '}',
    ),

    map_entry: $ => seq(
      $._expression,
      ':',
      $._expression,
    ),

    boolean: _ => choice('true', 'false'),

    string: $ => seq(
      '"',
      repeat(choice(
        $.escape_sequence,
        $.string_content,
      )),
      '"',
    ),

    string_content: _ => token.immediate(prec(-1, /[^"\\]+/)),

    escape_sequence: _ => token.immediate(seq(
      '\\',
      choice(
        /[\\'"nrtab0]/,
        /x[0-9a-fA-F]{2}/,
        /u[0-9a-fA-F]{4}/,
      ),
    )),

    number: _ => {
      const decimal = /\d+(\.\d+)?([eE][+-]?\d+)?/;
      const hex = /0[xX][0-9a-fA-F]+/;
      return token(choice(hex, decimal));
    },

    identifier: _ => /[A-Za-z_][A-Za-z0-9_]*/,

    comment: _ => token(seq('#', /.*/)),

    // Lowest-priority catch-all for injection contexts
    _fallback: _ => token(prec(-1, /[^\s]/)),
  },
});
