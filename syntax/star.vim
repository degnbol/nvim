" Vim syntax file
" Language:     STAR (mmCIF) file
" Maintainer:   Laszlo Kajan <lkajan@rostlab.org>
" License:      Distributed under the GNU General Public License (version 3 or above)
" Last Change:  2012 July 19
" URL:		ftp://rostlab.org/free/star.vim

" Standard syntax initialization
if version < 600
  syntax clear
elseif exists("b:current_syntax")
  finish
endif

" Should match case
syn case match

syn keyword	starRepeat	loop_
syn match	starKeyword	"\<\(save\|stop\|data\)_"
syn match	starComment	"#.*" contains=starTodo,@Spell
syn match	starIdentifier	"\<_\S\+"
syn match	starIdentifier	"\<_\S\+" nextgroup=starDataItem skipwhite
syn cluster	starTextString	contains=starRepeat,starKeyword,starIdentifier
syn match	starString	"\(^\|\s\@<=\)[^#[:space:]]\S*" contains=@starTextString
syn region	starString	matchgroup=starStringStartEnd start="'" skip="#" end="'" keepend contains=@starTextString
syn region	starString	matchgroup=starStringStartEnd start=+"+ skip="#" end=+"+ keepend contains=@starTextString
syn region	starString	matchgroup=starStringStartEnd start="^;" skip="#" end="^;" keepend contains=@starTextString

" synchronizing

" The default highlighting.
hi def link starRepeat		Repeat
hi def link starKeyword		Keyword
hi def link starComment		Comment
hi def link starString		String
hi def link starStringStartEnd	starString
hi def link starIdentifier	Identifier

let b:current_syntax = "star"

" vim: ts=8 sw=2
