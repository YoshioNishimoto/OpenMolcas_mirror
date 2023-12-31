# -*- coding: utf-8 -*-
#***********************************************************************
# This file is part of OpenMolcas.                                     *
#                                                                      *
# OpenMolcas is free software; you can redistribute it and/or modify   *
# it under the terms of the GNU Lesser General Public License, v. 2.1. *
# OpenMolcas is distributed in the hope that it will be useful, but it *
# is provided "as is" and without any express or implied warranties.   *
# For more details see the full text of the license in the file        *
# LICENSE or in <http://www.gnu.org/licenses/>.                        *
#                                                                      *
# Copyright (C) 2015,2017,2019,2022, Ignacio Fdez. Galván              *
#***********************************************************************

########################################################################
### Patches for pybtex ###

from pybtex.richtext import Text, Protected, String, Tag
from pybtex.scanner import Literal, PybtexSyntaxError
from pybtex.markup import LaTeXParser
from pybtex_docutils import Backend

# Create a new Math element, like Protected but uses a different render
class Math(Protected):
  def __repr__(self):
    reprparts = ', '.join(repr(part) for part in self.parts)
    return 'Math({})'.format(reprparts)
  def render(self, backend):
    text = super(Math, self).render(backend)
    try:
      return backend.format_math(text)
    except AttributeError:
      return backend.format_protected(text)

# Create a new LaTeX parser to output Math elements too
# (only works with $, and does not handle escaping)
# Modifications between the  #vvvvv  and  #^^^^^  lines
class LaTeXMathParser(LaTeXParser):
  #vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  DOLLAR = Literal(u'$')
  #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  def iter_string_parts(self, level=0, in_math=False):
    while True:
      #vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      # there will be no Protected inside Math,
      # since we need to preserve all braces
      if in_math:
        token = self.skip_to([self.DOLLAR])
      else:
        token = self.skip_to([self.LBRACE, self.RBRACE, self.DOLLAR])
      #token = self.skip_to([self.LBRACE, self.RBRACE])
      #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      if not token:
        remainder = self.get_remainder()
        if remainder:
          yield String(remainder)
        if level != 0:
          raise PybtexSyntaxError('unbalanced braces', self)
        break
      #vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      # yield Math element
      elif token.pattern is self.DOLLAR:
        if in_math:
          yield String(token.value[:-1])
          if level == 0:
            raise PybtexSyntaxError('unbalanced math', self)
          break
        else:
          yield String(token.value[:-1])
          yield Math(*self.iter_string_parts(level=level + 1, in_math=True))
      #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      elif token.pattern is self.LBRACE:
        yield String(token.value[:-1])
        yield Protected(*self.iter_string_parts(level=level + 1))
      else:  # brace.pattern is self.RBRACE
        yield String(token.value[:-1])
        if level == 0:
          raise PybtexSyntaxError('unbalanced braces', self)
        break

# Patch for Text.from_latex
@classmethod
def from_latex(cls, latex):
  import codecs
  import latexcodec  # noqa
  return LaTeXMathParser(codecs.decode(latex, 'ulatex')).parse()

Text.from_latex = from_latex

# Patch for Backend
def format_math(self, text):
  import docutils.nodes
  return [docutils.nodes.math('', '', *text)]

# plaintext is used for citation tooltips...
if Backend.default_suffix != '.txt':
  Backend.format_math = format_math

# Function to convert some LaTeX commands
def process_latex(text):
  parts = text.parts
  for i,part in enumerate(parts):
    if isinstance(part, String):
      # Format \emph as Tag('em')
      if (part.endswith(r'\emph ') and (type(parts[i+1]) == Protected)):
        parts[i] = parts[i][:-6]
        parts[i+1] = Tag('em', *parts[i+1].parts)
      # Handle \textsuperscript
      elif (part.endswith(r'\textsuperscript ') and (type(parts[i+1]) == Protected)):
        parts[i] = parts[i][:-17]
        if (str(parts[i+1]) == '2'):
          parts[i+1] = '²'
    elif (type(part) == Protected):
      parts[i] = process_latex(part)
  if (isinstance(text, Text)):
    return Text(*parts)
  elif (type(text) == Protected):
    return Protected(*parts)
  else:
    return text

########################################################################

# Force updating references if anything else changes
def update_bib(app, env, added, changed, removed):
  updated = added.union(changed).union(removed)
  if (app.config.ref_file and (len(updated) > 0)):
    if (app.config.ref_file not in updated):
      return [app.config.ref_file]
  return []

# Put the references file at the end of the list
def bib_at_end(app, env, docnames):
  if (app.config.ref_file):
    try:
      docnames.remove(app.config.ref_file)
      docnames.append(app.config.ref_file)
    except ValueError:
      pass

# Try to detect older sphinxcontrib-bibtex versions
try:
  import sphinxcontrib.bibtex.cache
  manual_sort = True
  # Set global citation order
  citeorder = {}
  # Build "citeorder" according to the citation order,
  # following the document order in sphinx.
  # This has to be done at source-read, because the
  # bibliography is sorted before the doctree-read stage
  def sort_citations(app, docname, source):
    if (docname == app.config.ref_file):
      env = app.builder.env
      cited = env.bibtex_cache._cited
      rel = env.collect_relations()
      # each element is [up, prev, next],
      # so start with the one with prev=None, and follow next
      doc = [x for x in rel if rel[x][1] is None][0]
      i = 0
      while (doc is not None):
        for cite in cited[doc]:
          if (cite not in citeorder):
            citeorder[cite] = i
            i = i+1
        doc = rel[doc][2]
except ImportError:
  manual_sort = False

def setup(app):
  app.connect('env-get-outdated', update_bib)
  app.connect('env-before-read-docs', bib_at_end)
  if manual_sort:
    app.connect('source-read', sort_citations)
  app.add_config_value('ref_file', 'references', 'html')


### Bibliography Styles ###

from pybtex.style.names import BaseNameStyle, name_part
from pybtex.style.formatting.unsrt import dashify, date, Style as UnsrtStyle
from pybtex.style.formatting import toplevel
from pybtex.style.template import join, words, field, optional, first_of, names, sentence, tag, optional_field, href
from pybtex.plugin import register_plugin

if (manual_sort):
  from pybtex.style.sorting import BaseSortingStyle
  # Trivial sorting style that just uses "citeorder" as key
  class CiteStyle(BaseSortingStyle):
    def sorting_key(self, entry):
      try:
        return citeorder[entry.key]
      except:
        return 9999999
  register_plugin('pybtex.style.sorting', 'cite', CiteStyle)

# Formatting style
class MolcasStyle(UnsrtStyle):

  if (manual_sort):
    default_sorting_style = 'cite'

  date = words [optional_field('month'), field('year')]

  def format_title(self, e, which_field, as_sentence=True):
    formatted_title = field(
      which_field, apply_func=lambda text: process_latex(text).capitalize()
    )
    if as_sentence:
      return sentence [ formatted_title ]
    else:
      return formatted_title

  def format_btitle(self, e, which_field, as_sentence=True):
    formatted_title = tag('em') [ field(
      which_field, apply_func=lambda text: process_latex(text)
    ) ]
    if as_sentence:
      return sentence [ formatted_title ]
    else:
      return formatted_title

  def format_names(self, role, as_sentence=True):
    formatted_names = names(role, sep=', ', sep2 = ', ', last_sep=', ')
    if as_sentence:
      return sentence [ formatted_names ]
    else:
      return formatted_names

  def get_article_template(self, e):
    pages = first_of [
      # article id with total pages
      optional [
        join [
          field('articleid'),
          optional [
            '(',join[u'1–', optional_field('pagetotal')],')'
          ]
        ]
      ],
      # pages only
      field('pages', apply_func=dashify),
    ]
    volume_and_pages = first_of [
      # volume and pages, with optional issue number
      optional [
        join [
          tag('strong') [field('volume')],
          optional [ '[', field('number'),']' ],
          ' (', field('year'), ')',
          ' ', pages
        ],
      ],
      # pages only
      words ['pages', pages],
    ]
    template = toplevel [
      self.format_names('author'),
      self.format_title(e, 'title'),
      sentence [
        tag('em') [field('journal')],
        optional [ volume_and_pages ],
      ],
      sentence [ optional_field('note') ],
      self.format_web_refs(e),
    ]
    return template

  def get_inbook_template(self, e):
    template = toplevel [
      self.format_author_or_editor(e),
      sentence [
        self.format_btitle(e, 'title', as_sentence=False),
        self.format_chapter_and_pages(e),
      ],
      self.format_volume_and_series(e),
      sentence [
        field('publisher'),
        optional_field('address'),
        optional [
          words [field('edition'), 'edition']
        ],
        date,
        optional_field('note'),
      ],
      self.format_web_refs(e),
    ]
    return template

  def get_incollection_template(self, e):
    template = toplevel [
      sentence [ self.format_names('author') ],
      self.format_title(e, 'title'),
      words [
        'In',
        sentence [
          optional [ self.format_editor(e, as_sentence=False) ],
          self.format_btitle(e, 'booktitle', as_sentence=False),
          self.format_volume_and_series(e, as_sentence=False),
          self.format_chapter_and_pages(e),
        ],
      ],
      sentence [
        optional_field('publisher'),
        optional_field('address'),
        self.format_edition(e),
        date,
      ],
      self.format_web_refs(e),
    ]
    return template

  def get_book_template(self, e):
    template = toplevel [
      self.format_author_or_editor(e),
      self.format_btitle(e, 'title'),
      self.format_volume_and_series(e),
      sentence [
        field('publisher'),
        optional_field('address'),
        self.format_edition(e),
        date
      ],
      optional [ sentence [ self.format_isbn(e) ] ],
      sentence [ optional_field('note') ],
      self.format_web_refs(e),
    ]
    return template

  def get_misc_template(self, e):
    template = toplevel [
      optional [ sentence [ self.format_names('author') ] ],
      optional [ self.format_title(e, 'title') ],
      sentence [
        optional [ field('howpublished') ],
        optional [ date ],
      ],
      sentence [ optional_field('note') ],
      self.format_web_refs(e),
    ]
    return template

  def get_phdthesis_template(self, e):
    template = toplevel [
      sentence [ self.format_names('author') ],
      self.format_btitle(e, 'title'),
      sentence [
        'PhD thesis',
        field('school'),
        optional_field('address'),
        date,
      ],
      sentence [ optional_field('note') ],
      self.format_web_refs(e),
    ]
    return template

  def get_manual_template(self, e):
    template = toplevel [
      optional [ sentence [ self.format_names('author') ] ],
      self.format_btitle(e, 'title'),
      sentence [
        optional_field('organization'),
        optional_field('address'),
        self.format_edition(e),
        optional [ date ],
      ],
      sentence [ optional_field('note') ],
      self.format_web_refs(e),
    ]
    return template

  def format_doi(self, e):
    return href [
      join [
        'https://doi.org/',
        field('doi')
      ],
      join [
        'doi:',
        field('doi')
      ]
    ]

# Shorter style for tooltips
class ShortStyle(MolcasStyle):

  def __init__(self, **kwargs):
    super(MolcasStyle, self).__init__(abbreviate_names=True, **kwargs)

  def format_web_refs(self, e):
    return []

  def format_names(self, role, as_sentence=True):
    formatted_names = names(role, sep=', ', sep2 = ', ', last_sep=', ')
    if as_sentence:
      return sentence [ formatted_names ]
    else:
      return formatted_names

  def get_article_template(self, e):
    pages = first_of [
      # article id with total pages
      optional [
        join [
          field('articleid'),
          optional [
            '(',join[u'1–', optional_field('pagetotal')],')'
          ]
        ]
      ],
      # pages only
      field('pages', apply_func=dashify),
    ]
    volume_and_pages = first_of [
      # volume and pages, with optional issue number
      optional [
        join [
          tag('strong') [field('volume')],
          optional [ '[', field('number'),']' ],
          ' (', field('year'), ')',
          ' ', pages
        ],
      ],
      # pages only
      words ['pages', pages],
    ]
    template = toplevel [
      self.format_names('author'),
      sentence [
        tag('em') [field('journal')],
        optional [ volume_and_pages ],
      ],
      sentence [ optional_field('note') ],
      self.format_web_refs(e),
    ]
    return template

  def get_inbook_template(self, e):
    template = toplevel [
      self.format_author_or_editor(e),
      sentence [
        self.format_btitle(e, 'title', as_sentence=False),
        self.format_chapter_and_pages(e),
      ],
      self.format_volume_and_series(e),
      sentence [
        field('publisher'),
        date,
      ],
      self.format_web_refs(e),
    ]
    return template

  def get_incollection_template(self, e):
    template = toplevel [
      sentence [ self.format_names('author') ],
      words [
        'In',
        sentence [
          self.format_btitle(e, 'booktitle', as_sentence=False),
          self.format_volume_and_series(e, as_sentence=False),
          self.format_chapter_and_pages(e),
        ],
      ],
      sentence [
        optional_field('publisher'),
        date,
      ],
      self.format_web_refs(e),
    ]
    return template

  def get_book_template(self, e):
    template = toplevel [
      self.format_author_or_editor(e),
      self.format_btitle(e, 'title'),
      self.format_volume_and_series(e),
      sentence [
        field('publisher'),
        date
      ],
      sentence [ optional_field('note') ],
      self.format_web_refs(e),
    ]
    return template

  def get_phdthesis_template(self, e):
    template = toplevel [
      sentence [ self.format_names('author') ],
      sentence [
        'PhD thesis',
        field('school'),
        optional_field('address'),
        date,
      ],
      sentence [ optional_field('note') ],
      self.format_web_refs(e),
    ]
    return template

register_plugin('pybtex.style.formatting', 'molcas', MolcasStyle)
register_plugin('pybtex.style.formatting', 'short', ShortStyle)
