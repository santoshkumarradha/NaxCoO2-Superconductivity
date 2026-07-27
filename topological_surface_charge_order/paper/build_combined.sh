#!/bin/bash
# Build the single review PDF: the article followed by its Supplemental Material.
# Both documents carry the same pastel hyperref palette; the supplement prefixes
# every PDF destination name (\HyperDestNameFilter) so that internal links
# survive concatenation without colliding with the article's.
set -e
cd "$(dirname "$0")"
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex >/dev/null
latexmk -pdf -interaction=nonstopmode -halt-on-error supplement.tex >/dev/null
if command -v qpdf >/dev/null 2>&1; then
  qpdf --empty --pages main.pdf supplement.pdf -- main_with_sm.pdf
  echo "wrote main_with_sm.pdf ($(du -h main_with_sm.pdf | cut -f1))"
else
  echo "qpdf not found; main.pdf and supplement.pdf built separately"
fi
