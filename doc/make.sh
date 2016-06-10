#!/bin/bash

### listings
echo "import pygments.formatters as f; print f.LatexFormatter().get_style_defs()" | python > pygments.tex

pygmentize -l python -o lst/plot.py.tex  examples/plot.py
pygmentize -l idl    -o lst/plot.pro.tex examples/plot.pro

### figures
#rsvg-convert -f pdf -o img/example_1.pdf     ../libmpdataxx/build/tests/paper_example_1/out.svg
#pdfcrop img/example_1.pdf
#rm img/example_1.pdf

#for opt in abs iga iga_tot iga_fct iga_tot_fct; do
#  rsvg-convert -f pdf -o img/example_2_${opt}.pdf     ../libmpdataxx/build/tests/paper_example_2/out_${opt}.svg
#  pdfcrop img/example_2_${opt}.pdf
#  rm img/example_2_${opt}.pdf
#done;

### LaTeX
pdflatex doc.tex
bibtex doc.aux
pdflatex doc.tex
pdflatex doc.tex

