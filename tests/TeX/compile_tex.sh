#!/bin/bash

# compile latex
pdflatex -shell-escape rand_poly.tex
pdflatex -shell-escape unity.tex
pdflatex -shell-escape spec_poly_results.tex
pdflatex -shell-escape spec_poly_list.tex
pdflatex -shell-escape methods.tex
pdflatex -shell-escape module.tex

# copy png files into figures
cp *.png ../figures/

# clean up
rm *.log
rm *.aux
rm *.pdf
rm *.png