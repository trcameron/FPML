#!/bin/bash

# compile latex
pdflatex -shell-escape *.tex
cp *.png ../figures/

# clean up
rm *.log
rm *.aux
rm *.pdf
rm *.png