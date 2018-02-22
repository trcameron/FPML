#!/bin/bash

# run test_modules
./test_modules2 1000 1000000 100

# compile test_modules.tex
pdflatex -shell-escape test_start_timing.tex
pdflatex -shell-escape test_start_accuracy.tex
pdflatex -shell-escape test_correction_timing.tex
cp *.png figures/

# clean up
rm *.log
rm *.aux
rm *.pdf
rm *.png