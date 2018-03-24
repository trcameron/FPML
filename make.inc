#############################################
#           FPML Make Include File          #
#############################################

# fortran compiler and flags
FC := gfortran
FFLAGS := -O2

# tests run variables
SD := 80		# start degree

ED := 2560		# end degree

IT := 10		# iteration number

# latex compiler and flags
LC := pdflatex
LFLAGS := -shell-escape