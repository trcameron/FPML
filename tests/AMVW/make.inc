############################################################
#
# EISCOR (EIgenSolvers based on unitary CORe transformations)
#
############################################################

# directory where the library will be installed
INSTALLDIR := $(CURDIR)

# fortran compiler and flags
FC := gfortran
FFLAGS := -O2

# libs required by compiler
LIBS := 