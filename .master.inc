#############################################
#           FPML Master Make Include File   #
#############################################
# These variables should not be changed as
# they are specific to the installation of
# AMVW eiscor library.
#############################################

# directory where the eiscor library will be installed
# $(shell pwd) returns the current directory of the AMVW makefile
INSTALLDIR := $(shell pwd)

# current version of AMVW eiscor
LIBNAME := eiscor
MAJOR := 0
MINOR := 2
PATCH := 0
VERSION := $(MAJOR).$(MINOR).$(PATCH)
UNAME := $(shell uname)

# change library extension based on OS
ifeq ($(findstring Windows_NT,$(OS)),Windows_NT)
	SLIB := dll
else ifeq ($(UNAME),Darwin)
	SLIB := dylib
else
	SLIB := so
endif