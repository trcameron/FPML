##############################################################################
#  example/Makefile
#
#  This software is provided for research purposes only.
#  Commercial usage requires license agreement.
#
#  This work was supported by the Director, Office of Science, Division
#  of Mathematical, Information, and Computational Sciences of the
#  U.S. Department of Energy under contract number DE-AC03-76SF00098.
# 
#  Copyright (c) 2001-2002-2003-2004
#  
#  Modified:       2008-08-25
#  Version:        1.0
#
##############################################################################

TARGET = testmp90 tmpmod90 pslq1 pslq2 pslq3 pslqm1 pslqm2 pslqm3 quaderf quadgs quadts roots

all: $(TARGET)

###### Fortran Test Programs ######
F90 = gfortran

F90FLAGS = -O2 -ffree-form

LD = gfortran

F90MOD	= mpfun90.o mpmod90.o mpmodm90.o mpmodx90.o

testmp90:  $(F90MOD) testmp90.o second.o
	$(LD) -o $@ testmp90.o second.o $(F90MOD)

tmpmod90:  $(F90MOD) tmpmod90.o
	$(LD) -o $@ tmpmod90.o $(F90MOD)

pslq1:  $(F90MOD) tpslq1.o second.o
	$(LD) -o $@ tpslq1.o second.o $(F90MOD)

pslq2:  $(F90MOD) tpslq2.o second.o
	$(LD) -o $@ tpslq2.o second.o $(F90MOD)

pslq3:  $(F90MOD) tpslq3.o second.o
	$(LD) -o $@ tpslq3.o second.o $(F90MOD)

pslqm1:  $(F90MOD) tpslqm1.o second.o
	$(LD) -o $@ tpslqm1.o second.o $(F90MOD)

pslqm2:  $(F90MOD) tpslqm2.o second.o
	$(LD) -o $@ tpslqm2.o second.o $(F90MOD)

pslqm3:  $(F90MOD) tpslqm3.o second.o
	$(LD) -o $@ tpslqm3.o second.o $(F90MOD)

quadgs:  $(F90MOD) tquadgs.o second.o
	$(LD) -o $@ tquadgs.o second.o $(F90MOD)

quaderf:  $(F90MOD) tquaderf.o second.o
	$(LD) -o $@ tquaderf.o second.o $(F90MOD)

quadts:  $(F90MOD) tquadts.o second.o
	$(LD) -o $@ tquadts.o second.o $(F90MOD)

roots:  $(F90MOD) roots.o second.o
	$(LD) -o $@ roots.o second.o $(F90MOD)

###### General Rules ######
.f.o:
	$(F90) $(F90FLAGS) -c $<

###### Clean Up ######
clean: 
	@rm -f *.o *.mod *~ *.d work.pc work.pcl $(TARGET)


