#############################################
#           FPML Make File TeX Directory    #
#############################################
include ../../make.inc

all: conv nat_poly rand_poly rand_unity rev_nat_poly spec_poly_list spec_poly_results trunc_exp unity clean
	
conv:
	$(LC) conv.tex
	
nat_poly:
	$(LC) nat_poly.tex
	
rand_poly:
	$(LC) rand_poly.tex
	
rand_unity:
	$(LC) rand_unity.tex

rev_nat_poly:
	$(LC) rev_nat_poly.tex
	
spec_poly_list:
	$(LC) spec_poly_list.tex
	
spec_poly_results:
	$(LC) spec_poly_results.tex
	
trunc_exp:
	$(LC) trunc_exp.tex
	
unity:
	$(LC) unity.tex
	
clean: copy
	@rm -r *.log
	@rm -r *.aux
	@rm -r *.pdf
	
copy:
	@cp -r *.pdf ../figures/
