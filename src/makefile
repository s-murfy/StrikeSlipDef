# type of compilation
TYPE = debug
#TYPE = run
# compiler type
FC = gfortran
#FC = ifort
# flags for debugging
ifeq ($(TYPE),debug)
 ifeq ($(FC),ifort)
   FCFLAGS = -debug -g -openmp -parallel -save-temps -traceback -fpe:0 -Wl, -fpic -warn all -check all -assume byterecl
 else
   FCFLAGS = -g -fopenmp -fbounds-check -fbacktrace
 endif
endif
ifeq ($(TYPE),run)
 ifeq ($(FC),ifort)
# FCFLAGS = -O2 -openmp -parallel -axSSE4.2 -assume byterecl
    FCFLAGS = -O2 -openmp -parallel -axAVX -assume byterecl
  else
    FCFLAGS = -O2 -fopenmp
 endif
endif
 FCFLAGS += -I/usr/local/include
# Libraries to be included
 LIB = #-lfftw3
# Executable name
EXEC = main
all:$(EXEC)

main.o: georef.o initiation.o okada_sub.o stiffness.o
main: georef.o initiation.o okada_sub.o stiffness.o
#	main.o:constants.o georef.o utils.o initiation.o okada_sub.o stiffness.o
#	main:constants.o georef.o utils.o initiation.o okada_sub.o stiffness.o

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LIB)
%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

# clean up
clean:
	rm -f *.o *.mod *.s
# clean up
vclean:
	rm -f *.o *.mod err  *.csv  fort.*
