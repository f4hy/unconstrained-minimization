#
# ifort -C -CB -132 implicitpde.f

compiler	= gfortran
mpicompiler	= mpif90
optflags	= -O3 -march=native
flags		= -llapack -ffixed-line-length-132 -fbounds-check -W -Wall -Wextra -fdefault-real-8

intelcc		= ifort
intelflags	= -CB -r8 -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -lpthread

openmp	= -fopenmp


rosen: rosen.exe rosen.intel.exe

tests:	 runtests.exe


runtests.exe: testmethods.f90 methods.f90
	${compiler} ${flags} -o runtests.exe testmethods.f90 methods.f90

rosen.exe: rosen.f90  main.f90 methods.f90  linesearch.f90
	${compiler} ${flags} -o rosen.exe rosen.f90 methods.f90 main.f90  linesearch.f90

rosen.intel.exe: rosen.f90  main.f90 methods.f90  linesearch.f90
	${intelcc} ${intelflags} -o rosen.intel.exe rosen.f90 methods.f90 main.f90  linesearch.f90

clean:
	-rm *.exe
