#
# ifort -C -CB -132 implicitpde.f

compiler	= gfortran
mpicompiler	= mpif90
optflags	= -O3 -march=native
flags		= -llapack -ffixed-line-length-132 -fbounds-check -W -Wall -Wextra
openmp	= -fopenmp


tests:	 runtests.exe


runtests.exe: testmethods.f90 methods.f90
	${compiler} ${flags} -o runtests.exe testmethods.f90 methods.f90


clean:
	-rm *.exe
