#
# ifort -C -CB -132 implicitpde.f

# compiler	= gfortran
# mpicompiler	= mpif90
# optflags	= -O3 -march=native
# flags		= -llapack -ffixed-line-length-132 -fbounds-check -W -Wall -Wextra -fdefault-real-8

compiler		= ifort
flags	= -warn -warn nointerfaces  -O3 -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -lpthread
double = -r8
quad = -r16
openmp	= -fopenmp

main = methods.f90 main.f90 linesearch.f90 dogleg.f90 finite.f90

all: rosen powell wood

rosen: rosen.exe rosen.16.exe

wood: wood.exe wood.16.exe

powell: powell.exe powell.16.exe

rosen.exe: rosen.f90  ${main}
	${compiler} ${flags} ${double} -o rosen.exe rosen.f90 ${main}

powell.exe: powell.f90  ${main}
	${compiler} ${flags} ${double} -o powell.exe powell.f90 ${main}

wood.exe: wood.f90  ${main}
	${compiler} ${flags} ${double} -o wood.exe wood.f90 ${main}

rosen.16.exe: rosen.f90  ${main}
	${compiler} ${flags} ${quad} -o rosen.16.exe rosen.f90 ${main}

powell.16.exe: powell.f90  ${main}
	${compiler} ${flags} ${quad} -o powell.16.exe powell.f90 ${main}

wood.16.exe: wood.f90  ${main}
	${compiler} ${flags} ${quad} -o wood.16.exe wood.f90 ${main}


clean:
	-rm *.exe
