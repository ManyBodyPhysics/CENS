#XLF= ifort -O3 -static  -heap-arrays 1024
XLF= gfortran -O3 

vrenorm.exe :   renorm-modules.o renorm-main.o renorm-potentials.o renorm-vkrg.o renorm-gmatrix.o  renorm-nocore.o renorm-library.o renorm-vlowk.o 
	       ${XLF} -o  vrenorm.exe  renorm-modules.o renorm-main.o renorm-potentials.o renorm-vkrg.o renorm-gmatrix.o renorm-nocore.o renorm-library.o renorm-vlowk.o 

renorm-modules.o :   renorm-modules.f90
		${XLF} -c  renorm-modules.f90

renorm-main.o :   renorm-main.f90
		${XLF} -c renorm-main.f90

renorm-gmatrix.o :   renorm-gmatrix.f90
		${XLF} -c   renorm-gmatrix.f90

renorm-nocore.o :   renorm-nocore.f90
		${XLF} -c   renorm-nocore.f90

renorm-vlowk.o :   renorm-vlowk.f90
		${XLF} -c   renorm-vlowk.f90

renorm-vkrg.o :   renorm-vkrg.f90
		${XLF} -c   renorm-vkrg.f90

renorm-potentials.o :   renorm-potentials.f
		${XLF} -c    renorm-potentials.f



renorm-library.o :   renorm-library.f
		${XLF} -c    renorm-library.f

