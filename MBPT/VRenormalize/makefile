#XLF= ifort -O3 -static  -heap-arrays 1024 -g -save-temps -fpic -Wl,-no_pie
XLF= ifort -O3 -static  -heap-arrays 1024

vrenorm.exe :   renorm-modules.o chiral-twobody-potentials.o chp_presets.o renorm-main.o renorm-potentials.o chinn6.o renorm-vkrg.o renorm-gmatrix.o  renorm-nocore.o renorm-library.o renorm-vlowk.o
	       ${XLF} -o  vrenorm.exe  renorm-modules.o chiral-twobody-potentials.o chp_presets.o renorm-main.o renorm-potentials.o chinn6.o renorm-vkrg.o renorm-gmatrix.o renorm-nocore.o renorm-library.o renorm-vlowk.o

chinn6.o :   chinn6.f
		${XLF} -c    chinn6.f

chiral-twobody-potentials.o :	chiral-twobody-potentials.f90
		${XLF} -c    chiral-twobody-potentials.f90

chp_presets.o :   chp_presets.f90
		${XLF} -c    chp_presets.f90

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
