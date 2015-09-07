


#    GCC=mpicc  -fpack-struct
#    GCC=mpicc -ccl icc -fpack-struct

GCC=mpicc   -fpack-struct


##GCC=mpicc
#GCC = icc -g -I/site/VERSIONS/openmpi-1.2.8.intel/include  -L/site/VERSIONS/openmpi-1.2.8.intel/lib -lmpi
#GCC  = gcc -g -I/opt/scali/include64 -L/opt/scali/lib64 -lmpi 

# GCC= scampicc -O3 

Lanczo :		PAR-shell-model-main.o  PAR-shell-model-input.o\
			PAR-id-sp-basis.o PAR-id-SD-basis.o\
			PAR-id-two-part-angJ.o PAR-id-two-part-intJ.o\
			PAR-id-two-part-veffM.o PAR-id-three-part-veffM.o\
			PAR-idStoreSD-MatrElem.o PAR-id-nondiag-calc.o\
			PAR-lanc-id-process.o PAR-id-MatrVectorCalc.o\
		 	PAR-id-lanc-result.o PAR-shell-model-lib.o

			$(GCC) PAR-shell-model-main.o PAR-shell-model-input.o\
			PAR-id-sp-basis.o PAR-id-SD-basis.o\
			PAR-id-two-part-angJ.o PAR-id-two-part-intJ.o\
			PAR-id-two-part-veffM.o PAR-id-three-part-veffM.o\
			PAR-idStoreSD-MatrElem.o PAR-id-nondiag-calc.o\
			PAR-lanc-id-process.o PAR-id-MatrVectorCalc.o\
			PAR-id-lanc-result.o PAR-shell-model-lib.o -lm -o PAR-lanczo 

PAR-shell-model-main.o :	PAR-shell-model-main.c PAR-shell-model.h 
				${GCC} -c  PAR-shell-model-main.c

PAR-shell-model-input.o :	PAR-shell-model-input.c PAR-shell-model.h 
				${GCC} -c PAR-shell-model-input.c

PAR-id-sp-basis.o :		PAR-id-sp-basis.c PAR-shell-model.h
				${GCC} -c PAR-id-sp-basis.c

PAR-id-SD-basis.o:		PAR-id-SD-basis.c PAR-shell-model.h 
				$(GCC) -c PAR-id-SD-basis.c

PAR-id-two-part-angJ.o:		PAR-id-two-part-angJ.c PAR-shell-model.h 
				$(GCC) -c PAR-id-two-part-angJ.c

PAR-id-two-part-intJ.o:		PAR-id-two-part-intJ.c PAR-shell-model.h 
				$(GCC) -c PAR-id-two-part-intJ.c

PAR-id-two-part-veffM.o:	PAR-id-two-part-veffM.c PAR-shell-model.h 
				$(GCC) -c PAR-id-two-part-veffM.c

PAR-id-three-part-veffM.o:	PAR-id-three-part-veffM.c PAR-shell-model.h 
				$(GCC) -c PAR-id-three-part-veffM.c

PAR-idStoreSD-MatrElem.o:	PAR-idStoreSD-MatrElem.c PAR-shell-model.h 
				$(GCC) -c PAR-idStoreSD-MatrElem.c

PAR-id-nondiag-calc.o:		PAR-id-nondiag-calc.c PAR-shell-model.h 
				$(GCC) -c PAR-id-nondiag-calc.c

PAR-lanc-id-process.o:		PAR-lanc-id-process.c PAR-shell-model.h 
				$(GCC) -c PAR-lanc-id-process.c

PAR-id-MatrVectorCalc.o:	PAR-id-MatrVectorCalc.c PAR-shell-model.h 
				$(GCC) -c PAR-id-MatrVectorCalc.c

PAR-id-lanc-result.o:		PAR-id-lanc-result.c PAR-shell-model.h 
				$(GCC) -c PAR-id-lanc-result.c

PAR-shell-model-lib.o :		PAR-shell-model-lib.c PAR-shell-model.h
				${GCC} -c PAR-shell-model-lib.c
