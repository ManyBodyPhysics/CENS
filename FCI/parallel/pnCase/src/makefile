
#    GCC=mpicc  -fpack-struct
#    GCC=mpicc -ccl icc -fpack-struct

  GCC=mpicc -fpack-struct


PAR-lanczo :		PAR-pnShellModelMain.o PAR-pnShellModelInput.o\
			PAR-pnSP-Basis.o PAR-pnSD-Basis.o\
			PAR-idTwoPartIntJ.o PAR-idTwoPartAngJ.o\
			PAR-pnTwoPartIntJ.o PAR-pnTwoPartAngJ.o\
			PAR-pnStoreSD_MatrElem.o PAR-idNondiagCalc.o\
			PAR-pnLancProcess.o PAR-pnMatrVectorCalc.o\
			PAR-pnLancResult.o\
			PAR-shellModel-Lib.o

			$(GCC) PAR-pnShellModelMain.o PAR-pnShellModelInput.o\
			PAR-pnSP-Basis.o PAR-pnSD-Basis.o\
			PAR-idTwoPartIntJ.o PAR-idTwoPartAngJ.o\
			PAR-pnTwoPartIntJ.o PAR-pnTwoPartAngJ.o\
			PAR-pnStoreSD_MatrElem.o PAR-idNondiagCalc.o\
			PAR-pnLancProcess.o PAR-pnMatrVectorCalc.o\
			PAR-pnLancResult.o\
			PAR-shellModel-Lib.o -lm -o PAR-lanczo


PAR-pnShellModelMain.o :	PAR-pnShellModelMain.c PAR-pnShellModel.h
				${GCC} -c PAR-pnShellModelMain.c

PAR-pnShellModelInput.o :	PAR-pnShellModelInput.c PAR-pnShellModel.h
				${GCC} -c PAR-pnShellModelInput.c

PAR-idTwoPartIntJ.o :		PAR-idTwoPartIntJ.c  PAR-pnShellModel.h
				${GCC} -c  PAR-idTwoPartIntJ.c

PAR-idTwoPartAngJ.o :		PAR-idTwoPartAngJ.c  PAR-pnShellModel.h
				${GCC} -c  PAR-idTwoPartAngJ.c

PAR-pnTwoPartIntJ.o :		PAR-pnTwoPartIntJ.c  PAR-pnShellModel.h
				${GCC} -c  PAR-pnTwoPartIntJ.c

PAR-pnTwoPartAngJ.o :		PAR-pnTwoPartAngJ.c  PAR-pnShellModel.h
				${GCC} -c  PAR-pnTwoPartAngJ.c

PAR-pnSP-Basis.o :		PAR-pnSP-Basis.c PAR-pnShellModel.h
				${GCC} -c  PAR-pnSP-Basis.c

PAR-pnSD-Basis.o :		PAR-pnSD-Basis.c PAR-pnShellModel.h
				${GCC} -c  PAR-pnSD-Basis.c

PAR-pnStoreSD_MatrElem.o :	PAR-pnStoreSD_MatrElem.c PAR-pnShellModel.h
				${GCC} -c  PAR-pnStoreSD_MatrElem.c

PAR-idNondiagCalc.o :		PAR-idNondiagCalc.c PAR-pnShellModel.h
				${GCC} -c  PAR-idNondiagCalc.c


PAR-pnLancProcess.o :		PAR-pnLancProcess.c PAR-pnShellModel.h
				${GCC} -c PAR-pnLancProcess.c

PAR-pnMatrVectorCalc.o :	PAR-pnMatrVectorCalc.c PAR-pnShellModel.h
				${GCC} -c PAR-pnMatrVectorCalc.c

PAR-pnLancResult.o :		PAR-pnLancResult.c PAR-pnShellModel.h
				${GCC} -c  PAR-pnLancResult.c

PAR-shellModel-Lib.o :		PAR-shellModel-Lib.c PAR-pnShellModel.h
				${GCC} -c  PAR-shellModel-Lib.c
