GCC= cc -O3 -Wall

Lanczo :		shell-model-main.o  shell-model-input.o\
			id-sp-basis.o id-SD-basis.o\
			id-two-part-angJ.o id-two-part-intJ.o\
			id-two-part-veffM.o id-three-part-veffM.o\
			idStoreSD_MatrElem.o id-nondiag-calc.o\
			lanc-id-process.o lanc-id-matrix-vector-calc.o\
			id-lanczo-matrix-element.o\
			id-lanc-result.o shell-model-lib.o

			$(GCC) shell-model-main.o shell-model-input.o\
			id-sp-basis.o id-SD-basis.o\
			id-two-part-angJ.o id-two-part-intJ.o\
			id-two-part-veffM.o id-three-part-veffM.o\
			idStoreSD_MatrElem.o id-nondiag-calc.o\
			lanc-id-process.o lanc-id-matrix-vector-calc.o\
			id-lanczo-matrix-element.o\
			id-lanc-result.o shell-model-lib.o -lm -o lanczo


shell-model-main.o :	shell-model-main.c shell-model.h 
			${GCC} -c  shell-model-main.c

shell-model-input.o :	shell-model-input.c shell-model.h 
			${GCC} -c shell-model-input.c

id-sp-basis.o :		id-sp-basis.c shell-model.h
			${GCC} -c id-sp-basis.c

id-SD-basis.o:			id-SD-basis.c shell-model.h 
				$(GCC) -c id-SD-basis.c

id-two-part-angJ.o:		id-two-part-angJ.c shell-model.h 
				$(GCC) -c id-two-part-angJ.c

id-two-part-intJ.o:		id-two-part-intJ.c shell-model.h 
				$(GCC) -c id-two-part-intJ.c

id-two-part-veffM.o:		id-two-part-veffM.c shell-model.h 
				$(GCC) -c id-two-part-veffM.c

id-three-part-veffM.o:		id-three-part-veffM.c shell-model.h 
				$(GCC) -c id-three-part-veffM.c

idStoreSD_MatrElem.o:		idStoreSD_MatrElem.c shell-model.h 
				$(GCC) -c idStoreSD_MatrElem.c

id-nondiag-calc.o:		id-nondiag-calc.c shell-model.h 
				$(GCC) -c id-nondiag-calc.c

lanc-id-process.o:		lanc-id-process.c shell-model.h 
				$(GCC) -c lanc-id-process.c
lanc-id-matrix-vector-calc.o:	lanc-id-matrix-vector-calc.c shell-model.h 
				$(GCC) -c lanc-id-matrix-vector-calc.c

id-lanczo-matrix-element.o:	id-lanczo-matrix-element.c  shell-model.h 
				$(GCC) -c id-lanczo-matrix-element.c

id-lanc-result.o:		id-lanc-result.c shell-model.h 
				$(GCC) -c id-lanc-result.c

shell-model-lib.o :		shell-model-lib.c shell-model.h
				${GCC} -c shell-model-lib.c

