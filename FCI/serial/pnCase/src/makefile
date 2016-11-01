GCC= gcc -O2

lanczo :	   	lanc-main.o  lanc-pn-SD-basis.o lanc-input.o lanc-pn-sp-basis.o\
			lanc-pn-pot.o lanc-pn-store.o lanc-pn-process.o\
			lanc-pn-matrix-vector.o lanc-lib.o 
			${GCC} lanc-main.o lanc-pn-SD-basis.o lanc-input.o\
			lanc-pn-sp-basis.o lanc-pn-pot.o lanc-pn-store.o\
			lanc-pn-process.o lanc-pn-matrix-vector.o lanc-lib.o -lm -o lanczo

lanc-main.o :			lanc-main.c shell.h 
		        	${GCC} -c lanc-main.c

lanc-input.o :			lanc-input.c shell.h 
		        	${GCC} -c lanc-input.c

lanc-pn-sp-basis.o :		lanc-pn-sp-basis.c shell.h 
		        	${GCC} -c lanc-pn-sp-basis.c

lanc-pn-SD-basis.o :		lanc-pn-SD-basis.c shell.h 
		        	${GCC} -c lanc-pn-SD-basis.c

lanc-pn-pot.o :			lanc-pn-pot.c shell.h 
				${GCC} -c lanc-pn-pot.c

lanc-pn-store.o :		lanc-pn-store.c shell.h 
				${GCC} -c lanc-pn-store.c

lanc-pn-process.o :		lanc-pn-process.c shell.h 
				${GCC} -c lanc-pn-process.c

lanc-pn-matrix-vector.o :	lanc-pn-matrix-vector.c shell.h 
				${GCC} -c lanc-pn-matrix-vector.c

lanc-lib.o :			lanc-lib.c shell.h 
				${GCC} -c lanc-lib.c
