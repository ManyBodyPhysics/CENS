GCC= gcc  -O3 -Wall


trans-pn :          transition-main.o  id-transition.o pn-transition.o  transition-lib.o   
		    ${GCC} transition-main.o id-transition.o pn-transition.o\
                                                transition-lib.o -lm -o trans-pn

transition-main.o : 	transition-main.c shell.h sm.h
			${GCC} -c transition-main.c

id-transition.o :   	id-transition.c shell.h sm.h
		    	${GCC} -c id-transition.c

pn-transition.o :	pn-transition.c  shell.h sm.h
		    	${GCC} -c pn-transition.c

transition-lib.o: 	transition-lib.c shell.h sm.h
		    	${GCC} -c transition-lib.c
