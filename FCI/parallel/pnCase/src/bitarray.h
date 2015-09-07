     /*
     ** The program module                                   
     **        butarray.h
     ** for testBitArray.c calculation.
     */

typedef unsigned char BARR_ELTYPE;

#define BARR_ELBITS (8*(sizeof(BARR_ELTYPE)))

#define BARR_ARRAYSIZE(N) (( (N) + BARR_ELBITS-1) / BARR_ELBITS)

#define BARR_MALLOC(N)  \
         ((BARR_ELTYPE *)malloc(BARR_ARRAYSIZE(N)*sizeof(BARR_ELTYPE)))

#define BARR_REALLOC(barr,N) \
           ((BARR_ELTYPE *)remalloc(barr, BARR_ARRAYSIZE(N)*sizeof(BARR_ELTYPE)))            


#define BARR_FREE(BARR)  free(barr)

#define BARR_CLEARARRAY(barr,N)  \
        memset(barr,0,BARR_ARRAYSIZE(N) * sizeof(BARR_ELTYPE))

#define BARR_ELNUM(N)  ((N) / BARR_ELBITS)

#define BARR_BITNUM(N) ((N) % BARR_ELBITS)

#define BARR_SET(barr,N) \
  ((barr)[(BARR_ELNUM(N))] |= (BARR_ELTYPE)1 << BARR_BITNUM(N))
    

#define BARR_CLEAR(barr,N) \
  ((barr)[BARR_ELNUM(N)] &= ~((BARR_ELTYPE)1 << BARR_BITNUM(N)))
    
#define BARR_FLIP(barr,N) \
   ((barr)[(BARR_ELNUM(N)] ^=  (BARR_ELTYPE)1 << BARR_BITNUM(N)))

#define BARR_TEST(barr,N) \
  ((barr)[BARR_ELNUM(N)]  & ((BARR_ELTYPE)1 << BARR_BITNUM(N)))

