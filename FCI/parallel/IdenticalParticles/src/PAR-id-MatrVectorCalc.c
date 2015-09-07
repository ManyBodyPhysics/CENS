/****************  The module  PAR-id-MatrVectorCalc.c  ******************/

#include "PAR-shell-model.h"

   /*
   ** The entrance function
   **         id_matrix_vector_calc()     
   ** takes as input an sdBlock[Rank] orthonormalized Lanczos
   ** vector and performs the process
   **   OP(typeInt) * |initVec[sdBlock[Rank]]> --> |finalVec[]>
   ** OP(typeInt) is an effective two- or three-particle operator 
   ** and its matrix elements is stored in the table MATR_OP typeInt[].
   ** The process contains up to four steps;
   **  1. Contribution from diag <SD|VEFF|SD>, |<
   **     a) if(typeInt == VEFF) stored in memory
   **     b) If(typeInt == (ANG or CM) stored on file
   **  2. Contribution from nondiag <SD'|VEFF|SD> 
   **     stored in memory lancProc->memNondiagPtr 
   **     up to memory size determined by 
   **           lancProc->permMemSize
   **  3. If step 2 does not complete the calculation
   **     add contribution from nondiag <SD'|OP(typeInt)|SD> 
   **     stored on files up to file size determined by 
   **           lancProc->file_nondiag,
   **  4. If step 3 does not complete the calculation
   **     add contribution from nondiag <SD'|OP(typeInt)|SD>
   **     by explicit calculation
   */

 /****  local function declarations ****/

static void idDiagMemSD_Contr(int numDiagSD, SD_STORE *sdStore, 
                           double *initVec, double *finalVec);

   /*
   ** calculates and store the contribution to 
   **       H(OP) * |initVec> ---> |finalVec>
   ** from the precalculated diagonal matrix elements
   **  stored as float values in memory
   */

static nondiagMemContr(SD_STORE *sdStore, double *initVec,
                                              double *finalVec);
    /*
    ** calculates contribution from all precalculated 
    ** non-diagonal matrix elements 
    **       <SD'()|H|SD()>
    ** stored in memory to produce a new  |finalVec>
    ** For MJ = 0: |SD()> = |SD_asym> followed by |SD_sym>
    ** For MJ > 0 : |SD()> = |SD_nosym> only 
    */

static STORE_F *asymContr(int MJ,int initSD,int numInitSDasym,STORE_F *memPtr,
                                      double *initVec,double *finalVec);
    /*
    ** calculates the contribution H * |initVec> ---> |finalVec>
    ** from a single block of precalculates <finalSD'|H|initSD>
    ** if( MJ == 0) initSD = |initSD_asym>
    ** if( MJ > 0) initSD = |initSD_nosym>
    ** return address to next data element
    */

static STORE_F *symContr(int initSD,int numInitSDsym,STORE_F *memPtr,
                            double *initVec,double *finalVec);
    /*
    ** calculates the contribution 
    ** H * |initVec_sym> ---> |finalVec>
    ** from a single block of precalculates 
    ** <finalSD'|H|initSD_sym>
    */

               /****   end local function declarations ****/

               /**** The function definitions  ****/ 
 
   /*
   ** The entrance function
   **         id_matrix_vector_calc()     
   ** takes as input an sdBlock[Rank] orthonormalized Lanczos
   ** vector and performs the process
   **   OP(typeInt) * |initVec[sdBlock[Rank]]> --> |finalVec[]>
   ** OP(typeInt) is an effective two- or three-particle operator 
   ** and its matrix elements is stored in the table MATR_OP typeInt[].
   ** The process contains up to four steps;
   **  1. Contribution from diag <SD|VEFF|SD>, 
   **     a) if(typeInt == VEFF) stored in memory
   **     b) If(typeInt == (ANG or CM) stored on file
   **  2. Contribution from nondiagem <SD'|VEFF|SD> 
   **     stored in memory lancProc->memNondiagPtr 
   **     up to memory size determined by 
   **           lancProc->permMemSize
   **  3. If step 2 does not complete the calculation
   **     add contribution from nondiag <SD'|OP(typeInt)|SD> 
   **     stored on files up to file size determined by 
   **           lancProc->file_nondiag,
   **  4. If step 3 does not complete the calculation
   **     add contribution from nondiag <SD'|OP(typeInt)|SD>
   **     by explicit calculation
   */

void  id_matrix_vector_calc(SP_BAS *spBas,SD_BAS *sdBas, SD_STORE *sdStore,
                                          double *initVec, double *finalVec)
{
  int   nextInitSD, k;




  idDiagMemSD_Contr(sdStore->numDiagSD, sdStore,initVec,finalVec);



    // All <SD'|OP|SD> stored in memory
  
    nondiagMemContr(sdStore, initVec, finalVec);


} // End: function  id_matrix_vector_veff_calc()

     /*
     **        idDiagMemSD_Contr()
     ** calculates contribution to 
     **       H(op) * |initVec> ---> |initVec>
     ** from the precalculated diagonal matrix elements
     **  stored as float values in memory
     */

static void idDiagMemSD_Contr(int numDiagSD,SD_STORE *sdStore,
                             double *initVec, double *finalVec)
{
  int       k;
  double    *initPtr,*finalPtr;
  float     *memPtr;
  

  initPtr  = initVec  + Rank;
  finalPtr = finalVec + Rank;

  memPtr = (float *)sdStore->diagMemPtr;

  for(k = 0; k < numDiagSD; k++,memPtr++) {
    *finalPtr = ((double)(*memPtr))*(*initPtr);
    initPtr  += NumProcs;
    finalPtr += NumProcs;
  }

} // End: function idDiagMemSD_Contr() 

     /*
     ** The function 
     **         nondiagMemContr()
     ** calculates contribution from all precalculated 
     ** non-diagonal matrix elements 
     **       <SD'()|H|SD()>
     ** stored in memory to produce a new  |finalVec>
     ** For MJ = 0: |SD()> = |SD_asym> followed by |SD_sym>
     ** For MJ > 0 : |SD()> = |SD_nosym> only 
     */

static nondiagMemContr(SD_STORE *sdStore,double *initVec, double *finalVec)
{
  int     initSD;
  ULL     topBit;
  double  *initVecPtr, *finalVecPtr, initAmpl;
  STORE_F *memPtr, *localPtr;

  topBit = ULL_ONE << (8 * sizeof(ULL) - 1);

  memPtr      = (STORE_F *)sdStore->tempMemPtr;
  initVecPtr  = initVec  + Rank;
  finalVecPtr = finalVec + Rank;

  for(initSD = Rank; initSD <= sdStore->last_asymSD; initSD += NumProcs,
                    initVecPtr += NumProcs, finalVecPtr += NumProcs) { 
    localPtr = memPtr;
    initAmpl = *initVecPtr;

    if(sdStore->MJ == 0) {

      while(!(localPtr->final & topBit)) {
	switch((localPtr->final) & 3)   {
           case 0:
           case 3: finalVec[localPtr->final >> 2] 
                      += ((double)localPtr->value)*initAmpl;
                   break;
	   case 1: finalVec[localPtr->final >> 2] 
                     +=  2*((double)localPtr->value)*initAmpl;
                   break;
	   case 2:
                   break;
         }  // end switch()

	localPtr++;

      } // end while()

      while(!(memPtr->final & topBit)) {
	*finalVecPtr += ((double)memPtr->value)*initVec[memPtr->final>>2];
	memPtr++;
      } // end while()
       memPtr++;
    } // end MJ = 0
 
    else {  // MJ > 0

      while(!(localPtr->final & topBit)) {
	finalVec[localPtr->final] 
	       += (double)localPtr->value *initAmpl;
	localPtr++;
      } // end while() 
      while(!(memPtr->final & topBit)) {
	*finalVecPtr += ((double)memPtr->value)*initVec[memPtr->final];
	memPtr++;
      } // end while()
      memPtr++;
    } // end MJ > 0

  }  // end loop through all precalculated nondiag asym matrix elements

  if(sdStore->MJ > 0) return;  // no sym-case


  for( ; initSD <= sdStore->lastSD;initSD += NumProcs,
        initVecPtr += NumProcs, finalVecPtr += NumProcs) { 
    localPtr = memPtr;
    initAmpl = *initVecPtr;
    while(!(localPtr->final & topBit)) {
      if(((localPtr->final) & 3) == 0)  {
	finalVec[localPtr->final >> 2] 
          += ((double)localPtr->value)*initAmpl;
      }
      localPtr++;
    } // end while()

    while(!(memPtr->final & topBit)) {
      switch((memPtr->final) & 3)   {
	 case 0:
	 case 1: 
                *finalVecPtr += 
                    ((double)memPtr->value)*initVec[memPtr->final >> 2];
                 break;
      }  // end switch() loop

      memPtr++;
    } // end while() loop
    memPtr++;
  }  // end loop through all precalculated nondiag sym matrix elements
  
} // End: funtion nondiagMemContr()

     /*
     ** The function 
     **       asymContr()
     ** calculates the contribution H * |initVec> ---> |finalVec>
     ** from a single block of precalculates <finalSD'|H|initSD>
     ** if( MJ == 0) initSD = |initSD_asym>
     ** if( MJ > 0) initSD = |initSD_nosym>
     ** return address to next data element
     */

static STORE_F *asymContr(int MJ,int initSD, int nextInitSDasym,STORE_F *memPtr,
                                   double *initVec,double *finalVec)
{
  double   initAmpl, *initVecPtr, *finalVecPtr;
  ULL      topBit;
  STORE_F  *localPtr;

  topBit = ULL_ONE << (8 * sizeof(ULL) - 1);

  initVecPtr  = initVec  + initSD;
  finalVecPtr = finalVec + initSD;

  for( ; initSD < nextInitSDasym;initSD += NumProcs,
	 initVecPtr += NumProcs, finalVecPtr += NumProcs) {
    initAmpl = *initVecPtr;
    localPtr = memPtr;

    if(MJ == 0) {
      while(!(localPtr->final & topBit)) {
         switch((localPtr->final) & 3)   {
           case 0:  
           case 3: finalVec[localPtr->final >> 2] 
		    += ((double)localPtr->value)*initAmpl;
	           break;
           case 1: finalVec[localPtr->final >> 2] 
                     += 2 * ((double)localPtr->value)*initAmpl;
                  break;
           case 2: break;
         }  // end switch()

       localPtr++;
      } // end while() loop
      while(!(memPtr->final & topBit)) {
	*finalVecPtr += ((double)memPtr->value) * initVec[memPtr->final >> 2];

	memPtr++;
      } // end while()
      memPtr++;
    } // end loop for MJ == 0

    else {  // nosym case 
      while(!(localPtr->final & topBit)) {
	finalVec[localPtr->final] += ((double)localPtr->value) *initAmpl;
	localPtr++;
      } // end while() loop
      while(!(memPtr->final & topBit)) {
	*finalVecPtr += ((double)memPtr->value)*initVec[memPtr->final];
	memPtr++;
      } // end while() loop
      memPtr++;
    } // end loop for MJ > 0

  } // end contr from all precalc |asymSD()> in a mem block

  return memPtr;

} // End: function asymContr()

     /*
     ** The function 
     **       symContr()
     ** calculates the contribution 
     ** H * |initVec_sym> ---> |finalVec>
     ** from a single block of precalculates 
     ** <finalSD'|H|initSD_sym>
     ** return address to next data element
     */

static STORE_F *symContr(int initSD, int nextInitSDsym,STORE_F *memPtr,
                   double *initVec,double *finalVec)
{
  int      count;
  double   initAmpl, *initVecPtr, *finalVecPtr;
  ULL      topBit;
  STORE_F  *localPtr;

  topBit = ULL_ONE << (8 * sizeof(ULL) - 1);

  initVecPtr  = initVec  + initSD;
  finalVecPtr = finalVec + initSD;

  for( ; initSD < nextInitSDsym;initSD += NumProcs,
         initVecPtr += NumProcs, finalVecPtr += NumProcs) {
    initAmpl = *initVecPtr;
    localPtr = memPtr;

    while(!(localPtr->final & topBit)){
      if(((localPtr->final) & 3) == 0) {
	finalVec[(localPtr->final >> 2)] 
	  += ((double)localPtr->value)*initAmpl;
      }

      count++;

      localPtr++;
    } // end while()

    while(!(memPtr->final & topBit)) {
      switch((memPtr->final) & 3)   {
	    case 0: *finalVecPtr 
                     += ((double)memPtr->value)*initVec[memPtr->final >> 2];
                    break;
	    case 1: *finalVecPtr
                     += ((double)memPtr->value)*initVec[memPtr->final >> 2];
                    break;
      } // end switch() loop

      memPtr++;
    } /* end while() loop */
    memPtr++;
  } // end contr from all precalc |symSD()> in a mem block

  return memPtr;

} // End: function symContr()
