//********************  File: PAR-idStoreSD_MatrElem.c  *************/

#include "PAR-shell-model.h"

     /* 
     ** The module entrance function 
     **       storeSD_MatrElem()
     ** calculates and stores diagonal if(diagCalc == NO)
     ** and as many as possible of the nondiagonal SD
     ** matrix elements
     **              <SD'()|OP|SD()>
     */

	    /**** local function declarations ****/

static int twoPartDiagSD_MatrElem(int typeInt, SD_BAS *sdBas, MATR_OP *op_int,
                                                 float *diag_ptr);
     /*
     ** calculates all diagonal matrix elements 
     **           <SD()|OP()| SD()> 
     ** in current Rank process and return number of elements
     ** The two-particle matrix elements of OP() are stored in
     ** op_int->id_diag[] using the formula
     **        op(k,l) = ((num - k) * k)/2 + (l - 1)
     ** where num = 2 * m_orb - 3. 
     ** The result is stores in diag_SD[]. 
     ** The result is stored in memory pointed to by sdStore->freePermMemPtr.
     */

static int threePartDiagSD_MatrElem(SD_BAS *sdBas,MATR_OP *op_veff,
                                                    float *diag_ptr);
     /*
     ** calculates all diagonal matrix elements 
     **           <SD()|op_veff(three_part)()| SD()>
     ** in current Rank process and return number of elements.  
     ** The three-particle matrix elements of OP() are stored in
     ** op_veff->id_diag[] using the formula
     **  op(k,l,m) =  k*(num1 - k*(num2 -k))/6 + (l*(num3 - l))/2
     **             - num_orb + m     
     ** where num1 = 3*num_orb^2 -12*num_orb + 11
     **       num2 = 3*num_orb - 6
     **       num3 = 2*num_orb - 3
     ** The result is stored in memory pointed to by sdStore->freePermMemPtr.
     */

static void nondiagMatrElem(SP_BAS *spBas, SD_BAS *sdBas,
                                           SD_STORE *sdStore);
     /*
     ** calculates and stores on file as many as possible
     ** of the nondiagonal SD matrix elements
     **                <SD()'|OP|SD()>
     ** not previously calculated and stored in memory
     ** The interaction Op is defined through sdStore->typeInt
     */

static void blockOfNondiagMatrElem(SP_BAS *spBas, SD_BAS *sdBas,
				            SD_STORE *sdStore);

     /*
     ** calculates and stores in windMem as many as possible
     ** of the nondiagonal SD matrix elements
     **                <SD()'|OP|SD()>
     ** not previously calculated and stored in memory
     ** The interaction OP is defined through sdStore->typeInt
     */

void printMatrElemData(SD_BAS *sdBas, SD_STORE *sdStore,
                              TID wallTime, TID cpuTime);
     /*
     ** collect for all processes except MASTER
     ** info about precalculation of <SD'|OP|SD>
     ** and send data to outputData() 
     */
void outputData(int rank,SD_STORE *sdStore,int *data);
     /*
     ** send  info about precalculation of
     ** <SD'|OP|SD> to output file in MASTER 
     */

        /**** End: function declarations ****/


// The function definitions

     /* .
     ** The module entrance function 
     **       storeSD_MatrElem()
     ** calculates and stores diagonal 
     **     <SD|OP(typeInt)|SD> 
     ** and nondiag matrix elements
     **     <SD'()|OP|SD()>
     */

void storeSD_MatrElem(SP_BAS *spBas,SD_BAS *sdBas,
                                SD_STORE *sdStore)
{
  // calculate and store <SD|OP(typeInt)|SD> in memory 

  if(sdStore->typeInt == VEFF_INT3) {
    sdStore->numDiagSD = threePartDiagSD_MatrElem(sdBas,&sdStore->op,
                                        (float *)sdStore->diagMemPtr);
  }
  else {
    sdStore->numDiagSD = twoPartDiagSD_MatrElem(ANG_INT, sdBas, &sdStore->op,
			               (float *)sdStore->diagMemPtr);
  }

  // Calculate and store nondiag <SD'|OP|SD> 

  nondiagMatrElem(spBas, sdBas,sdStore);

} // End: function  storeSD_MatrElem()

     /*
     ** The function 
     **         twoPartDiagSD_MatrElem()
     ** calculates all diagonal matrix elements 
     **           <SD()|OP()| SD()> 
     ** in current Rank process and return number of elements 
     ** The two-particle matrix elements of OP() are stored in
     ** op_int->id_diag[] using the formula
     **        op(k,l) = ((num - k) * k)/2 + (l - 1)
     ** where num = 2 * m_orb - 3. 
     ** The result is stores in diag_SD[]. 
     ** The result is stored in memory pointed to by sdStore->freePermMemPtr.
     */

static int twoPartDiagSD_MatrElem(int typeInt, SD_BAS *sdBas, MATR_OP *op_int,
                                                              float *diag_ptr)
{
  char      *func = {"twoPartDiagSD_MatrElem(): "};
   int      num_part, *occ, *occ_ptr, par, orb, num, numSD, 
            loopSD, num_part_1, k, l, *k_ptr,*l_ptr;
   ULL      *sd, sd_value, pos;
   double   *table, *table_ptr, value; 

   table = op_int->id_diag;  /* <|H()|> */
   occ = MALLOC(sdBas->part, int, func, "occ[]"); // local memory

   num_part  = sdBas->part; // initialization
   num_part_1 = num_part - 1;
   num        = (sdBas->numm_orb << 1) - 3;
   sd         = sdBas->SD + Rank;
   numSD      = 0;

   for(loopSD = Rank; loopSD < sdBas->tot_dimSD;
                    loopSD += NumProcs, sd += NumProcs) {
     sd_value = *sd;
     occ_ptr = occ;
     for(par = 0, orb = 0, pos = ULL_ONE; par < num_part; par++,
                                                 orb++, pos <<= 1)  {
       for(;!(sd_value & pos); orb++, pos <<= 1); // particle found
       *occ_ptr = orb; // save orbit
       occ_ptr++;
     }
     value = D_ZERO;  // calculate <SD(Z)|OP(pp)|SD(Z)> */
     k     = num_part_1;
     k_ptr = occ;
     do {    // k-particle loop
       table_ptr = table + (((num - (*k_ptr)) * (*k_ptr)) >> 1) -1;
       l     = k;
       k_ptr++;
       l_ptr = k_ptr;
       do { // l-particle loop */ 
	 value += *(table_ptr + (*l_ptr));
	 l_ptr++;
 
       } while(--l);  // end l particle loop 
     } while(--k);  // end k particle loop
     *diag_ptr = (float) value;
     numSD++;
     diag_ptr++;

   }  // end loop through all |SD()> in current Rank process

   free(occ); // release local memory

   return  numSD;

} // End: function  twoPartDiagSD_MatrElem()

     /*
     ** The function 
     **         threePartDiagSD_MatrElem()
     ** calculates all diagonal matrix elements 
     **           <SD()|op_veff(three_part)()| SD()>
     ** in curren Rank process and  
     ** where the three-particle matrix elements of OP() are stored in
     ** op_veff->id_diag[] using the formula
     **  op(k,l,m) =  k*(num1 - k*(num2 -k))/6 + (l*(num3 - l))/2
     **             - num_orb + m     
     ** where num1 = 3*num_orb^2 -12*num_orb + 11
     **       num2 = 3*num_orb - 6
     **       num3 = 2*num_orb - 3
     ** The result is stored in memory pointed to by sdStore->freePermMemPtr.
     */

static int threePartDiagSD_MatrElem(SD_BAS *sdBas, MATR_OP *op_veff, 
                                                    float *diag_ptr)
{
  char    *func = {"idDiagSD_elem(): "};
   int    num_part, *occ, *occ_ptr, par, orb, numSD, 
          loopSD, num_part_2, num1, num2,num3, num_orb,
          k, l, m,  *k_ptr,*l_ptr, *m_ptr;
   ULL    *sd, sd_value, pos;
   double *table, *table_k_ptr, *table_kl_ptr, value; 

   table = op_veff->id_diag;  // <|H()|>
   occ   = MALLOC(sdBas->part, int, func, "occ[]");

   num_part  = sdBas->part;  // initialization
   num_part_2 = num_part - 2;

   num_orb = sdBas->numm_orb;
   num1    = 3*num_orb*num_orb -12*num_orb + 11;
   num2    = 3*num_orb - 6;
   num3    = 2*num_orb - 3;
   sd      = sdBas->SD + Rank;
   numSD   = 0;

   for(loopSD = Rank; loopSD < sdBas->tot_dimSD;
                    loopSD += NumProcs, sd += NumProcs) {
     sd_value = *sd;
     occ_ptr = occ;
     for(par = 0, orb = 0, pos = ULL_ONE; par < num_part;
                                   par++, orb++, pos <<= 1)  {
       for(;!(sd_value & pos); orb++, pos <<= 1);// particle found
       *occ_ptr = orb; // save orbit
       occ_ptr++;
     }
     value = D_ZERO; // <SD|OP()|SD> = 0
     k     = num_part_2;
     k_ptr = occ;
     do {      // k-particle loop
       table_k_ptr =   table 
                     + ((*k_ptr)*(num1 - (*k_ptr)*(num2 - (*k_ptr))))/6;
       l     = k;
       k_ptr++;
       l_ptr = k_ptr;
       do {    // l-particle loop
	 table_kl_ptr = table_k_ptr + ((*l_ptr)*(num3 - (*l_ptr)))/2;
	 m = l;
         l_ptr++;
	 m_ptr = l_ptr;
	 do {  //  m-particle loop
	   value += *(table_kl_ptr + (*m_ptr - num_orb));
	   m_ptr++;
	 } while(--m); // end m particle loop 
       } while(--l);   // end l particle loop 
     } while(--k);     // end k particle loop
     *diag_ptr = (float) value;
     numSD++;
     diag_ptr++;
   } // end  loop through all |SD()> in current Rank process

   free(occ);  // release local memory

   return numSD;

} // End: function  threePartDiagSD_MatrElem()

     /*
     ** The function
     **       nondiagMatrElem()
     ** calculates and stores in memory
     ** nondiagonal SD matrix elements
     **         <SD()'|OP|SD()>
     ** The interaction Op is defined 
     ** through sdStore->typeInt
     */

static void nondiagMatrElem(SP_BAS *spBas, SD_BAS *sdBas,
                                     SD_STORE *sdStore)
{ 
  char     *func = {"nondiagMatrElem(): "}; 
  STORE_F  *memPtr, *tempPtr;
  TID        wallTime, cpuTime;

  sdStore->storedSD_elem = 0;  // initialization

  // time controll for calculating and store <SD'|OP|SD>

  wallClock(3,0);   // initialization
  cpuClock(3,0);
  wallClock(3,1);    // start clock
  cpuClock(3,1);

  blockOfNondiagMatrElem(spBas, sdBas,sdStore);



/*************   REALLOC  ************

  tempPtr = (STORE_F *)sdStore->tempMemPtr;
  sdStore->tempMemPtr = (char *)REALLOC(tempPtr,numElem,STORE_F,func,
					 "realloc_tempMemory[]");

  if(tempPtr != (STORE_F *)sdStore->tempMemPtr) {
    printf("\n\nRank%d: Error in function nondiagMatrElem():", Rank);
    printf("\nRank%d: REALLOC procedure in error\n",Rank);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }

**************   END REALLOC  *******/

  wallClock(3,2);  // stop clock
  cpuClock(3,2);
  wallTime = wallClock(3,3);  // read clock
  cpuTime  = cpuClock(3,3);

  if(sdStore->typeInt != ANG_INT) {

    printMatrElemData(sdBas,sdStore,wallTime, cpuTime);
  }

} // End: function nondiagMatrElem()

     /*
     ** The function
     **        blockOfNondiagMatrElem()
     ** calculates and stores in temporary memory 
     ** the nondiagonal SD matrix elements
     **                <SD()'|OP|SD()>
     ** The interaction OP is defined through sdStore->typeInt
    */

void blockOfNondiagMatrElem(SP_BAS *spBas, SD_BAS *sdBas,
                                       SD_STORE *sdStore)
{
  HEAD      head;
  STORE_F   *memPtr;

  // calculate and store in temporary memory <SD'|typeInt|SD>

  head.memSize = (int)(sdStore->tempMemSize/sizeof(STORE_F));
  head.memPtr  =(STORE_F *)sdStore->tempMemPtr;// start of data field
  head.numSD   = 0;

  head.typeSD  = (sdBas->MJ == 0) ? 0 : 2;
  head.startSD = Rank;
  head.endSD   = sdBas->numSD[0];

  if(sdStore->typeInt == VEFF_INT3) {
    threePartNondiagSD_Calc(&head,spBas,sdBas,&sdStore->op);
  }
  else {
    twoPartNondiagSD_Calc(&head,spBas,sdBas,&sdStore->op);
  }
  sdStore->last_asymSD    = head.lastSD; // last calculated asym |initSD>
  sdStore->lastSD         = sdStore->last_asymSD;
  sdStore->storedSD_elem += head.numSD;

  if(sdBas->numSD[1])  {  // <SD'|typeInt|SD_ sym>
    head.memSize -= head.numSD;
    head.memPtr   = (STORE_F *)sdStore->tempMemPtr + head.numSD;
    head.typeSD   = 1;
    head.startSD  = head.lastSD + NumProcs;
    head.endSD    = sdBas->tot_dimSD;
    head.numSD    = 0;
    if(sdStore->typeInt == VEFF_INT3) {
      threePartNondiagSD_Calc(&head,spBas,sdBas,&sdStore->op);
    }
    else {
      twoPartNondiagSD_Calc(&head,spBas,sdBas,&sdStore->op);
    }
    sdStore->lastSD         = head.lastSD;
    sdStore->storedSD_elem += head.numSD;

  }   // end sym case

  return;
    
} // End: function  blockOfNondiagMatrElem()
  
  /*
  ** The function
  **  printMatrElemData()
  ** collect for all processes except MASTER
  ** info about precalculation of <SD'|OP|SD>
  ** and send data to outputData() 
  */

void printMatrElemData(SD_BAS *sdBas,SD_STORE *sdStore,
                       TID wall_time, TID cpu_time)
{ 
  int   data[9];
  
  data[0]  = (int)sdStore->storedSD_elem;
  data[1]  = sdStore->lastSD;
  data[2]  = sdBas->tot_dimSD;
  data[3]  = (int) wall_time.hour;
  data[4]  = (int) wall_time.min;
  data[5]  = (int) wall_time.sec;
  data[6]  = (int) cpu_time.hour;
  data[7]  = (int) cpu_time.min;
  data[8]  = (int) cpu_time.sec;
 
  if(Rank == MASTER) {
    int         receive;
    MPI_Status  status;

    outputData(MASTER,sdStore,data);


/*****************   Only Rank = MASTER **********

    for(receive = 1; receive < NumProcs; receive++) {
      MPI_Recv(data, 15, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG,
                                     MPI_COMM_WORLD, &status);
      outputData(status.MPI_SOURCE,sdStore,data);
    }

*************   end output  ****/ 

  } // end MASTER

/*****************************

  else {
    MPI_Send(data,15,MPI_INT,MASTER,Rank,MPI_COMM_WORLD);
  } // end slave nodes 
*********************************/

} // End: function printMatrElemData()

  /*
  ** The function 
  **      outputData()
  ** send  info about precalculation of
  ** <SD'|OP|SD> to output file in MASTER 
  */

void outputData(int rank,SD_STORE *sdStore, int *data)
{
  char  text[5];
  int   num;
  FILE  *file_ptr;

  switch(sdStore->typeInt) {
      case  VEFF_INT:  sprintf(text,"%s","Veff");
                       break;
      case  VEFF_INT3: sprintf(text,"%s","Veff");
                       break;
      case  ANG_INT:   sprintf(text,"%s","Ang");
	               break;
       case  CM_INT:   sprintf(text,"%s","CM");
                       break;
  } // end switch

  if((file_ptr = fopen(sdStore->result,"a"))== NULL) {
    printf("\n\nRank(%d): Error in function outputData():", rank);
    printf("\nWrong file = %s for the output data file\n", sdStore->result);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }
  if(sdStore->singleDataBlock == YES) {
    fprintf(file_ptr,"\n\nRank%d Stored nondiag matrix elements in memory",rank);
  }
  else {
    fprintf(file_ptr,"\n\nRank%d Stored nondiag matrix elements in memory",rank);
  }
  fprintf(file_ptr,"\nNumber of stored <SD'|%s|initSD>: %d",text,data[0]);
  fprintf(file_ptr,"\nUp to and including |initSD> = %d",data[1]); 
  fprintf(file_ptr," - Last |initSD> =  %d",data[2]);
  
  num = (text,data[0] * sizeof(STORE_F))/1000; 
  fprintf(file_ptr,"\n Memory requirement = %d KB",num);
  fprintf(file_ptr,"\nWall time used : %d hour %d min %d sec",
	                             data[3],data[4],data[5]);
  fprintf(file_ptr, "\nCPU time used : %d hour %d min %d sec",
            	                      data[6],data[7],data[8]);
  fclose(file_ptr);
  
} // End: function outputData()
