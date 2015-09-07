/****************  The module  id-matrix-vector-calc.c  ******************/

#include "shell-model.h"

   /*
   ** The entrance function
   **         id_matrix_vector_calc id_matrix_vector_veff_calc()     
   ** takes as input an orthonormalized Lanczos vector
   ** and performs the process
   **            H(veff) * |initVec[]> = |finalVec[]>
   ** H(veff) is a two- or three-particle operator and its matrix
   ** elements is stored in the table MATR_OP veff[].
   ** The process contains up to four steps;
   **  1. Contribution from diag <SD|VEFF|SD> 
   **     allways stored in memory
   **  2. Contribution from nondiag <SD'|VEFF|SD> 
   **     stored in memory lancProc->memNondiagPtr 
   **     up to memory size determined by 
   **           lancProc->permMemSize
   **  3. If step 2 does not complete the calculation
   **     add contribution from nondiag <SD'|VEFF|SD> 
   **     stored on files up to file size determined by 
   **           lancProc->file_nondiag,
   **  4. If step 3 does not complete the calculation
   **     add contribution from nondiag <SD'|VEFF|SD>
   **     by explicit calculation
   */

               /****  local function declarations ****/

static void idDiagMemSD_Contr(SD_STORE *sdStore, SD_BAS *sdBas, 
                             double *initVec, double *finalVec);
   /*
   ** calculates and store the contribution to 
   **       H * |initVec> ---> |finalVec>
   ** from the precalculated diagonal matrix elements
   **  stored as float values in memory
   */

static void nondiagMemContr(SD_STORE *sdStore, SD_BAS *sdBas, 
		     double *initVec, double *finalVec);
   /*
    ** calculates contribution to
    **    H(veff) * |vec.one[]> = |vec_two[]>
    ** from precalculated nondiag <SD'|VEFF|SD> stored
    ** in memory.
    */

static void idDiagFileSD_Contr(SD_STORE *sd_store,
	       SD_BAS *sdBas,double *initVec,  double *finalVec);
    /*
    ** calculates diagonal contributions to 
    **     <finalVec|OP(typeInt) |initVec>
    ** from the precalculated diagonal matrix elements
    **  stored as float values on file 
    */

static void nondiagFileContr(SD_STORE *sdStore, SD_BAS *sdBas, 
		      double *initVec, double *finalVec);
    /*
    ** calculates contribution to the matrix-vector
    **    H(VEFF) * |vec.one[]> = |vec_two[]>
    ** from precalculated nondiag <SD'|VEFF|SD> stored 
    ** on files
    */

static void nondiagExplicitCalc(int lastSD_asym, int lastSD_sym,
                        SD_STORE *sdStore, SP_BAS *spBas, SD_BAS *sdBas, 
				double *initVec, double *finalVec);
    /*
    ** calculates contribution to the matrix-vector 
    **    H(veff) * |vec.one[]> = |vec_two[]>
    ** from the precalculated nondiag <SD'|VEFF|SD> stored 
    ** on files
    */

static void  idNondiagSymMemContr(SD_BAS *sdBas, SD_STORE *sdStore, 
				       double *initVec, double *finalVec);
    /*
    ** calculates contribution from all precalculated 
    ** non-diagonal matrix elements of  
    **  <SD'()|H|asymSD()>  and  <SD'()|H|symSD()
    ** stored in memory to the new lanczo vector |finalVec>.
    */

static void idNondiagNosymMemContr(SD_BAS *sdBas, SD_STORE *sdStore,
				       double *initVec, double *finalVec);
    /*
    ** calculates contribution  to the new lanczo
    ** from all precalculated non-diagonal matrix 
    ** elements <SD'()|H|nosymSD()> stored in 
    ** permanent memory 
    */

static void idNondiagAsymSD_FileContr(SD_BAS *sdBas, int *sect_file_no, 
                   SD_STORE *sdStore, double *initVec, double *finalVec);
    /*
    ** calculates contribution from all precalculated 
    ** non-diagonal matrix elements 
    **                 <SD'()|H|asymSD()>
    ** stored on file to the new lanczo vector |finalVec>.
    ** The function returns the number of last initSD.
    */

static void idNondiagSymSD_FileContr(SD_BAS *sdBas, int *sect_file_no, 
		  SD_STORE *sdStpre, double *initVec, double *finalVec);
    /*
    ** calculates contribution from all precalculated 
    ** non-diagonal matrix elements 
    **                 <SD'()|H|asymSD()>
    ** to the new lanczo vector |finalVec>.
    ** The function returns the number of last initSD.
    */

static void idNondiagNosymSD_FileContr(SD_BAS *sdBas, int *sect_file_no, 
                   SD_STORE *sdStore, double *initVec, double *finalVec);
    /*
    ** calculates contribution from all precalculated 
    ** non-diagonal matrix elements 
    **                 <SD'()|H|nosymSD()>
    ** to the new lanczo vector |finalVec>.
    ** The function returns the number of last initSD.
    */

static void id_nondiag_asymSD_calc_contr(int initSD, SP_BAS *spBas, 
                                          SD_BAS *sdBas, SD_STORE *sdStore, 
                                           double *initVec, double *finalVec);
   /*
   ** calculates the contribution from the non-precalculated 
   ** non-diagonal matrix elements 
   **                <SD'()|H|asymSD()>
   ** to the new lanczo vector |vec[2]> 
   */

static void id_nondiag_symSD_calc_contr(int initSD, SP_BAS *spBas,
                                       SD_BAS *sdBas, SD_STORE *sdStore,
				     double *initVec, double *finalVec);
   /*
   ** calculates the contribution from the non-precalculated 
   ** non-diagonal matrix elements 
   **                 <SD'()|H|asymSD()>
   ** to the new lanczo vector |finalVec> 
   */

static void id_nondiag_nosymSD_calc_contr(int initSD, SP_BAS *spBas,
                                               SD_BAS *sdBas, SD_STORE *sdStore, 
					  double *initVec, double *finalVec);
   /*
   ** calculates the contribution from the non-precalculated 
   ** non-diagonal matrix elements 
   **                 <SD'()|H|nosymSD()>
   ** to the new lanczo vector |finalVec> 
   */
               /****   end local function declarations ****/

               /**** The function definitions  ****/ 
 
   /*
   ** The entrance function
   **         id_matrix_vector_calc()     
   ** takes as input an orthonormalized Lanczos vector
   ** and performs the process
   **          H(veff) * |initVec[]> = |finalVec[]>
   ** H(veff) is an effective two- or three-particle operator 
   ** and its matrix elements is stored in the table MATR_OP veff[].
   ** The process contains up to four steps;
   **  1. Contribution from diag <SD|VEFF|SD> 
   **     stored in memory
   **  2. Contribution from nondiag <SD'|VEFF|SD> 
   **     stored in memory lancProc->memNondiagPtr 
   **     up to memory size determined by 
   **           lancProc->permMemSize
   **  3. If step 2 does not complete the calculation
   **     add contribution from nondiag <SD'|VEFF|SD> 
   **     stored on files up to file size determined by 
   **           lancProc->file_nondiag,
   **  4. If step 3 does not complete the calculation
   **     add contribution from nondiag <SD'|VEFF|SD>
   **     by explicit calculation
   */

void  id_matrix_vector_calc(SP_BAS *spBas, SD_BAS *sdBas,
             SD_STORE *sdStore, double *initVec, double *finalVec)
{
  int  lastSD_asym, lastSD_sym;

  if(sdStore->diagMemCalc == YES) {
    idDiagMemSD_Contr(sdStore, sdBas, initVec, finalVec);
  }

  if(  (sdStore->lastSD_mem_asym > -1)    
     ||(sdStore->lastSD_mem_sym  > -1)) {

    // <SD'|OP|SD> stored in memory

    nondiagMemContr(sdStore, sdBas, initVec, finalVec);

  } // end contribution from memory

  if(sdStore->diagMemCalc == NO) {
    idDiagFileSD_Contr(sdStore,sdBas,initVec,finalVec);
  }

  if(  (sdStore->lastSD_file_asym > -1)
     ||(sdStore->lastSD_file_sym  > -1)) {
	
    //  <SD'|OP|SD> stored on files

    nondiagFileContr(sdStore, sdBas, initVec, finalVec);

  }  // end contribution from file

  lastSD_asym = MAX(sdStore->lastSD_mem_asym,
                    sdStore->lastSD_file_asym); 
  lastSD_sym  = MAX(sdStore->lastSD_mem_sym,
                    sdStore->lastSD_file_sym);

  if(  ((sdBas->numSD[0] > 0)&&(lastSD_asym < (sdBas->numSD[0] - 1)))
     ||((sdBas->numSD[1] > 0)&&(lastSD_sym < (sdBas->tot_dimSD - 1)))) {

    //  explicit calculation of <SD'|VEFF|SD> 

    nondiagExplicitCalc(lastSD_asym, lastSD_sym,
                sdStore, spBas, sdBas, initVec, finalVec);

  } // end contribution explicit calculation

 } // End: function  id_matrix_vector_veff_calc()

     /* The function
     **        idDiagMemSD_Contr()
     ** calculates contribution to 
     **       H * |initVec> ---> |initVec>
     ** from the precalculated diagonal matrix elements
     **  stored as float values in memory
     */

static void idDiagMemSD_Contr(SD_STORE  *sdStore, SD_BAS *sdBas,
                              double *initVec, double *finalVec)
{
  int       k;
  float     *mem_ptr;

  mem_ptr = (float *)sdStore->permMemPtr;
  for(k = 0; k< sdBas->tot_dimSD; k++) { 
    *(finalVec++) = (double)(*(mem_ptr++)) * (*(initVec++));
  }

} // End: function idDiagMemSD_Contr() 

     /*
     ** The function 
     **      nondiagMemContr()
     ** calculates contribution to
     **    OP * |vec.one[]> = |vec_two[]>
     ** from precalculated nondiag <SD'|OP|SD> stored
     ** in memory.
     */

static void nondiagMemContr(SD_STORE *sdStore, SD_BAS *sdBas, 
		    double *initVec, double *finalVec)
{
  if(sdBas->MJ == 0) {
    idNondiagSymMemContr(sdBas, sdStore,
                           initVec, finalVec);
  }
  else {
    idNondiagNosymMemContr(sdBas, sdStore,
                            initVec, finalVec);
  }
  
} // End function nondiagMemContr()

     /*
     ** The function 
     **        idDiagFileSD_Contr()
     ** calculates diagonal contributions to 
     **     <finalVec|OP(typeInt) |initVec>
     ** from the precalculated diagonal matrix elements
     **  stored as float values on file 
     */

static void idDiagFileSD_Contr(SD_STORE *sdStore,
       SD_BAS *sdBas,double *initVec,  double *finalVec)
{
  char    sect_file[ONE_LINE];
  int     k, header[2], num_elem;
  float   *mem_ptr;
  FILE    *file_ptr;

  sprintf(sect_file,"%s%s",sdStore->title, DIAG);

  if((file_ptr = fopen(sect_file,"rb")) == NULL) {
    printf("\n\nError in function idDiagFileSD_Contr():");
    printf("\nNot allowed to open the file %s\n", sect_file);
    exit(1);
  }
  rewind(file_ptr);

  if(fread((void *) &header,(size_t) sizeof(int),2,file_ptr) != 2)  {
    printf("\n\nError in function id_diag_file__SD_contr():");
    printf("\nIn reading diagonal header elements from file %s\n",
                                                        sect_file);
    exit(1);
  }
  num_elem = header[0] + header[1];

  if(num_elem != sdBas->tot_dimSD) {
    printf("\n\nError in function idDiagFileSD_Contr():");
    printf("\nNumber of basis |SD()> is inconsistent");
    printf("\nheader[0] = %d  header[1] = %d",
                         header[0], header[1]);
    printf("\n should be equal to numSD[0] = %d",
                                sdBas->numSD[0]);
    printf(" and numSD[1] = %d in file %s\n",
                   sdBas->numSD[1], sect_file);
    exit(1);
  } 

  // read diag matrix elements  from file

  mem_ptr = (float *)sdStore->windMemPtr;
  if(fread((void *)mem_ptr,(size_t) sizeof(float),num_elem, file_ptr)
                                                         != num_elem)  {
    printf("\n\nError in function idDiagFileSD_Contr():");
    printf("\nIn reading %d diag matrix elements from  file %s\n",
                                              num_elem, sect_file);
    exit(1);
  }
  for(k = 0; k < sdBas->numSD[0]; k++,initVec++, finalVec++) {
    *finalVec = ((double)(*mem_ptr++)) * (*initVec);
  }
  for(k = 0; k < sdBas->numSD[1]; k++,initVec++, finalVec++) {
    *finalVec = ((double)(*mem_ptr++)) * (*initVec);
  }

  fclose(file_ptr);  // close sect_file

} // End: function idDiagFileSD_Contr()

     /*
     ** The function 
     **      nondiagFileContri()
     ** calculates contribution to the matrix-vector
     **    H(VEFF) * |vec.one[]> = |vec_two[]>
     ** from precalculated nondiag <SD'|VEFF|SD> stored 
     ** on files
     */

static void nondiagFileContr(SD_STORE *sdStore,SD_BAS *sdBas, 
                             double *initVec, double *finalVec)
{
  int   sect_file_no;

  sect_file_no = 1; // initialize first nondiag section file number

  if((sdBas->MJ == 0) && (sdStore->lastSD_file_asym > -1)) {
    idNondiagAsymSD_FileContr(sdBas, &sect_file_no, sdStore,
                                              initVec, finalVec);
  } // end asym contribution

  if((sdBas->MJ == 0) && (sdStore->lastSD_file_sym > -1)) {
    idNondiagSymSD_FileContr(sdBas, &sect_file_no, sdStore,
                                             initVec, finalVec);
  } // end sym contribution

  if((sdBas->MJ != 0) && (sdStore->lastSD_file_asym > -1)) {
    idNondiagNosymSD_FileContr(sdBas, &sect_file_no, sdStore,
                                               initVec, finalVec);
  } // end nosym contribution

} // End: function  nondiagFileCalc()

     /*
     ** The function 
     **      nondiagExplicitCal()
     ** calculates explicitly contribution to the matrix-vector 
     **    H(veff) * |vec.one[]> --> |vec_two[]>
     ** from not precalculated nondiag <SD'|VEFF|SD>
     */

static void nondiagExplicitCalc(int lastSD_asym, int lastSD_sym,
                 SD_STORE *sdStore, SP_BAS *spBas, SD_BAS *sdBas, 
                                 double *initVec, double *finalVec)
{
  int initSD;

  if(   (sdBas->MJ == 0)
     && (lastSD_asym < (sdBas->numSD[0] - 1))) {

    initSD = lastSD_asym + 1;
    id_nondiag_asymSD_calc_contr(initSD, spBas,sdBas, sdStore,
                                               initVec,finalVec);
  } // end time-reversal asym contribution

  if(   (sdBas->MJ == 0)
     && (lastSD_sym < (sdBas->tot_dimSD - 1))) {
    initSD  = lastSD_sym + 1;
    id_nondiag_symSD_calc_contr(initSD, spBas, sdBas,
                                  sdStore, initVec, finalVec);
  }  // end time-reversal sym contribution 

 if(   (sdBas->MJ != 0)
      &&(lastSD_asym < (sdBas->numSD[0] - 1))) {
    initSD = lastSD_asym + 1;

    id_nondiag_nosymSD_calc_contr(initSD, spBas,sdBas, sdStore,
                                                 initVec,finalVec);
  }  // end no time-reversal symmetry contribution
  
} // End: function  nondiagExplicitCalc()

     /*
     ** The function 
     **         idNondiagSymMemContr()
     ** calculates contribution from all precalculated 
     ** non-diagonal matrix elements of  
     **  <SD'()|OP|asymSD()>  and  <SD'()|OP|symSD()
     ** stored in memory to the new lanczo vector |finalVec>.
     */

static void  idNondiagSymMemContr(SD_BAS *sdBas, SD_STORE *sdStore, 
                                        double *initVec, double *finalVec)
{
  int     initSD;
  ULL     topBit;
  double  *finalPos, initAmpl;
  STORE_F *memPtr, *localPtr;

  topBit = ULL_ONE << (8 * sizeof(ULL) - 1);
  memPtr = (STORE_F *)(sdStore->permMemPtr + (sdBas->tot_dimSD*sizeof(float)));
  for(initSD = 0; initSD <= sdStore->lastSD_mem_asym; initSD++) { 
    localPtr = memPtr;
    initAmpl = initVec[initSD];
    while(!(localPtr->final & topBit)) {
      switch((localPtr->final) & 3)   {
         case 0:
         case 3: finalVec[localPtr->final >> 2] 
                    += ((double)localPtr->value)*initAmpl;
                 break;
	 case 1: finalVec[localPtr->final >> 2] 
                    += 2*((double)localPtr->value)*initAmpl;
                 break;
	 case 2: break;
      }  // end switch()
      localPtr++;
    } // end while()
    finalPos = &finalVec[initSD];
    while(!(memPtr->final & topBit)) {
      *finalPos += ((double)memPtr->value)*initVec[memPtr->final >> 2];
      memPtr++;
    } // end while()
    memPtr++;
  }  // end loop through all precalculated nondiag asym matrix elements

  for(initSD = sdBas->numSD[0]; initSD <= sdStore->lastSD_mem_sym; 
                                                               initSD++) { 
    localPtr = memPtr;
    initAmpl = initVec[initSD];
    while(!(localPtr->final & topBit)) {
      if(((localPtr->final) & 3) == 0)  {
	finalVec[localPtr->final >> 2]
          += ((double)localPtr->value)*initAmpl;
      }
      localPtr++;
    } // end while()
    finalPos = &finalVec[initSD];
    while(!(memPtr->final & topBit)) {
      switch((memPtr->final) & 3)   {
	 case 0:
	 case 1: *finalPos += ((double)memPtr->value)
                             *initVec[memPtr->final >> 2];
                 break;
      }  // end switch() loop
      memPtr++;
    } // end while() loop
    memPtr++;
  }  // end loop through all precalculated nondiag sym matrix elements

} // End: funtion idNondiagSymMemContr()

     /*
     ** The function 
     **         idNondiagNosymMemContr()
     ** calculates contribution  to the new lanczo
     ** from all precalculated non-diagonal matrix 
     ** elements <SD'()|H|nosymSD()> stored in 
     ** permanent memory 
     */

static void idNondiagNosymMemContr(SD_BAS *sdBas, SD_STORE *sdStore,
                                   double *initVec, double *finalVec)
{
  int       initSD;
  ULL       topBit;
  double    *finalPos, initAmpl;
  STORE_F   *memPtr, *localPtr;

  topBit = ULL_ONE << (8 * sizeof(ULL) - 1);
  memPtr = (STORE_F *)(sdStore->permMemPtr + (sdBas->tot_dimSD*sizeof(float))); 

  for(initSD = 0; initSD <= sdStore->lastSD_mem_asym; initSD++) {
    localPtr = memPtr;
    initAmpl = initVec[initSD];
    while(!(localPtr->final & topBit)) {
      finalVec[localPtr->final] += (double)localPtr->value * initAmpl;
      localPtr++;
    } // end while() 
    finalPos = &finalVec[initSD];
    while(!(memPtr->final & topBit)) {
      *finalPos += (double)memPtr->value*initVec[memPtr->final];
      memPtr++;
    } // end while()
    memPtr++;
  }  // end loop through all precalculated nondiag nosym matrix elements

} // End: function idNondiagNosymMemContr()

     /*
     ** The function 
     **         idNondiagAsymSD_FileContr()
     ** calculates contribution from all precalculated 
     ** non-diagonal matrix elements 
     **                 <SD'()|H|asymSD()>
     ** stored on file to the new lanczo vector |finalVec>.
     ** The function returns the number of last initSD.
     */

static void idNondiagAsymSD_FileContr(SD_BAS *sdBas,
                                int *sectFileNo, SD_STORE *sdStore,
                                double *initVec, double *finalVec)
{
  char     basFile[ONE_LINE], sectFile[ONE_LINE];
  int      lastSectFileNo, lastSD_File,
           initSD, header[3], num;
  double   initAmpl, *finalPos;
  ULL      topBit;
  STORE_F  *memPtr, *localPtr;
  FILE     *filePtr;

  sprintf(basFile,"%s%s",sdStore->title, NONDIAG);

  sprintf(sectFile,"%s%d",basFile,(*sectFileNo)++);
  initSD = (sdStore->lastSD_mem_asym+1);
  lastSectFileNo = sdStore->sect_file_no;
  lastSD_File    = sdStore->lastSD_file_asym;

  if((filePtr = fopen(sectFile,"rb")) == NULL) {
    printf("\n\nError in function idNondiagAsymSD_FileContr():");
    printf("\nNot allowed to open the file %s\n", sectFile);
    exit(1);
  }
  rewind(filePtr);

  topBit = ULL_ONE << (8 * sizeof(ULL) - 1);

  for( ; ; ) {  // loop through blocks all precalculated <SD'()|OP|asymSD()>
  
   // read header[] elements from current data block

    if(fread((void *)&header,(size_t) sizeof(int), 3, filePtr) != 3)  {
      printf("\n\nError in function idNondiagAsymSD_FileContr():");
      printf("\nIn reading nondiagonal header elements from file %s\n",
                                                             sectFile);
      exit(1);
    }
    if(header[2] == -1) {     
      fclose(filePtr);
      if((initSD - 1) != lastSD_File) {
	printf("\n\nError in function idNondiagAsymSD_FileContr():");
	printf("\n %s: initSD = %d  lastSD_file= %d\n", 
   	               sectFile, initSD, lastSD_File);
	exit(1);
      }
      break; 
    } // no more asym matr elem
    
    if(header[2] == -2) {  // change to new sectFile
      fclose(filePtr);
      if((*sectFileNo) > lastSectFileNo) break;  // no more section files
      sprintf(sectFile,"%s%d",basFile, (*sectFileNo)++);
      if((filePtr = fopen(sectFile,"rb")) == NULL) {
	printf("\n\nError in function idNondiagAsymSD_FileContr():");
	printf("\nNot allowed to open the file %s\n", sectFile);
	exit(1);
      }
      continue;  // read header file from new sect_file
    }
    if(initSD != header[0]) {
      printf("\n\nError in function idNondiagAsymSD_FileContr():");
      printf("\n %s: initSD = %d  header[0] = %d\n", 
   	               sectFile, initSD, header[0]);
      exit(1);
    }
    // read list of nondiag matrix elements in format float

    memPtr = (STORE_F *)(sdStore->windMemPtr);

    num = fread((void *)memPtr,(size_t) sizeof(STORE_F), header[2], filePtr);

    if(num != header[2])  {
      printf("\n\nError in function idNondiagAsymSD_FileContr();");
      printf("\n num elements read tall = %d",num);
      printf("\nIn reading %d nondiag asym matrix elements from file %s\n",
                                                      header[2], sectFile);
      exit(1);
    }
    for( ; initSD <= header[1]; initSD++) {
      initAmpl = initVec[initSD];
      localPtr = memPtr;
      while(!(localPtr->final & topBit)) {
	switch((localPtr->final) & 3)   {
	  case 0:  
	  case 3: finalVec[localPtr->final >> 2] 
                     += ((double)localPtr->value) * initAmpl;
                  break;
          case 1: finalVec[localPtr->final >> 2] 
                     += 2 * ((double)localPtr->value) * initAmpl;
                  break;
          case 2: break;
	}  // end switch()
	localPtr++;
      } // end while() loop

      finalPos = &finalVec[initSD];
      while(!(memPtr->final & topBit)) {
	*finalPos += ((double)memPtr->value) * initVec[memPtr->final >> 2];
	memPtr++;
      } // end while()
      memPtr++;
    }   // end contribution from all precalculated |asymSD()> in a block
  }  // end loop through all precalculated nondiag asym matrix elements

} // End: funtion idNondiagAsymSD_FileContr()

     /*
     ** The function 
     **         idNondiagSymSD_FileContr()
     ** calculates contribution from all precalculated 
     ** non-diagonal matrix elements 
     **                 <SD'()|H|symSD()>
     ** to the new lanczo vector |finalVec>.
     */

static void idNondiagSymSD_FileContr(SD_BAS *sdBas,
                                int *sectFileNo, SD_STORE *sdStore,
                                double *initVec, double *finalVec)
{
  char     basFile[ONE_LINE], sectFile[ONE_LINE];
  int      lastSectFileNo , lastSD_File, 
           initSD, header[3], num;
  double   initAmpl, *finalPos;
  ULL      topBit;
  STORE_F  *memPtr, *localPtr;
  FILE     *filePtr;

  sprintf(basFile,"%s%s",sdStore->title,NONDIAG);

  sprintf(sectFile,"%s%d",basFile,(*sectFileNo)++);
  initSD = MAX(sdBas->numSD[0],(sdStore->lastSD_mem_sym+1));
  lastSectFileNo = sdStore->sect_file_no;
  lastSD_File    = sdStore->lastSD_file_sym;

  if((filePtr = fopen(sectFile,"rb")) == NULL) {
    printf("\n\nError in function  idNondiagSymSD_FileContr:");
    printf("\nNot allowed to open the file %s\n", sectFile);
    exit(1);
  }
  rewind(filePtr);

  topBit = ULL_ONE << (8 * sizeof(ULL) - 1);

  for( ; ; ) {  // loop through blocks all precalculated <SD'()|OP|asymSD()>
  
   // read header[] elements from current data block

    if(fread((void *)&header,(size_t) sizeof(int), 3, filePtr) != 3)  {
      printf("\n\nError in function  idNondiagSymSD_FileContr():");
      printf("\nIn reading nondiagonal header elements from file %s\n",
                                                             sectFile);
      exit(1);
    }
    if(header[2] == -1) {     
      fclose(filePtr);
      if((initSD - 1) != lastSD_File) {
	printf("\n\nError in function idNondiagSymSD_FileContr():");
	printf("\n %s: initSD = %d  lastSD_file_asym = %d\n", 
   	               sectFile, initSD, lastSD_File);
	exit(1);
      }
      break; 
    }   // no more asym matr elem
    
    if(header[2] == -2) {  // change to new sect_file
      fclose(filePtr);
      if((*sectFileNo) > lastSectFileNo) break;  // no more section files
      sprintf(sectFile,"%s%d",basFile, (*sectFileNo)++);
      if((filePtr = fopen(sectFile,"rb")) == NULL) {
	printf("\n\nError in function idNondiagSymSD_FileContr():");
	printf("\nNot allowed to open the file %s\n", sectFile);
	exit(1);
      }
      continue;  // read header file from new sect_file
    }
    if(initSD != header[0]) {
      printf("\n\nError in function idNondiagSymSD_FileContr():");
      printf("\n %s: initSD = %d  header[0] = %d\n", 
   	               sectFile, initSD, header[0]);
      exit(1);
    }
    // read list of nondiag matrix elements in format float

    memPtr = (STORE_F *)(sdStore->windMemPtr);

    num = fread((void *)memPtr,(size_t) sizeof(STORE_F), header[2], filePtr);

    if(num != header[2])  {
      printf("\n\nError in function idNondiagSymSD_FileContr();");
      printf("\n num elements read tall = %d",num);
      printf("\nIn reading %d nondiag asym matrix elements from file %s\n",
                                                      header[2], sectFile);
      exit(1);
    }

    for( ; initSD <= header[1]; initSD++) {
      initAmpl = initVec[initSD];
      localPtr = memPtr;
      while(!(localPtr->final & topBit)) {
	if(((localPtr->final) & 3) == 0) {
	  finalVec[(localPtr->final >> 2)] 
             += ((double)localPtr->value)*initAmpl;
	}
	localPtr++;
      } // end while()
      finalPos = &finalVec[initSD];  
      while(!(memPtr->final & topBit)) {
	switch((memPtr->final) & 3)   {
	    case 0: *finalPos 
                     += ((double)memPtr->value)*initVec[memPtr->final >> 2];
                    break;
	    case 1: *finalPos 
                     += ((double)memPtr->value)*initVec[memPtr->final >> 2];
                    break;
	} // end switch() loop
	memPtr++;
      } /* end while() loop */
      memPtr++;
    } // end contribution from current block of calculated |symSD()>

  }  // end loop through all precalculated nondiag sym matrix elements

} // End: funtion idNondiagSymSD_FileContr()

     /*
     ** The function 
     **         idNondiagNosymSD_FileContr()
     ** calculates contribution from all precalculated 
     ** non-diagonal matrix elements 
     **                 <SD'()|H|nosymSD()>
     ** to the new lanczo vector |finalVec>.
     */

static void idNondiagNosymSD_FileContr(SD_BAS *sdBas,
                                int *sectFileNo, SD_STORE *sdStore,
                                double *initVec, double *finalVec)
{
  char     basFile[ONE_LINE], sectFile[ONE_LINE];
  int      lastSectFileNo , lastSD_File, 
           initSD, header[3], num;
  double   initAmpl, *finalPos;
  ULL      topBit;
  STORE_F  *memPtr, *localPtr;
  FILE     *filePtr;

  sprintf(basFile,"%s%s",sdStore->title, NONDIAG);
  sprintf(sectFile,"%s%d",basFile,(*sectFileNo)++);
  initSD = (sdStore->lastSD_mem_asym+1);
  lastSectFileNo = sdStore->sect_file_no;
  lastSD_File    = sdStore->lastSD_file_asym;

  if((filePtr = fopen(sectFile,"rb")) == NULL) {
    printf("\n\nError in function idNondiagNosymSD_FileContr:");
    printf("\nNot allowed to open the file %s\n", sectFile);
    exit(1);
  }
  rewind(filePtr);

  topBit = ULL_ONE << (8 * sizeof(ULL) - 1);

  for( ; ; ) {  // loop through blocks all precalculated <SD'()|OP|asymSD()>
  
   // read header[] elements from current data block

    if(fread((void *)&header,(size_t) sizeof(int), 3, filePtr) != 3)  {
      printf("\n\nError in function idNondiagNosymSD_FileContr:");
      printf("\nIn reading nondiagonal header elements from file %s\n",
                                                             sectFile);
      exit(1);
    }
    if(header[2] == -1) {     
      fclose(filePtr);
      if((initSD - 1) != lastSD_File) {
	printf("\n\nError in function idNondiagNosymSD_FileContr():");
	printf("\n %s: initSD = %d  lastSD_file_asym = %d\n", 
   	               sectFile, initSD, lastSD_File);
	exit(1);
      }
      break; 
    }   // no more asym matr elem
    
    if(header[2] == -2) {  // change to new sect_file
      fclose(filePtr);
      if((*sectFileNo) > lastSectFileNo) break;  // no more section files
      sprintf(sectFile,"%s%d",basFile, (*sectFileNo)++);
      if((filePtr = fopen(sectFile,"rb")) == NULL) {
	printf("\n\nError in function idNondiagNosymSD_FileContr():");
	printf("\nNot allowed to open the file %s\n", sectFile);
	exit(1);
      }
      continue;  // read header file from new sect_file
    }
    if(initSD != header[0]) {
      printf("\n\nError in function idNondiagNosymSD_FileContr():");
      printf("\n %s: initSD = %d  header[0] = %d\n", 
   	               sectFile, initSD, header[0]);
      exit(1);
    }
    // read list of nondiag matrix elements in format float

    memPtr = (STORE_F *)(sdStore->windMemPtr);

    num = fread((void *)memPtr,(size_t) sizeof(STORE_F), header[2], filePtr);

    if(num != header[2])  {
      printf("\n\nError in function idNondiagNosymSD_FileContr()");
      printf("\n num elements read tall = %d",num);
      printf("\nIn reading %d nondiag nosym matrix elements from file %s\n",
                                                      header[2], sectFile);
      exit(1);
    }
    for( ; initSD <= header[1]; initSD++) {
      initAmpl = initVec[initSD]; 
      localPtr = memPtr;
      while(!(localPtr->final & topBit)) {
	finalVec[localPtr->final] += ((double)(localPtr->value)*initAmpl);
	localPtr++;
      } // end while() loop
      finalPos = &finalVec[initSD]; 
      while(!(memPtr->final & topBit)) {
	*finalPos += ((double)memPtr->value)*initVec[memPtr->final];
	memPtr++;
      } // end while() loop
      memPtr++;
    } // end contribution from all precalculated |nosymSD()> in a block
  }  // end loop through all precalculated nondiag nosym matrix elements

} // End: funtion idNondiagNosymSD_FileContr()

     /*
     ** The function 
     **         id_nondiag_asymSD_calc_contr()
     ** calculates the contribution from the non-precalculated 
     ** non-diagonal matrix elements 
     **                 <SD'()|H|asymSD()>
     ** to the new lanczo vector |finalVec> 
     */

static void id_nondiag_asymSD_calc_contr(int initSD, SP_BAS *spBas,
                                           SD_BAS *sdBas, SD_STORE *sdStore,
                                            double *initVec, double *finalVec)
{
  ULL      topBit = ULL_ZERO;
  HEAD     head;
  double   constAmpl, *finalPos; 
  STORE_F  *memPtr, *localPtr;

  topBit      = ULL_ONE << (8 * sizeof(ULL) - 1); //initialization
  head.typeSD  = 0; // asym case
  head.memSize = sdStore->windMemSize;
  head.memPtr  = sdStore->windMemPtr;
  head.endSD   = sdBas->numSD[0];

  while(initSD < sdBas->numSD[0]) { // more contribution
    head.startSD = initSD;
    head.numSD   = 0;
    if(sdStore->typeInt == VEFF_INT3) {
      threePartNondiagSD_Calc(&head, spBas, sdBas, &sdStore->op);
    }
    else {
      twoPartNondiagSD_Calc(&head, spBas, sdBas, &sdStore->op);
    }
    memPtr = (STORE_F *)head.memPtr;

    for( ; initSD <= head.endSD; initSD++) {
      constAmpl = initVec[initSD];
      localPtr = memPtr;
      while(!(localPtr->final & topBit)) {
	switch((localPtr->final) & 3)   {
	    case 0:  
	    case 3: finalVec[localPtr->final >> 2] 
                       += ((double)localPtr->value) * constAmpl;
                    break;
	    case 1: finalVec[localPtr->final >> 2] 
                       += 2 * ((double)localPtr->value) * constAmpl;
                     break;
	    case 2:  break;
	}  // end switch()
	localPtr++;
      } // end while() loop

      finalPos = &finalVec[initSD];
      while(!(memPtr->final & topBit)) {
	*finalPos += ((double)memPtr->value) * initVec[memPtr->final >> 2];
	memPtr++;
      } // end while()
      memPtr++;
    }   // end contribution from all <SD'|VEFF|SD> in current block
  }  // end loop through all asym matrix elements

} // End: funtion id_nondiag_asymSD_calc_contr()

     /*
     ** The function 
     **         id_nondiag_symSD_calc_contr()
     ** calculates the contribution from the non-precalculated 
     ** non-diagonal matrix elements 
     **                 <SD'()|H|asymSD()>
     ** to the new lanczo vector |finalVec> 
     */

static void id_nondiag_symSD_calc_contr(int initSD, SP_BAS *spBas,
                                            SD_BAS *sdBas, SD_STORE *sdStore,
                                             double *initVec, double *finalVec)
{
  ULL      topBit = ULL_ZERO;
  HEAD     head;
  double   constAmpl, *finalPos; 
  STORE_F  *memPtr, *localPtr;

  topBit   = ULL_ONE << (8 * sizeof(ULL) - 1); // initialization
  head.typeSD  = 1; // sym case
  head.memSize = sdStore->windMemSize;
  head.memPtr  = sdStore->windMemPtr;
  head.endSD   = sdBas->numSD[0];

  while(initSD < sdBas->tot_dimSD) {  // more contribution to calculate
    head.startSD = initSD;
    head.numSD   = 0;
    if(sdStore->typeInt == VEFF_INT3) {
      threePartNondiagSD_Calc(&head, spBas, sdBas, &sdStore->op);
    }
    else {
      twoPartNondiagSD_Calc(&head, spBas, sdBas, &sdStore->op);
    }
    memPtr = (STORE_F *)head.memPtr;
    for( ; initSD <= head.endSD; initSD++) {
      constAmpl = initVec[initSD];
      localPtr = memPtr;
      while(!(localPtr->final & topBit)) {
	if(((localPtr->final) & 3) == 0) {
	  finalVec[(localPtr->final >> 2)] 
             += ((double)localPtr->value)*constAmpl;
	}
	localPtr++;
      } // end while()
      finalPos = &finalVec[initSD];  
      while(!(memPtr->final & topBit)) {
	switch((memPtr->final) & 3)   {
	    case 0: *finalPos += ((double)memPtr->value)
                                *initVec[memPtr->final >> 2];
                     break;
	    case 1: *finalPos +=  ((double)memPtr->value)
                                 *initVec[memPtr->final >> 2];
                     break;
	} // end switch() loop
	memPtr++;
      } /* end while() loop */
      memPtr++;
    }   // end contribution from all <SD'|VEFF|SD> in current block
  }  // end loop through all sym matrix elements
} // End: function id_nondiag_symSD_contr()

     /*
     ** The function 
     **         id_nondiag_nosymSD_calc_contr()
     ** calculates the contribution from the non-precalculated 
     ** non-diagonal matrix elements 
     **                 <SD'()|H|nosymSD()>
     ** to the new lanczo vector |finalVec> 
     */

static void id_nondiag_nosymSD_calc_contr(int initSD, SP_BAS *spBas,
                                               SD_BAS *sdBas, SD_STORE *sdStore, 
                                               double *initVec, double *finalVec)
{
  ULL      topBit = ULL_ZERO;
  HEAD     head;
  double   constAmpl, *finalPos; 
  STORE_F  *memPtr, *localPtr;

  topBit   = ULL_ONE << (8 * sizeof(ULL) - 1);  /* initialization */
  head.typeSD  = 2; // nosym case
  head.memSize = sdStore->windMemSize;
  head.memPtr  = sdStore->windMemPtr;
  head.endSD   = sdBas->numSD[0];

  while(initSD < sdBas->numSD[0]) { // more contribution

    head.startSD = initSD;
    head.numSD   = 0;
    if(sdStore->typeInt == VEFF_INT3) {
      threePartNondiagSD_Calc(&head, spBas, sdBas, &sdStore->op);
    }
    else {
      twoPartNondiagSD_Calc(&head, spBas, sdBas, &sdStore->op);
    }
    memPtr = (STORE_F *)head.memPtr;

    for( ; initSD <= head.endSD; initSD++) {
      constAmpl = initVec[initSD]; 
      localPtr = memPtr;
      while(!(localPtr->final & topBit)) {
	finalVec[localPtr->final] += ((double)(localPtr->value)*constAmpl);
	localPtr++;
      } // end while() loop
      finalPos = &finalVec[initSD]; 
      while(!(memPtr->final & topBit)) {
	*finalPos += ((double)memPtr->value)*initVec[memPtr->final];
	memPtr++;
      } // end while() loop
      memPtr++;
    }   // end contribution from all <SD'|VEFF|SD> in current block
  }  // end loop through all nosym matrix elements
}  // End: function id_nondiag_nosymSD_contr()
