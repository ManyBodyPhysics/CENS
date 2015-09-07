  /****************  The module id-lanczo-matrix-element.c  ******************/

#include "shell-model.h"

   /*
   ** The entrance function
   **         id_lanczo_matrix_element()
   ** calculates and return 
   **     <lanc_vec1|OP1lanc_vec2>
   ** where the operator OP = ANG or CM  
   */

               /****  local function declarations ****/

static double idDiagMemSD_Contr(SD_STORE  *sdStore, SD_BAS *sdBas,
                              double *initVec, double *finalVec);
   /*
   ** calculates contribution to 
   **       result = <initVec|OP|finalVec>
   ** from the precalculated diagonal matrix elements
   **  stored as float values in memory
   */

static double nondiagMemContr(SD_STORE *sdStore, SD_BAS *sdBas, 
			      double *initVec, double *finalVec);
   /*
   ** calculates contribution to
   **      result = <initVec|OP|finalVec>
   ** from precalculated nondiag <SD'|OP|SD> stored
   ** in memory.
   */

static double idDiagFileSD_Contr(SD_STORE *sdStore,SD_BAS *sdBas,
                               double *initVec, double *finalVec);
     /*
     ** The function 
     **        idDiagFileSD_Contr()
     ** calculates diagonal contributions to 
     **     result =   <finalVec| H |initVec>
     ** from the precalculated diagonal matrix elements
     **  stored as float values on file 
     */

static double nondiagFileContr(SD_STORE *sdStore,SD_BAS *sdBas, 
                           double *initVec, double *finalVec);
   /*
   ** calculates nondiagonal contributions to 
   **     result = <finalVec_SD'|OP|initVec_SD>
   ** from the precalculated nondiagonal matrix elements
   ** stored as float values on file 
   */
 
static double  nondiagExplicitCalc(SD_STORE *sdStore, 
                                SP_BAS *spBas, SD_BAS *sdBas, 
				double *initVec, double *finalVec);
   /*
   ** calculates explicitly nondiagonal contributions to 
   **      result =  <finalVec|OP|initVec>
   ** not stored on file 
   */

static double idNondiagSymMemContr(SD_BAS *sdBas, SD_STORE *sdStore, 
				   double *initVec, double *finalVec);
   /*
   ** calculates contribution from all precalculated 
   ** non-diagonal matrix elements of  
   **  <SD'()|OP|asymSD()>  and  <SD'()|OP|symSD()
   ** stored in memory to the matrix elements
   **    result = <finalVec|OP|initVec>
   */

static double idNondiagNosymMemContr(SD_BAS *sdBas, SD_STORE *sdStore,
				     double *initVec, double *finalVec);
   /*
   ** calculates contribution from all precalculated 
   ** non-diagonal matrix elements of  
   **         <SD'()|OP|nosymSD()>
   ** stored in memory to the matrix elements
   **    result = <finalVec|OP|initVec>
   */

static double idNondiagAsymSD_FileContr(SD_BAS *sdBas,
                          int *sectFileNo,SD_STORE *sdStore,
                           double *initVec, double *finalVec);
     /*
     ** calculates contribution from all precalculated 
     ** non-diagonal matrix elements to 
     **       result = <finalVec|OP|initVec(asymSD)>
     ** stored on file to the new lanczo vector |finalVec>.
     */

static double idNondiagSymSD_FileContr(SD_BAS *sdBas,
                                int *sectFileNo, SD_STORE *sdStore,
				  double *initVec, double *finalVec);
     /*
     ** calculates contribution from all precalculated 
     ** non-diagonal matrix elements to 
     **       result = <finalVec|OP|initVec(symSD)>
     ** to the new lanczo vector |finalVec>.
     */

static double idNondiagNosymSD_FileContr(SD_BAS *sdBas,
                                int *sectFileNo, SD_STORE *sdStore,
				       double *initVec, double *finalVec);
     /*
     ** calculates contribution from the 
     ** non-diagonal matrix elements to 
     **       result = <finalVec|OP|initVec(nosymSD)>
     ** to the new lanczo vector |finalVec>.
     */

static double idNondiagAsymSD_CalcContr(int initSD, SP_BAS *spBas,SD_BAS *sdBas,
                           SD_STORE *sdStore,double *initVec, double *finalVec);
     /*
     ** calculates the contribution from the
     ** non-diagonal matrix elements to
     **     result = <finalVec|OP|initVec(asymSD)>
     ** not store on file 
     */

static double idNondiagSymSD_CalcContr(int initSD, SP_BAS *spBas,
                                          SD_BAS *sdBas, SD_STORE *sdStore,
				       double *initVec, double *finalVec);
     /*
     ** calculates the contribution from the
     ** non-diagonal matrix elements 
     **     result = <finalVec|OP|initVec(symSD)>
     ** not store on file 
     */

static double idNondiagNosymSD_CalcContr(int initSD, SP_BAS *spBas,
                                          SD_BAS *sdBas, SD_STORE *sdStore,
					 double *initVec, double *finalVec);
     /*
     ** calculates the contribution from the non-precalculated 
     ** non-diagonal matrix elements 
     **     result = <finalVec|OP|initVec(nosymSD)>
     ** not store on file 
     */
               /****   end local function declarations ****/

               /**** The function definitions  ****/ 
 
   /*
   ** The entrance function
   **         id_lanczo_matrix_element()
   ** calculates and return 
   **     <lanc_vec1|OP1lanc_vec2>
   ** where the operator OP = ANG or CM  
   */

double id_lanczo_matrix_element(SP_BAS *spBas, SD_BAS *sdBas,
              SD_STORE *sdStore, double *initVec, double *finalVec)
{
  double  result;

  result = 0,0;

  if(sdStore->diagMemCalc == YES) {
    result += idDiagMemSD_Contr(sdStore, sdBas, initVec, finalVec);
  }

  if(  (sdStore->lastSD_mem_asym > -1)    
     ||(sdStore->lastSD_mem_sym  > -1)) {

    // <SD'|OP|SD> stored in memory

    result += nondiagMemContr(sdStore, sdBas, initVec, finalVec);

  } // end contribution from memory

  if(sdStore->diagMemCalc == NO) {
    result += idDiagFileSD_Contr(sdStore, sdBas, initVec, finalVec);
  }

  if(   (sdStore->lastSD_file_asym > -1)
     || (sdStore->lastSD_file_sym  > -1)) {
	
    //  <SD'|VEFF|SD> stored on files

    result += nondiagFileContr(sdStore, sdBas, initVec, finalVec);

  }  // end contribution from file

  if(  (   (sdBas->numSD[0] > 0) 
	&& (sdStore->lastSD_file_asym < (sdBas->numSD[0] - 1)))
     ||(   (sdBas->numSD[1] > 0) 
        && (sdStore->lastSD_file_sym  < (sdBas->tot_dimSD - 1)))) {

    //  explicit calculation of <SD'|OP|SD> 

    result += nondiagExplicitCalc(sdStore, spBas, sdBas,
                                       initVec, finalVec);

  } // end contribution explicit calculation

  return  result;

} // End: function   id_lanczo_matrix_element()

     /* The function
     **        idDiagMemSD_Contr()
     ** calculates contribution to 
     **       result = <initVec|OP|finalVec>
     ** from the precalculated diagonal matrix elements
     **  stored as float values in memory
     */

static double idDiagMemSD_Contr(SD_STORE *sdStore, SD_BAS *sdBas,
                              double *initVec, double *finalVec)
{
  int       k;
  float     *mem_ptr;
  double    result, factor;

  result = 0.0;
  mem_ptr = (float *)sdStore->permMemPtr;

  factor = (sdBas->MJ == 0) ? 2.0 : 1.0;

  for(k = 0; k < sdBas->numSD[0]; k++) {
    result +=  factor*((double)(*mem_ptr++)) 
              * (*(initVec++)) * (*(finalVec++));
  }
  for(k = 0; k < sdBas->numSD[1]; k++) {
    result +=  (double)(*mem_ptr++) 
             * (*(initVec++)) * (*(finalVec++));
  }

  return result;

} // End: function idDiagMemSD_Contr() 

     /*
     ** The function 
     **      nondiagMemContr()
     ** calculates contribution to
     **   result = <initVec|OP|finalVec>
     ** from precalculated nondiag <SD'|OP|SD> stored
     ** in memory.
     */

static double nondiagMemContr(SD_STORE *sdStore, SD_BAS *sdBas, 
		    double *initVec, double *finalVec)
{
  double  result;

  result = 0.0;

  if(sdBas->MJ == 0) {
    result += idNondiagSymMemContr(sdBas, sdStore,
                               initVec, finalVec);
  }
  else {
    result += idNondiagNosymMemContr(sdBas, sdStore,
                                 initVec, finalVec);
  }
  return result;

} // End function nondiagMemContr()

     /*
     ** The function 
     **        idDiagFileSD_Contr()
     ** calculates diagonal contributions to 
     **       result = <finalVec| H |initVec>
     ** from the precalculated diagonal matrix elements
     **  stored as float values on file 
     */

static double idDiagFileSD_Contr(SD_STORE *sdStore, SD_BAS *sdBas,
                                double *initVec,  double *finalVec)
{
 char    sect_file[ONE_LINE];
  int     k, header[2], num_elem;
  float   *mem_ptr;
  double  result, factor;
  FILE    *file_ptr;

  result = 0.0; // initialization

  sprintf(sect_file,"%s%s",sdStore->title,DIAG);

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

  factor = (sdBas->MJ == 0) ? 2.0 : 1.0;

  for(k = 0; k < sdBas->numSD[0]; k++,initVec++, finalVec++) {
    result +=  factor*((double)(*mem_ptr++)) 
              * (*initVec) * (*finalVec);
  }
  for(k = 0; k < sdBas->numSD[1]; k++,initVec++, finalVec++) {
    result +=  (double)(*mem_ptr++) 
             * (*initVec) * (*finalVec);
  }

  fclose(file_ptr);  // close sect_file

  return result;

} // End: function idDiagFileSD_Contr()

     /*
     ** The function 
     **      nondiagFileContri()
     ** calculates nondiagonal contributions to 
     **      result = <finalVec_SD'|OP|initVec_SD>
     ** from the precalculated nondiagonal matrix elements
     ** stored as float values on file 
     */

static double nondiagFileContr(SD_STORE *sdStore,SD_BAS *sdBas, 
                              double *initVec, double *finalVec)
{
  int      lastSD_asym, lastSD_sym,  sect_file_no;
  double   result;
  
  result = 0.0; // initialization

  lastSD_asym =  sdStore->lastSD_file_asym;
  lastSD_sym  =  sdStore->lastSD_file_sym;

  sect_file_no = 1; // initialize first nondiag section file number

  if(   (sdBas->MJ == 0) && (lastSD_asym > -1)) {
    result = idNondiagAsymSD_FileContr(sdBas,&sect_file_no,sdStore,
                                                initVec, finalVec);
  } // end asym contribution

  if(   (sdBas->MJ == 0) && (lastSD_sym > -1)) {
    result = idNondiagSymSD_FileContr(sdBas, &sect_file_no,sdStore,
                                                initVec, finalVec);
  } // end sym contribution

  if(   (sdBas->MJ != 0) && (lastSD_asym > -1)) {
    result += idNondiagNosymSD_FileContr(sdBas, &sect_file_no, sdStore,
				                   initVec, finalVec);
  } // end nosym contribution

  return result;

} // End: function  nondiagFileContr()

    /*
    ** The function 
    **      nondiagExplicitCal()
    ** calculates explicitly nondiagonal contributions to 
    **     result = <finalVec|OP|initVec>
    ** not stored on file 
    */

static double nondiagExplicitCalc(SD_STORE *sdStore, 
                                SP_BAS *spBas, SD_BAS *sdBas, 
                                 double *initVec, double *finalVec)
{
  int     initSD, lastSD_asym, lastSD_sym;
  double  result;

  lastSD_asym = MAX(sdStore->lastSD_mem_asym, 
                    sdStore->lastSD_file_asym);
;
  lastSD_sym  = MAX(sdStore->lastSD_mem_sym, 
                    sdStore->lastSD_file_sym);

  result = 0.0; // initialization 

  if(   (sdBas->MJ == 0)
     && (lastSD_asym < (sdBas->numSD[0] - 1))) {

    initSD = lastSD_asym + 1;
    result += idNondiagAsymSD_CalcContr(initSD, spBas, sdBas, 
                                          sdStore, initVec,finalVec);
  } // end time-reversal asym contribution

  if(   (sdBas->MJ == 0)
     && (lastSD_sym < (sdBas->tot_dimSD - 1))) {
    initSD  = lastSD_sym + 1;
    result += idNondiagSymSD_CalcContr(initSD, spBas, sdBas,
                                        sdStore, initVec, finalVec);
  }  // end time-reversal sym contribution 

 if(   (sdBas->MJ != 0)
      &&(lastSD_asym < (sdBas->numSD[0] - 1))) {
    initSD = lastSD_asym + 1;
    result += idNondiagNosymSD_CalcContr(initSD, spBas,sdBas,
                                         sdStore, initVec,finalVec);
  }  // end no time-reversal symmetry contribution
  
 return  result;

} // End: function  nondiagExplicitCalc()

     /*
     ** The function 
     **         idNondiagSymMemContr()
     ** calculates contribution from all precalculated 
     ** non-diagonal matrix elements of  
     **  <SD'()|OP|asymSD()>  and  <SD'()|OP|symSD()
     ** stored in memory to the matrix elements
     **    result = <finalVec|OP|initVec>
     */

static double  idNondiagSymMemContr(SD_BAS *sdBas, SD_STORE *sdStore, 
                                        double *initVec, double *finalVec)
{
  int     initSD, vecPos;
  ULL     topBit;
  double  result, value, constAmpl;
  STORE_F *memPtr, *localPtr;

  topBit = ULL_ONE << (8 * sizeof(ULL) - 1);
  result = 0.0;

  memPtr = (STORE_F *)(sdStore->permMemPtr + sdBas->tot_dimSD);
  for(initSD = 0; initSD <= sdStore->lastSD_mem_asym; initSD++) { 
      constAmpl = initVec[initSD];
      localPtr = memPtr;
      while(!(localPtr->final & topBit)) {
	vecPos = localPtr->final >> 2;
	  switch((localPtr->final) & 3)   {
	  case 0:  
	  case 3: value =  ((double)localPtr->value) 
                         * constAmpl * finalVec[vecPos]; 
	          if(vecPos < sdBas->numSD[0]) {
		    result += 2.0 * value;
	          }
		  else {
		    result += value;
		  }
                  break;
	  case 1: value =   2 * ((double)localPtr->value) 
                          * constAmpl * finalVec[vecPos];
	          if(vecPos < sdBas->numSD[0]) {
		    result += 2.0 * value;
		  }
		  else {
		    result += value;
		  }
		  break;
	  case 2:  break;
	}  // end switch()
	localPtr++;
      }  // end while() loop

      constAmpl = finalVec[initSD];
      while(!(memPtr->final & topBit)) {
	vecPos = memPtr->final >> 2;
	switch((memPtr->final) & 3) {
	  case 0: value = ((double)memPtr->value) 
                         *initVec[vecPos] * constAmpl;
	          if(vecPos < sdBas->numSD[0]) {
		    result += 2.0 * value;
		  }
		  else {
		    result += value;
		  }
		  break;	
	default: result +=  2.0*((double)memPtr->value) 
                              * initVec[vecPos] 
                              * constAmpl;
	          break;
	} // end switch()
	memPtr++;
      } // end while()
      memPtr++;
      
  }  // end loop through all precalculated nondiag asym matrix elements

  for(initSD = sdBas->numSD[0]; initSD <= sdStore->lastSD_mem_sym;
                                                          initSD++) { 
     constAmpl = initVec[initSD];
      localPtr = memPtr;
      while(!(localPtr->final & topBit)) {
	vecPos = localPtr->final >> 2; 
	if(((localPtr->final) & 3) == 0) {
	  value =  ((double)localPtr->value) 
                  *constAmpl * finalVec[vecPos];
          if(vecPos < sdBas->numSD[0]) {
	    result += 2.0 *value;
	  }
	  else  {
	    result += value;
	  }
	}
	localPtr++;
      } // end while()
      constAmpl = finalVec[initSD];
      while(!(memPtr->final & topBit)) {
	vecPos = memPtr->final >> 2;
	switch((memPtr->final) & 3)   {
	    case 0: result += ((double)memPtr->value)
                              *initVec[vecPos] * constAmpl;
                    break;
	    case 1: result += ((double)memPtr->value)
                             *initVec[vecPos] * constAmpl;
                    break;
	} // end switch() loop
	memPtr++;
      } /* end while() loop */
      memPtr++;
  }  // end loop through all precalculated nondiag sym matrix elements

  return result;

} // End: funtion idNondiagSymMemContr()

     /*
     ** The function 
     **         idNondiagNosymMemContr()
     ** calculates contribution from all precalculated 
     ** non-diagonal matrix elements of  
     **         <SD'()|OP|nosymSD()>
     ** stored in memory to the matrix elements
     **    result = <finalVec|OP|initVec>
     */

static double idNondiagNosymMemContr(SD_BAS *sdBas, SD_STORE *sdStore,
                                   double *initVec, double *finalVec)
{
  int       initSD;
  ULL       topBit;
  double    result, constAmpl;
  STORE_F   *memPtr, *localPtr;

  topBit = ULL_ONE << (8 * sizeof(ULL) - 1);
  memPtr = (STORE_F *)(sdStore->permMemPtr + sdBas->tot_dimSD); 

  result = 0.0;

  for(initSD = 0; initSD <= sdStore->lastSD_mem_asym; initSD++) {
    constAmpl = initVec[initSD]; 
    localPtr = memPtr;
    while(!(localPtr->final & topBit)) {
      result +=  ((double)(localPtr->value)
                  * constAmpl) * finalVec[localPtr->final];
      localPtr++;
    } // end while() loop
    constAmpl = finalVec[initSD]; 
    while(!(memPtr->final & topBit)) {
      result +=   ((double)memPtr->value)
                  * constAmpl * initVec[memPtr->final];
      memPtr++;
    } // end while() loop
    memPtr++;
  } // end contribution from all precalculated |nosymSD()> in a block
    // end loop through all precalculated nondiag nosym matrix elements

  return result;

}// End: function idNondiagNosymMemContr()

     /*
     ** The function 
     **         idNondiagAsymSD_FileContr()
     ** calculates contribution from all precalculated 
     ** non-diagonal matrix elements to 
     **       result = <finalVec|OP|initVec(asymSD)>
     ** stored on file to the new lanczo vector |finalVec>.
     */

static double idNondiagAsymSD_FileContr(SD_BAS *sdBas,int *sectFileNo,
                  SD_STORE *sdStore,double *initVec, double *finalVec)
{
  char     basFile[ONE_LINE], sectFile[ONE_LINE];
  int      lastSectFileNo, lastSD_File,
           initSD, vecPos, header[3], num;
  double   result, value, constAmpl;
  ULL      topBit;
  STORE_F  *memPtr, *localPtr;
  FILE     *filePtr;

  sprintf(basFile,"%s%s",sdStore->title, NONDIAG);

  sprintf(sectFile,"%s%d",basFile,(*sectFileNo)++);
  initSD         = (sdStore->lastSD_mem_asym+1);
  lastSectFileNo = sdStore->sect_file_no;
  lastSD_File    = sdStore->lastSD_file_asym;

  if((filePtr = fopen(sectFile,"rb")) == NULL) {
    printf("\n\nError in function idNondiagAsymSD_FileContr():");
    printf("\nNot allowed to open the file %s\n", sectFile);
    exit(1);
  }
  rewind(filePtr);

  topBit = ULL_ONE << (8 * sizeof(ULL) - 1);
  result = 0.0;

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


      constAmpl = initVec[initSD];
      localPtr = memPtr;
      while(!(localPtr->final & topBit)) {
	vecPos = localPtr->final >> 2;
	  switch((localPtr->final) & 3)   {
	  case 0:  
	  case 3: value =  ((double)localPtr->value) 
                         * constAmpl * finalVec[vecPos]; 
	          if(vecPos < sdBas->numSD[0]) {
		    result += 2.0 * value;
	          }
		  else {
		    result += value;
		  }
                  break;
	  case 1: value =   2 * ((double)localPtr->value) 
                          * constAmpl * finalVec[vecPos];
	          if(vecPos < sdBas->numSD[0]) {
		    result += 2.0 * value;
		  }
		  else {
		    result += value;
		  }
		  break;
	  case 2:  break;
	}  // end switch()
	localPtr++;
      }  // end while() loop

      constAmpl = finalVec[initSD];
      while(!(memPtr->final & topBit)) {
	vecPos = memPtr->final >> 2;
	switch((memPtr->final) & 3) {
	  case 0: value = ((double)memPtr->value) 
                         *initVec[vecPos] * constAmpl;
	          if(vecPos < sdBas->numSD[0]) {
		    result += 2.0 * value;
		  }
		  else {
		    result += value;
		  }
		  break;	
	default: result +=  2.0*((double)memPtr->value) 
                              * initVec[vecPos] 
                              * constAmpl;
	          break;
	} // end switch()
	memPtr++;
      } // end while()
      memPtr++;
      
    } // end contribution from all precalculated |asymSD()> in a block
  }  // end loop through all precalculated nondiag asym matrix elements

  return result;

} // End: funtion idNondiagAsymSD_FileContr()

     /*
     ** The function 
     **         idNondiagSymSD_FileContr()
     ** calculates contribution from all precalculated 
     ** non-diagonal matrix elements to 
     **       result = <finalVec|OP|initVec(symSD)>
     ** to the new lanczo vector |finalVec>.
     */

static double idNondiagSymSD_FileContr(SD_BAS *sdBas,
                                int *sectFileNo, SD_STORE *sdStore,
                                double *initVec, double *finalVec)
{
  char     basFile[ONE_LINE], sectFile[ONE_LINE];
  int      lastSectFileNo , lastSD_File, 
           initSD, vecPos, header[3], num;
  double   result, value, constAmpl;
  ULL      topBit;
  STORE_F  *memPtr, *localPtr;
  FILE     *filePtr;

  sprintf(basFile,"%s%s",sdStore->title, NONDIAG);

  sprintf(sectFile,"%s%d",basFile,(*sectFileNo)++);
  initSD         = (sdStore->lastSD_mem_asym+1);
  lastSectFileNo = sdStore->sect_file_no;
  lastSD_File    = sdStore->lastSD_file_asym;

  if((filePtr = fopen(sectFile,"rb")) == NULL) {
    printf("\n\nError in function idNondiagSymSD_FileContr():");
    printf("\nNot allowed to open the file %s\n", sectFile);
    exit(1);
  }
  rewind(filePtr);

  topBit = ULL_ONE << (8 * sizeof(ULL) - 1);
  result = 0.0;
  initSD = sdBas->numSD[0];

  for( ; ; ) {  // loop through blocks all precalculated <SD'()|OP|symSD()>
  
   // read header[] elements from current data block

    if(fread((void *)&header,(size_t) sizeof(int), 3, filePtr) != 3)  {
      printf("\n\nError in function idNondiagSymSD_FileContr():");
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
    }   // no more sym matr elem
    
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
      constAmpl = initVec[initSD];
      localPtr = memPtr;
      while(!(localPtr->final & topBit)) {
	vecPos = localPtr->final >> 2; 
	if(((localPtr->final) & 3) == 0) {
	  value =  ((double)localPtr->value) 
                  *constAmpl * finalVec[vecPos];
          if(vecPos < sdBas->numSD[0]) {
	    result += 2.0 *value;
	  }
	  else  {
	    result += value;
	  }
	}
	localPtr++;
      } // end while()
      constAmpl = finalVec[initSD];
      while(!(memPtr->final & topBit)) {
	vecPos = memPtr->final >> 2;
	switch((memPtr->final) & 3)   {
	    case 0: result += ((double)memPtr->value)
                              *initVec[vecPos] * constAmpl;
                    break;
	    case 1: result += ((double)memPtr->value)
                             *initVec[vecPos] * constAmpl;
                    break;
	} // end switch() loop
	memPtr++;
      } /* end while() loop */
      memPtr++;
    } // end contribution from current block of calculated |symSD()>
  }  // end contribution from all block of calculated |symSD()>

  return result;

} // End: funtion idNondiagSymSD_FileContr()

     /*
     ** The function 
     **         idNondiagNosymSD_FileContr()
     ** calculates contribution from all precalculated 
     ** non-diagonal matrix elements to 
     **       result = <finalVec|OP|initVec(nosymSD)>
     ** to the new lanczo vector |finalVec>.
     */

static double idNondiagNosymSD_FileContr(SD_BAS *sdBas,
                                int *sectFileNo, SD_STORE *sdStore,
                                double *initVec, double *finalVec)
{
  char    sectFile[ONE_LINE], basFile[ONE_LINE];
  int     lastSectFileNo, lastSD_File, 
          initSD, header[3], num;
  double  result, constAmpl;
  ULL     topBit;
  STORE_F *memPtr, *localPtr;
  FILE    *filePtr;

  sprintf(basFile,"%s%s",sdStore->title, NONDIAG);

  sprintf(sectFile,"%s%d",basFile,(*sectFileNo)++);
  initSD         = (sdStore->lastSD_mem_asym+1);
  lastSectFileNo = sdStore->sect_file_no;
  lastSD_File    = sdStore->lastSD_file_asym;

  if((filePtr = fopen(sectFile,"rb")) == NULL) {
    printf("\n\nError in function idNondiagNosymSD_FileContr():");
    printf("\nNot allowed to open the file %s\n", sectFile);
    exit(1);
  }
  rewind(filePtr);

  topBit = ULL_ONE << (8 * sizeof(ULL) - 1);
  result = 0.0;
  initSD = 0;

  for( ; ; ) {  // loop through blocks all precalculated <SD'()|OP|asymSD()>
  
   // read header[] elements from current data block

    if(fread((void *)&header,(size_t) sizeof(int), 3, filePtr) != 3)  {
      printf("\n\nError in function idNondiagNosymSD_FileContr():");
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
      printf("\n\nError in function idNondiagNosymSD_FileContr();");
      printf("\n num elements read tall = %d",num);
      printf("\nIn reading %d nondiag nosym matrix elements from file %s\n",
                                                      header[2], sectFile);
      exit(1);
    }
    for( ; initSD <= header[1]; initSD++) {
      constAmpl = initVec[initSD]; 
      localPtr = memPtr;
      while(!(localPtr->final & topBit)) {
	result +=  ((double)(localPtr->value)
                  * constAmpl) * finalVec[localPtr->final];
	localPtr++;
      } // end while() loop
      constAmpl = finalVec[initSD]; 
      while(!(memPtr->final & topBit)) {
	result +=   ((double)memPtr->value)
                  * constAmpl * initVec[memPtr->final];
	memPtr++;
      } // end while() loop
      memPtr++;
    } // end contribution from all precalculated |nosymSD()> in a block
  }  // end loop through all precalculated nondiag asym matrix elements

  return result;

} // End: funtion idNondiagNosymSD_FileContr()

     /*
     ** The function 
     **         idNondiagAsymSD_CalcContr()
     ** calculates the contribution from the
     ** non-diagonal matrix elements to
     **     result = <finalVec|OP|initVec(asymSD)>
     ** not store on file 
     */

static double idNondiagAsymSD_CalcContr(int initSD, SP_BAS *spBas,
                                          SD_BAS *sdBas, SD_STORE *sdStore,
                                          double *initVec, double *finalVec)
{
  ULL      topBit;
  HEAD     head;
  MATR_OP  *op_int;
  double   result, value, constAmpl; 
  STORE_F  *memPtr, *localPtr;

  result       = 0.0;   //initialization
  op_int       = &sdStore->op;
  topBit       = ULL_ONE << (8 * sizeof(ULL) - 1);
  head.typeSD  = 0; // asym case
  head.memSize = sdStore->windMemSize;
  head.memPtr  = sdStore->windMemPtr;
  head.endSD   = sdBas->numSD[0];

  while(initSD < sdBas->numSD[0]) { // more contribution
    head.startSD = initSD;
    twoPartNondiagSD_Calc(&head, spBas, sdBas, op_int);
    memPtr = (STORE_F *)head.memPtr;
    for( ; initSD <=  head.endSD; initSD++) {
      constAmpl = initVec[initSD];
      localPtr = memPtr;
      while(!(localPtr->final & topBit)) {
	switch((localPtr->final) & 3)   {
	  case 0:  
	  case 3: value =  ((double)localPtr->value) 
                         * constAmpl
		         * finalVec[localPtr->final >> 2]; 
	          if(finalVec[localPtr->final >> 2] < sdBas->numSD[0]) {
		    result += 2.0 * value;
	          }
		  else {
		    result += value;
		  }
                  break;
	  case 1: value =   2 * ((double)localPtr->value) 
                          * constAmpl
		          * finalVec[localPtr->final >> 2];
	          if(finalVec[localPtr->final >> 2] < sdBas->numSD[0]) {
		    result += 2.0 * value;
		  }
		  else {
		    result += value;
		  }
		  break;
	  case 2:  break;
	}  // end switch()
	localPtr++;
      } // end while() loop

      constAmpl = finalVec[initSD];
      while(!(memPtr->final & topBit)) {
	switch((localPtr->final) & 3) {
	  case 0: value = ((double)memPtr->value) 
                         *initVec[memPtr->final >> 2]
                         *constAmpl;
	          if(initVec[memPtr->final >> 2] < sdBas->numSD[0]) {
		    result += 2.0 * value;
		  }
		  else {
		    result += value;
		  }
		  break;	
	  default: result +=  2.0*((double)memPtr->value) 
                              * initVec[memPtr->final >> 2] 
                              * constAmpl;
	          break;
	} // end switch()
	memPtr++;
      } // end while()
      memPtr++;
      
    } // end contribution from all precalculated |asymSD()> in a block
  }  // end loop through all precalculated nondiag asym matrix elements

  return result;

} // End: funtion idNondiagAsymSD_CalcContr()

     /*
     ** The function 
     **    idNondiagSymSD_CalcContr()
     ** calculates the contribution from the
     ** non-diagonal matrix elements 
     **                 <SD'()|H|asymSD()>
     *' not stored on file
     */

static double idNondiagSymSD_CalcContr(int initSD, SP_BAS *spBas,
                                          SD_BAS *sdBas, SD_STORE *sdStore,
                                          double *initVec, double *finalVec)
{
  ULL      topBit;
  HEAD     head;
  MATR_OP  *op_int;
  double   result, value, constAmpl; 
  STORE_F  *memPtr, *localPtr;

  result       = 0.0;   //initialization
  op_int       = &sdStore->op;
  topBit       = ULL_ONE << (8 * sizeof(ULL) - 1);
  head.typeSD  = 1; // sym case
  head.memSize = sdStore->windMemSize;
  head.memPtr  = sdStore->windMemPtr;
  head.endSD   = sdBas->tot_dimSD;

  while(initSD < sdBas->numSD[0]) { // more contribution
    head.startSD = initSD;
    twoPartNondiagSD_Calc(&head, spBas, sdBas, op_int);
    memPtr = (STORE_F *)head.memPtr;

    for( ; initSD <= head.endSD; initSD++) {
      constAmpl = initVec[initSD];
      localPtr = memPtr;
      while(!(localPtr->final & topBit)) {
	if(((localPtr->final) & 3) == 0) {
	  value =  ((double)localPtr->value) 
                  *constAmpl
	          *finalVec[(localPtr->final >> 2)];
          if(localPtr->final >> 2 < sdBas->numSD[0]) {
	    result += 2.0 *value;
	  }
	  else  {
	    result += value;
	  }
	}
	localPtr++;
      } // end while()
      constAmpl = finalVec[initSD];
      while(!(memPtr->final & topBit)) {
	switch((memPtr->final) & 3)   {
	    case 0: result += ((double)memPtr->value)
                              *initVec[memPtr->final >> 2]
                              *constAmpl;
                    break;
	    case 1: result += ((double)memPtr->value)
                             *initVec[memPtr->final >> 2]
                             *constAmpl;
                    break;
	} // end switch() loop
	memPtr++;
      } /* end while() loop */
      memPtr++;
    } // end contribution from current block of calculated |symSD()>
  }  // end contribution from all block of calculated |symSD()>

  return result;

} // End: funtion idNondiagSymSD_CalcContr()

     /*
     ** The function 
     **         id_nondiag_nosymSD_calc_contr()
     ** calculates exlicitly contribution from all 
     ** non-diagonal matrix elements to
     **       result = <final_SD'|OP|init_nosymSD>
     ** not stored on file.
     */

static double idNondiagNosymSD_CalcContr(int initSD, SP_BAS *spBas,
                                          SD_BAS *sdBas, SD_STORE *sdStore,
                                          double *initVec, double *finalVec)
{
  ULL      topBit;
  HEAD     head;
  MATR_OP  *op_int;
  double   result, constAmpl; 
  STORE_F  *memPtr, *localPtr;

  result       = 0.0;   //initialization
  op_int       = &sdStore->op;
  topBit       = ULL_ONE << (8 * sizeof(ULL) - 1);
  head.typeSD  = 2; // asym case
  head.memSize = sdStore->windMemSize;
  head.memPtr  = sdStore->windMemPtr;
  head.endSD   = sdBas->numSD[0];

  while(initSD < sdBas->numSD[0]) { // more contribution
    head.startSD = initSD;
    twoPartNondiagSD_Calc(&head, spBas, sdBas, op_int);
    memPtr = (STORE_F *)head.memPtr;
    for( ; initSD <= head.endSD; initSD++) {
      constAmpl = initVec[initSD]; 
      localPtr = memPtr;
      while(!(localPtr->final & topBit)) {
	result +=  ((double)(localPtr->value)
                  * constAmpl) * finalVec[localPtr->final];
	localPtr++;
      } // end while() loop
      constAmpl = finalVec[initSD]; 
      while(!(memPtr->final & topBit)) {
	result +=   ((double)memPtr->value)
                  * constAmpl * initVec[memPtr->final];
	memPtr++;
      } // end while() loop
      memPtr++;
    } // end contribution from all precalculated |nosymSD()> in a block
  }  // end loop through all precalculated nondiag asym matrix elements

  return  result;
 
} // End: funtion idNondiagNosymSD_FileContr()
