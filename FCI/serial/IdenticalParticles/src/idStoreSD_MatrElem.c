/********************  File: idStoreSD_MatrElem.c  *************/

#include "shell-model.h"

  /* 
  **    Two entrance function
  */

     /* 1.
     ** The entrance function
     **       storeSD_MatrElemInMem()
     ** calculates and stores diagonal <SD|OP|SD> and as
     ** many as possible of the nondiagonal <SD()'|OP|SD()>
     ** The interaction is defined through typeInt
     */

     /* 2.
     ** The module entrance function 
     **       storeSD_MatrElemOnFile()
     ** calculates and stores diagonal if(diagCalc == NO)
     ** and as many as possible of the nondiagonal SD
     ** matrix elements
     **              <SD'()|OP|SD()>
     */

	    /**** local function declarations ****/

static void twoPartDiagSD_MatrElem(SD_BAS *sdBas, MATR_OP *op_int, 
				                    float *diag_ptr);
     /*
     ** calculates all diagonal matrix elements 
     **           <SD()|OP()| SD()> 
     ** where the two-particle matrix elements of OP() are stored in
     ** op_int->id_diag[] using the formula
     **        op(k,l) = ((num - k) * k)/2 + (l - 1)
     ** where num = 2 * m_orb - 3. 
     ** The result is stores in diag_SD[]. 
     ** The result is stored in memory pointed to by sdStore->freePermMemPtr.
     */

static void threePartDiagSD_MatrElem(SD_BAS *sdBas, MATR_OP *op_veff, 
				                    float *diag_ptr);
     /*
     ** calculates all diagonal matrix elements 
     **           <SD()|op_veff(three_part)()| SD()> 
     ** where the three-particle matrix elements of OP() are stored in
     ** op_veff->id_diag[] using the formula
     **  op(k,l,m) =  k*(num1 - k*(num2 -k))/6 + (l*(num3 - l))/2
     **             - num_orb + m     
     ** where num1 = 3*num_orb^2 -12*num_orb + 11
     **       num2 = 3*num_orb - 6
     **       num3 = 2*num_orb - 3
     ** The result is stored in memory pointed to by sdStore->freePermMemPtr.
     */

static void writeDiagSD_MatrElem(SD_BAS *sdBas, SD_STORE *sdStore);
     /*
     ** create a new file "file_name" and write diagonal
     ** matrix elements <SD()|OP|SD()>
     ** to file: 
     **     a) Header elements: num_asymSD or num_nosymSD  and num_symSD
     **     b) list of matrix elements 
     */

static void nondiagMemMatrElem(SP_BAS *spBas,SD_BAS *sdBas,SD_STORE *sdStore);
     /*
     **           twoPartNondiagMemMatrElem()
     ** calculates and stores in memory as many as
     ** possible of the nondiagonal SD matrix elements
     **                <SD()'|OP|SD()>
     */

static void  nondiagFileMatrElem(SP_BAS *spBas,SD_BAS *sdBas,SD_STORE *sdStore);
     /*
     ** calculates and stores on file as many as possible
     ** of the nondiagonal SD matrix elements
     **                <SD()'|OP|SD()>
     ** not previously calculated and stored in memory
     ** The interaction Op is defined through typeInt
     */

static int new_section_file(SD_STORE *sdStore);
     /*
     ** creates a new section file with a corresponding
     ** section number and increase section file number
     ** if enough file space is available
     ** The function return TRUE if new section file is 
     ** is created, otherwise FALSE
     */

static void writeNondiagSD_MatrElem(HEAD *head, SD_STORE *sdStore);
     /*
     ** write nondiagonal matrix elements <SD()|OP|SD()>
     ** to file: 
     **     a) Header elements
     **     b) list of matrix elements 
     */

void printMatrElemData(SD_BAS *sdBas, int storeType, SD_STORE *sdStore,
                                                           TID ex_time);
     /*
     **  printMatrElemData()
     ** prints to output file info about precalculation
     ** of <SD'|OP|SD>
     */

        /**** End: function declarations ****/


                // The function definitions
     /* 1.
     ** The entrance function
     **       storeSD_MatrElemInMem()
     ** calculates and stores diagonal <SD|OP|SD> and as
     ** many as possible of the nondiagonal <SD()'|OP|SD()>
     ** The interaction is defined through typeInt
     */

void storeSD_MatrElemInMem(SP_BAS *spBas,SD_BAS *sdBas, SD_STORE *sdStore)
{
  TID  ex_time;

  if(sdStore->diagMemCalc == YES) {

    // calculate and store <SD|OP(typeInt)|SD> in memory 

    if(sdStore->typeInt == VEFF_INT3) {
      threePartDiagSD_MatrElem(sdBas, &sdStore->op, (float *)sdStore->freePermMemPtr);
    }
    else {
      twoPartDiagSD_MatrElem(sdBas, &sdStore->op,(float *)sdStore->freePermMemPtr);
    }
  } // end diagMemCalc

  if(sdStore->memCase == 1)    {

    // calculate and store <SD|OP(typeInt)|SD> in memory 

    sdStore->freePermMemPtr  += sdBas->tot_dimSD * sizeof(float);
    sdStore->freePermMemSize -= (UL)(sdBas->tot_dimSD*sizeof(float));

    time_step(10,1);  // start clock
    nondiagMemMatrElem(spBas, sdBas, sdStore);
    ex_time = time_step(10,2);  // read clock

    printMatrElemData(sdBas, STORE_MEMORY, sdStore, ex_time);
    
  }

} // End: function idStoreSD_MatrElemInMem()

     /* 2.
     ** The module entrance function 
     **       storeSD_MatrElemOnFile()
     ** calculates and stores diagonal if(diagCalc == NO)
     ** and as many as possible of the nondiagonal SD
     ** matrix elements
     **              <SD'()|OP|SD()>
     */

void storeSD_MatrElemOnFile(SP_BAS *spBas, SD_BAS *sdBas, SD_STORE *sdStore)
{
  TID      ex_time;

  if(sdStore->diagMemCalc == NO) {

    // calculate and store <SD|OP(typeInt)|SD> pn file 

    if(sdStore->typeInt == VEFF_INT3) {
      threePartDiagSD_MatrElem(sdBas, &sdStore->op, (float *)sdStore->windMemPtr);
    }
    else {
      twoPartDiagSD_MatrElem(sdBas, &sdStore->op,(float *)sdStore->windMemPtr);
    }

   writeDiagSD_MatrElem(sdBas, sdStore);

  } // end diag Calc on file

       /*
       ** if not all nondiag <SD'|OP|SD>
       ** are calculated and stored in memory
       ** calculate and store on file as many as 
       ** possible of the remaining <SD'|OP|SD> 
       */

  if(   (sdStore->lastSD_mem_asym < (sdBas->numSD[0] - 1))
     || (   (sdStore->MJ == 0)
         && (sdStore->lastSD_mem_sym < (sdBas->tot_dimSD - 1)))) {

    time_step(10,1);  // start clock
    nondiagFileMatrElem(spBas, sdBas, sdStore);
    ex_time = time_step(10,2);  // read clock

    printMatrElemData(sdBas, STORE_FILE, sdStore,ex_time);
  }

} // End: function storeSD_MatrElemOnFile()
   
     /*
     ** The function 
     **         twoPartDiagSD_MatrElem()
     ** calculates all diagonal matrix elements 
     **           <SD()|OP()| SD()> 
     ** where the two-particle matrix elements of OP() are stored in
     ** op_int->id_diag[] using the formula
     **        op(k,l) = ((num - k) * k)/2 + (l - 1)
     ** where num = 2 * m_orb - 3. 
     ** The result is stores in diag_SD[]. 
     ** The result is stored in memory pointed to by sdStore->freePermMemPtr.
     */

static void twoPartDiagSD_MatrElem(SD_BAS *sdBas, MATR_OP *op_int,
                                                    float *diag_ptr)
{
  char      *func = {"twoPartDiagSD_MatrElem(): "};
   int      num_part, *occ, *occ_ptr, par, orb, num, numSD, num_part_1,
            k, l, *k_ptr,*l_ptr;
   ULL      *sd, sd_value, pos;
   double   *table, *table_ptr, value; 

   table = op_int->id_diag;  /* <|H()|> */
   occ = MALLOC(sdBas->part, int, func, "occ[]"); // local memory

   num_part  = sdBas->part; // initialization
   num_part_1 = num_part - 1;
   num        = (sdBas->numm_orb << 1) - 3;
   numSD      = sdBas->tot_dimSD; 
   sd         = sdBas->SD;
   do {
    sd_value = *(sd++);
     occ_ptr = occ;
     for(par = 0, orb = 0, pos = ULL_ONE; par < num_part; par++,
                                                 orb++, pos <<= 1)  {
       for(;!(sd_value & pos); orb++, pos <<= 1); // particle found
       *(occ_ptr++) = orb; // save orbit
     }
     value = D_ZERO;  // calculate <SD(Z)|OP(pp)|SD(Z)> */
     k     = num_part_1;
     k_ptr = occ;
     do {    // k-particle loop
       table_ptr = table + (((num - (*k_ptr)) * (*k_ptr)) >> 1) -1;
       l     = k;
       l_ptr = (++k_ptr);
       do { // l-particle loop */ 
	 value += *(table_ptr + (*(l_ptr++))); 
       } while(--l);  // end l particle loop 
     } while(--k);  // end k particle loop
     *(diag_ptr++) = (float) value; 
   } while(--numSD); // end  loop through all |SD()>

   free(occ); // release local memory

} // End: function  twoPartDiagSD_MatrElem()

     /*
     ** The function 
     **         threePartDiagSD_MatrElem()
     ** calculates all diagonal matrix elements 
     **           <SD()|op_veff(three_part)()| SD()> 
     ** where the three-particle matrix elements of OP() are stored in
     ** op_veff->id_diag[] using the formula
     **  op(k,l,m) =  k*(num1 - k*(num2 -k))/6 + (l*(num3 - l))/2
     **             - num_orb + m     
     ** where num1 = 3*num_orb^2 -12*num_orb + 11
     **       num2 = 3*num_orb - 6
     **       num3 = 2*num_orb - 3
     ** The result is stored in memory pointed to by sdStore->freePermMemPtr.
     */

static void threePartDiagSD_MatrElem(SD_BAS *sdBas, MATR_OP *op_veff, 
                                                    float *diag_ptr)
{
  char    *func = {"idDiagSD_elem(): "};
   int    num_part, *occ, *occ_ptr, par, orb, numSD, num_part_2,
          num1, num2,num3, num_orb, k, l, m,  *k_ptr,*l_ptr, *m_ptr;
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

   numSD = sdBas->tot_dimSD; 
   sd    = sdBas->SD;
   do {
    sd_value = *(sd++);
     occ_ptr = occ;
     for(par = 0, orb = 0, pos = ULL_ONE; par < num_part;
                                   par++, orb++, pos <<= 1)  {
       for(;!(sd_value & pos); orb++, pos <<= 1);// particle found
       *(occ_ptr++) = orb; // save orbit
     }
     value = D_ZERO; // <SD|OP()|SD> = 0
     k     = num_part_2;
     k_ptr = occ;
     do {      // k-particle loop
       table_k_ptr =   table 
                     + ((*k_ptr)*(num1 - (*k_ptr)*(num2 - (*k_ptr))))/6;
       l     = k;
       l_ptr = (++k_ptr);
       do {    // l-particle loop
	 table_kl_ptr = table_k_ptr + ((*l_ptr)*(num3 - (*l_ptr)))/2;
	 m = l;
	 m_ptr = (++l_ptr);
	 do {  //  m-particle loop
	   value += *(table_kl_ptr + (*(m_ptr++) - num_orb));
	 } while(--m); // end m particle loop 
       } while(--l);   // end l particle loop 
     } while(--k);     // end k particle loop
     *(diag_ptr++) = (float) value; 
   } while(--numSD); // end  loop through all |SD()>

   free(occ);  // release local memory

} // End: function  threePartDiagSD_MatrElem()

   /*
   ** The function 
   **         writeDiagSD_MatrElem()
   ** create a new file "file_name" and write diagonal
   ** matrix elements <SD()|OP|SD()>
   ** to file: 
   **     a) Header elements: num_asymSD or num_nosymSD  and num_symSD
   **     b) list of matrix elements 
   */

static void writeDiagSD_MatrElem(SD_BAS *sdBas, SD_STORE *sdStore) 
{
  char    filename[ONE_LINE];
  int     header[2], num_elem;
  float   *mem_ptr;
  FILE    *file_ptr;

  sprintf(filename,"%s%s",sdStore->title,DIAG);
  if((file_ptr = fopen(filename,"wb")) == NULL) {
    printf("\n\nError in function writeDiagSD_MatrElem():");
    printf("\nNot allowed to create the file %s to store diag_elem\n\n",
                                                              filename);
    exit(1);
  }
  header[0] = sdBas->numSD[0];  // write header elements 
  header[1] = sdBas->numSD[1];

  if( fwrite((const void *)header,(size_t) sizeof(int), 2, file_ptr) != 2) { 
    printf("\n\nError in function writeDiagSD_MatrElem():");
    printf("\nIn writing the header elements to file %s\n\n", filename);
    exit(1);
  }
  // write list of diag matrix elements

  mem_ptr = (float *)sdStore->windMemPtr;
  num_elem =  header[0] +  header[1];
  if(fwrite((const void *)mem_ptr,(size_t) sizeof(float), (size_t) num_elem,
                                                 file_ptr) != num_elem) { 
    printf("\n\nError in function writeDiagSD_MatrElem():");
    printf("\nIn writing %d diag matrix elements to file %s\n\n", 
                                              num_elem, filename);
    exit(1);
  }

  fclose(file_ptr);       // close file

} // End: function writeDiagSD_MatrElem()

     /*
     ** The function
     **           nondiagMemMatrElem()
     ** calculates and stores in memory as many as
     ** possible of the nondiagonal SD matrix elements
     **                <SD()'|VEFF|SD()>
     */

void nondiagMemMatrElem(SP_BAS *spBas, SD_BAS *sdBas, SD_STORE *sdStore)
{
  HEAD   head;

  // calculate and store in memory nondiag matrix elements

  if(!(sdBas->MJ)) { // Time-reversal symmetry 

    // calculate and store in memory nondiag asym matr elem 

    if(sdBas->numSD[0])  {
      head.typeSD  = 0;
      head.startSD = 0;
      head.endSD   = sdBas->numSD[0];
      head.numSD   = 0;
      head.memSize = sdStore->freePermMemSize;
      head.memPtr  = sdStore->freePermMemPtr;

      if(sdStore->typeInt == VEFF_INT3) {
	threePartNondiagSD_Calc(&head, spBas, sdBas, &sdStore->op);
      }
      else {
	twoPartNondiagSD_Calc(&head, spBas, sdBas, &sdStore->op);
      }

      sdStore->freePermMemPtr   +=  head.numSD * sizeof(STORE_F);
      sdStore->freePermMemSize  -=  head.numSD * sizeof(STORE_F);

      if(head.numSD > 0) {
	sdStore->lastSD_mem_asym = head.lastSD;
      }
      sdStore->storedSD_mem_elem = head.numSD;
    }

    // calculate and store in memory nondiag sym matr elem

    if(sdBas->numSD[1]) {
      head.typeSD  = 1;
      head.startSD = sdBas->numSD[0];
      head.endSD   = sdBas->tot_dimSD;
      head.numSD   = 0;
      head.memSize = sdStore->freePermMemSize;
      head.memPtr  = sdStore->freePermMemPtr;

      if(sdStore->typeInt == VEFF_INT3) {
	threePartNondiagSD_Calc(&head, spBas, sdBas,&sdStore-> op);
      }
      else {
	twoPartNondiagSD_Calc(&head, spBas, sdBas, &sdStore->op);
      }
      if(head.numSD > 0) {
	sdStore->lastSD_mem_sym = head.lastSD;
      }
      sdStore->storedSD_mem_elem += head.numSD;
    }
  }
  else { // no time-reversal symmetry 

    // calculate and store in memory nondiag nosym matr elem

    head.typeSD  = 2;
    head.startSD = 0;
    head.endSD   = sdBas->numSD[0];
    head.numSD   = 0;
    head.memSize = sdStore->freePermMemSize;
    head.memPtr  = sdStore->freePermMemPtr;

    if(sdStore->typeInt == VEFF_INT3) {
      threePartNondiagSD_Calc(&head, spBas, sdBas, &sdStore->op);
    }
    else {
      twoPartNondiagSD_Calc(&head, spBas, sdBas, &sdStore->op);
    }
    if(head.numSD > 0) {
      sdStore->lastSD_mem_asym = head.lastSD;
    }
    sdStore->storedSD_mem_elem = head.numSD;

  }  // end nosym case

} // End: function  nondiagMemMatrElem()

     /*
     ** The function
     **       nondiagFileMatrElem()
     ** calculates and stores on file as many as possible
     ** of the nondiagonal SD matrix elements
     **                <SD()'|OP|SD()>
     ** not previously calculated and stored in memory
     ** The interaction Op is defined through typeInt
    */

void  nondiagFileMatrElem(SP_BAS *spBas,SD_BAS *sdBas, SD_STORE *sdStore)
{
  HEAD   head;

  // calculate and store on file <SD'|typeInt|SD>

  sdStore->storedSD_file_elem = 0;  // initialization
  sdStore->sect_file_no = 0;
  if(new_section_file(sdStore) == FALSE) return; // no more sect files

  head.memSize = sdStore->windMemSize;  // memory initialization
  head.memPtr  = sdStore->windMemPtr;

  if(!(sdBas->MJ)) { // time-reversal symmetry 
    if(sdBas->numSD[0]) { // <SD'|typeInt|SD_asym>
      head.typeSD  = 0;
      head.startSD = sdStore->lastSD_mem_asym + 1;
      head.endSD   = sdBas->numSD[0];
      while(head.startSD < head.endSD) {
	head.numSD = 0;
	if(sdStore->sect_file_size < sdStore->windMemSize) {
	  head.numSD   = -2;
	  writeNondiagSD_MatrElem(&head,sdStore);
 	  if(new_section_file(sdStore) == FALSE) {
	    head.numSD   = -1;
	    writeNondiagSD_MatrElem(&head,sdStore);
	    return;
	  }
	}
	if(sdStore->typeInt == VEFF_INT3) {
	  threePartNondiagSD_Calc(&head, spBas, sdBas,&sdStore-> op);
	  if(head.numSD == -1)  {
	    writeNondiagSD_MatrElem(&head,sdStore);
	    return;
	  }
	}
	else {
	  twoPartNondiagSD_Calc(&head, spBas, sdBas,&sdStore-> op);
	  if(head.numSD == -1)  {
	    writeNondiagSD_MatrElem(&head,sdStore);
	    return;
	  }
	}
	writeNondiagSD_MatrElem(&head,sdStore);
	sdStore->storedSD_file_elem += head.numSD;
	head.startSD = head.lastSD + 1;
      } // end while()

     //  terminate <SD'|OP|sd_asym> file storage


/***********************/
      head.numSD   = -1;
      writeNondiagSD_MatrElem(&head,sdStore);
/*********************/

      sdStore->lastSD_file_asym = head.startSD - 1;

    } // end asym case 

    // <SD'|typeInt|SD_sym>

    if(sdBas->numSD[1])  {  // <SD'|typeInt|SD_ sym>

      /***********************************/
      if(new_section_file(sdStore) == FALSE) return;  // no more sect files
/****************************/

      head.typeSD  = 1;
      head.startSD = MAX((sdStore->lastSD_mem_sym + 1),sdBas->numSD[0]);
      head.endSD   = sdBas->tot_dimSD;
      while(head.startSD < head.endSD) {
	head.numSD = 0;
	if(sdStore->sect_file_size < sdStore->windMemSize) {
	  head.numSD   = -2;
	  writeNondiagSD_MatrElem(&head,sdStore);
 	  if(new_section_file(sdStore) == FALSE) {
	    head.numSD   = -1;
	    writeNondiagSD_MatrElem(&head,sdStore);
	    return;
	  }
	}
	if(sdStore->typeInt == VEFF_INT3) {
	  threePartNondiagSD_Calc(&head, spBas, sdBas,&sdStore-> op);
	  if(head.numSD == -1)  {
	    writeNondiagSD_MatrElem(&head,sdStore);
	    return;
	  }
	}
	else {
	  twoPartNondiagSD_Calc(&head, spBas, sdBas,&sdStore-> op);
	  if(head.numSD == -1)  {
	    writeNondiagSD_MatrElem(&head,sdStore);
	    return;
	  }
	}
	writeNondiagSD_MatrElem(&head,sdStore);
	sdStore->storedSD_file_elem += head.numSD;
	head.startSD = head.lastSD + 1;
      } // end while()


      //  terminate <SD'|OP|sd_sym> file storage

      head.numSD   = -1;
      writeNondiagSD_MatrElem(&head,sdStore);

      sdStore->lastSD_file_sym = head.startSD - 1;

    } // end sym case 

  } // end sdBas->MJ == 0

  else { // no time-reversal symmetry 
    if(sdBas->numSD[0])  {   // <SD'|typeInt|SD_nosym>
      head.startSD  = sdStore->lastSD_mem_asym + 1;
      head.typeSD  = 2;
      head.startSD  = sdStore->lastSD_mem_asym + 1;
      head.endSD = sdBas->numSD[0];
      while(head.startSD < head.endSD) {
	head.numSD = 0;
	if(sdStore->sect_file_size < sdStore->windMemSize) {
	  head.numSD   = -2;
	  writeNondiagSD_MatrElem(&head,sdStore);
 	  if(new_section_file(sdStore) == FALSE) {
	    head.numSD   = -1;
	    writeNondiagSD_MatrElem(&head,sdStore);
	    return;
	  }
	}
	if(sdStore->typeInt == VEFF_INT3) {
	  threePartNondiagSD_Calc(&head, spBas, sdBas,&sdStore-> op);
	  if(head.numSD == -1)  {
	    writeNondiagSD_MatrElem(&head,sdStore);
	    return;
	  }
	}
	else {
	  twoPartNondiagSD_Calc(&head, spBas, sdBas,&sdStore-> op);
	  if(head.numSD == -1)  {
	    writeNondiagSD_MatrElem(&head,sdStore);
	    return;
	  }
	}
	writeNondiagSD_MatrElem(&head,sdStore);
	sdStore->storedSD_file_elem += head.numSD;
	head.startSD = head.lastSD + 1;
      } // end while()

      //  terminate <SD'|OP|sd_nosym> file storage

      head.numSD   = -1;
      writeNondiagSD_MatrElem(&head,sdStore);

      sdStore->lastSD_file_asym = head.startSD - 1;

    } // end nosym case 

  } // end sdBas->MJ != 0
} // End: function nondiagFileMatrElem()
  
     /* 
     ** The function
     **      new_section_file()
     ** creates a new section file with a corresponding
     ** section number and increase section file number
     ** if enough file space is available
     ** The function return TRUE if new section file is 
     ** is created, otherwise FALSE
     */

static int new_section_file(SD_STORE *sdStore)
{
  static int   fileBytes;
  int          numMbytes;
  FILE         *file_ptr;

  if( sdStore->sect_file_no > 0) {
    sdStore->tot_file_size -= (fileBytes - sdStore->sect_file_size)/M_BYTES + 1;
  }
  if(sdStore->tot_file_size < (sdStore->windMemSize/M_BYTES + 1)) {

    // create a termination file

    sdStore->sect_file_no++;      // new section file number
    sprintf(sdStore->sect_file,"%s%s%d",sdStore->title, NONDIAG,
	                                  sdStore->sect_file_no);
    if((file_ptr = fopen(sdStore->sect_file,"wb")) == NULL) {
      printf("\n\nError in function new_section_file():");
      printf("\nNot allowed to create the file %s\n", sdStore->sect_file);
      exit(1);
    }
    fclose(file_ptr);

    return FALSE;
  } // no more section files

  sdStore->sect_file_no++;      // new section file number
  sprintf(sdStore->sect_file,"%s%s%d",sdStore->title, NONDIAG,
	                                    sdStore->sect_file_no);
  numMbytes = MIN(sdStore->tot_file_size, MAX_FILE_SIZE);  // in Mbytes
  fileBytes = numMbytes * M_BYTES;
  sdStore->sect_file_size = fileBytes;  // in bytes

             // create section file to store nondiag <..|OP|SD) matrix elements 

  if((file_ptr = fopen(sdStore->sect_file,"wb")) == NULL) {
    printf("\n\nError in function new_section_file():");
    printf("\nNot allowed to create the file %s\n", sdStore->sect_file);
    exit(1);
  }
  fclose(file_ptr);

  return  TRUE;

} // End: function  new_section_file()

   /*
   ** The function 
   **         writeNondiagSD_MatrElem()
   ** write nondiagonal matrix elements <SD()|OP|SD()>
   ** to file: 
   **     a) Header elements
   **     b) list of matrix elements 
   */

static void writeNondiagSD_MatrElem(HEAD *head, SD_STORE *sdStore) 
{
  int       header[3], file_size; 
  STORE_F   *mem_ptr; 
  FILE      *file_ptr;

  if((file_ptr = fopen(sdStore->sect_file,"ab")) == NULL) {
    printf("\n\nError in function writeNondiagSD_MatrElem():");
    printf("\nNot allowed to open the file %s to store nondiag_elem\n\n",
                                                    sdStore->sect_file);
    exit(1);
  }
  // write header elements
  
  header[0] = head->startSD;
  header[1] = head->lastSD;
  header[2] = head->numSD;

  if( fwrite((const void *)header,(size_t) sizeof(int), 3, file_ptr) != 3) { 
    printf("\n\nError in function writeNondiagSD_MatrElem():");
    printf("\nIn writing the header elements to file\n\n");
    exit(1);
  }
  file_size = 3 * sizeof(int);   // number of header elements

        // write list of nondiag matrix elements

  if(header[2] > 0) {
    mem_ptr  = (STORE_F *)sdStore->windMemPtr;
    if(fwrite((const void *)mem_ptr,(size_t) sizeof(STORE_F), (size_t) header[2],
                                                          file_ptr) != header[2]) { 
      printf("\n\nError in function writeNondiagSD_MatrElem():");
      printf("\nIn writing %d nondiag STORE_F matrix elements to file %s\n\n",
                                             header[2], sdStore->sect_file);
      exit(1);
    }
    file_size += header[2] * sizeof(STORE_F);

  }   // end write list of STORE_F matrix elements 

  fclose(file_ptr);

  sdStore->sect_file_size -= file_size; // reduce sect_file_size

} // End: function idWriteNondiagSDMatrElem()

  /*
  ** The function
  **        printMatrElemData()
  ** prints to output file info about precalculation
  ** of <SD'|OP|SD>
  */

void printMatrElemData(SD_BAS *sdBas, int storeType, SD_STORE *sdStore,
                                                          TID ex_time)
{
  char  text[5];
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
  }
  if((file_ptr = fopen(sdStore->result,"a"))== NULL) {
    printf("\n\nError in function lanc-main.c():");
    printf("\nWrong file = %s for the output data\n", sdStore->result);
    exit(1);
  }

  if(storeType == STORE_MEMORY) {
    fprintf(file_ptr,"\n\nStored nondiag <SD'|%s|SD> in memory", text);
    fprintf(file_ptr,"\nNumber of stored <SD'|%s|initSD>: %d",text,
                                     (int)sdStore->storedSD_mem_elem); 
    if(!(sdBas->MJ)) {
      fprintf(file_ptr,"\nUp to and including |initSD_asym> = %d",
                                        sdStore->lastSD_mem_asym + 1);
      fprintf(file_ptr,"   Last |initSD_asym> =  %d", 
		                      sdBas->numSD[0]);
      fprintf(file_ptr,"\nUp to and including |initSD_sym>  = %d", 
	                                 sdStore->lastSD_mem_sym +1);
      fprintf(file_ptr," - Last |initSD_sym>  =  %d", 
                                         sdBas->tot_dimSD);
    }
    else {
      fprintf(file_ptr,"\nUp to and including |initSD> = %d", 
                                        sdStore->lastSD_mem_asym + 1);
      fprintf(file_ptr," - Last |initSD> =  %d",
                                  sdBas->numSD[0]);
    }
    fprintf(file_ptr, "\nTime used : %lld hour %lld min %lld sec\n",
	    ex_time.hour, ex_time.min, ex_time.sec);
  } // end memory storage
  else if(storeType == STORE_FILE) {
    fprintf(file_ptr,"\n\nStored nondiag  <SD'|%s|SD>stored on file",text);
    fprintf(file_ptr,"\nNumber of stored <SD'|%s|initSD>: %d", text,
                               (int)sdStore->storedSD_file_elem); 
    if(!(sdBas->MJ)) {
      fprintf(file_ptr,"\nUp to and including |initSD_asym> = %d", 
                                       sdStore->lastSD_file_asym + 1);
      fprintf(file_ptr," - Last |initSD_asym> =  %d",
                                          sdBas->numSD[0]);
      fprintf(file_ptr,"\nUp to and including |initSD_sym>  = %d", 
	                              sdStore->lastSD_file_sym +1);
      fprintf(file_ptr," - Last |initSD_sym>  =  %d", 
                                     sdBas->tot_dimSD);
    }
    else {
      fprintf(file_ptr,"\nUp to and including |initSD> = %d", 
                                sdStore->lastSD_file_asym + 1);
      fprintf(file_ptr," - Last |initSD> =  %d",
                                 sdBas->numSD[0]);
    }
    fprintf(file_ptr, "\nTime used : %lld hour %lld min %lld sec\n",
    	                  ex_time.hour, ex_time.min, ex_time.sec);
  }  // end calculate and store remaining <SD'|OP|SD> on file

  fclose(file_ptr);

} // End: function printMatrElemData()
