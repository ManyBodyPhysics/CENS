
/********************  File PAR-id-nondiag-calc.c *************/

#include "PAR-shell-model.h"

// The file contains two entrance functions 

    /*  
    ** 1. The module entrance function
    **           twoPartNondiagSD_Calc()
    ** calculates nondiagonal identical particle matrix elements 
    **             <SD()'|OP|SD()>
    ** from |SD()> = |head->startSD> and store as many as possible
    ** in head->memSize. |SD()> are many-particle basis state.
    ** Pointers to all two-particle matrix elements are stored in 
    ** structure id_diag by the formula
    **        op(k,l) = ((num - k) * k)/2 -1 +l
    ** where num = 2 * m_orb - 3.
    ** All <SD'|OP|SD< are stored temporary in memory
    ** The three first elements are reserved for special purpose 
    ** Element two and three are here calculated to 
    ** STORE_F elem2.final = initSD number
    ** STORE_F elem3.final = number of |initSD> groups in 
    ** current mem block
    */

    /*
    ** 2. The module entrance function
    **            threePartNondiagSD_Calc()
    ** calculates all nondiagonal identical particle matrix elements 
    **             <SD()'|op(veff(three_part)| SD()> 
    ** where |SD()> is a many-particle basis state.
    ** Pointers to all three-particle matrix elements are stored in 
    ** structure MATR_OP op_veff->id_nondiag_table by the formula
    **  op(k,l,m) =  (num1 - k*(num2 -k))/6 + (l*(num3 - l))/2
    **             - num_orb + m     
    ** where num1 = 3*num_orb^2 -12*num_orb + 11
    **       num2 = 3*num_orb - 6
    **       num3 = 2*num_orb - 3
    ** All <SD'|OP|SD< are stored temporary in memory
    ** The three first elements are reserved for special purpose 
    ** Element two and three are here calculated to 
    ** STORE_F elem2.final = initSD number
    ** STORE_F elem3.final = number of |initSD> groups in 
    ** current mem block
    */
	    /**** local function declarations ****/

static int id_sort_nondiag(const STORE_F *one, const STORE_F *two);
    /*
    ** is a utility function for the library function qsort() in order to
    ** sort nondiagSD matrix elements of type STORE store[]
    ** after increasing store[].final.
    */

static int id_compress_nondiag_SD_elem(int num, STORE_F *store, int loopSD);
    /*
    ** takes num calculated nondiag SD elem which are sorted
    ** after increasing store[].final and add together all
    ** elements  with same store[].final.
    ** Note: The function must have num > 1;
    ** The final number of elements is returned.
    */

static void id_binary_search_nondiag_asym_elem(HEAD *head, 
				        const SD_BAS *sd_bas);
    /*
    ** performs a binary search on all time-reversal
    ** asymmetric |final_SD()> and replace them with the
    ** corresponding sequence number.
    */

static void id_binary_search_nondiag_sym_elem(HEAD *head,
				         const SD_BAS *sd_bas);
    /*
    ** performs a binary search on all time-reversal symmetric 
    ** |final_SD()> and replace them with the corresponding
    ** sequence number.
    */

static void id_binary_search_nondiag_nosym_elem(HEAD *head,
				          const SD_BAS *sd_bas);
    /*
    ** and replace them with the corresponding sequence number.
    ** Then all calculated  <SD()_final|OP|SD()_init> are 
    ** stored on file.
    */

                /**** End: function declarations ****/

     // The function definitions

    /*
    ** The entrance function
    **           twoPartNondiagSD_Calc()
    ** calculates nondiagonal identical particle matrix elements 
    **             <SD()'|OP|SD()>
    ** from |SD()> = |head->startSD> and store as many as possible
    ** in head->memSize. |SD()> are many-particle basis state.
    ** Pointers to all two-particle matrix elements are stored in 
    ** structure id_diag by the formula
    **        op(k,l) = ((num - k) * k)/2 -1 +l
    ** where num = 2 * m_orb - 3.
    ** All <SD'|OP|SD< are stored temporary in memory
    ** The three first elements are reserved for special purpose 
    ** Element two and three are here calculated to 
    ** STORE_F elem2.final = initSD number
    ** STORE_F elem3.final = number of |initSD> groups in 
    ** current mem block
    */

void twoPartNondiagSD_Calc(HEAD *head, SP_BAS *sp_bas,
   		         SD_BAS *sd_bas, MATR_OP *op_int)
{
  char    *func = {"id_nondiag_SD_calc(): "};
  int     num_part, num_part_1, num, num_stored,
          number, k, l, orb,  *occ, kl_phase, phase, count,
          *k_ptr, *l_ptr, loopSD, maxMemElem;
  ULL     *initSD, pos, sd_k, sd_kl, new_sd, low, high,
          topBit;
  STORE_F *matr_ptr, *curr_ptr;
  ID_INT  **table, **tab_k, *ij_ptr;

  int newNum;

  topBit     = ULL_ONE << (8 * sizeof(ULL) - 1);
  matr_ptr   = (STORE_F *)head->memPtr;
  maxMemElem = head->memSize; // available memory storage

  table      = op_int->id_nondiag_table;  // <|H|>
  num_part   = sd_bas->part;
  num_part_1 = num_part - 1;
  num        = (sd_bas->numm_orb << 1) - 3;
  num_stored = 0;  // initialization

  occ = MALLOC(sd_bas->part, int, func, "occ[]"); // local memory

  initSD     = sd_bas->SD + head->startSD;
  for(loopSD = head->startSD; loopSD < head->endSD;
                     loopSD += NumProcs, initSD += NumProcs) {
    number   = 0;   // initialization
    curr_ptr = matr_ptr;
    k_ptr    = occ;
    for(k = 0, orb = 0, pos = ULL_ONE; k < num_part; 
                              k++, orb++, pos <<= 1)  {
      for(;!((*initSD) & pos); orb++, pos <<= 1); // particle found
      *k_ptr = orb; //save orbit
      k_ptr++;
    }

    k_ptr  = occ;
    for(k = 0; k < num_part_1; k++, k_ptr++) { 
      tab_k    = table + (((num - (*k_ptr)) * (*k_ptr)) >> 1) -1;
      sd_k     = ((*initSD) ^ (ULL_ONE << (*k_ptr)));
      l_ptr    = k_ptr + 1;
      kl_phase = +1;
      for(l = k + 1; l <= num_part_1; l++, l_ptr++, kl_phase = -kl_phase) {
	if((ij_ptr = *(tab_k+(*l_ptr))) == NULL_PTR) continue; 
	sd_kl = sd_k ^ (ULL_ONE << *l_ptr); 
	for(  ; ij_ptr->one; ij_ptr++) {

	  for(; sd_kl & ij_ptr->one; ij_ptr++);  
	  if(!ij_ptr->one) break;  // no more ij-pairs 
	  new_sd = sd_kl ^ ij_ptr->one; // contr found, new SD-state

             /* 
	     ** check that |new_sd> has allowed number of 
	     ** particles in all j-orbits - only active if
	     **  0 < allowed number < 2*j + 1
             */
	  
	  if(sp_bas->mask.num > 0) {    /* check masking */
	    for(count = 0; count < sp_bas->mask.num; count++) {
	      high = new_sd & sp_bas->mask.list[count];
	      phase = 0;
	      while(high) {
		high &= high - 1;
		phase++;
	      }
	     if(  (phase > sp_bas->mask.lim[2*count + 1])
		  ||(phase < sp_bas->mask.lim[2*count])) {
               break; // |new_sd> not accepted
	     }
	    }
            if(count < sp_bas->mask.num) continue; 

	  } //end check masking

	  low   = new_sd & ij_ptr->two; // permutation phase
	  phase = kl_phase;
	  while(low) {
	    low &= low - 1;
	    phase = -phase;
	  }
	  curr_ptr->value     = (float)(phase * ij_ptr->val);
	  curr_ptr->final = new_sd;
	  curr_ptr++;
	  number++;

            /* 
	    ** if not more memory is available
	    ** to save another <SD|OP|SD> 
	    ** no mote contribution from interaction
	    */ 

	  if(number >= maxMemElem) {

	           //  terminate the caculation: No more memory

	      FILE      *filePtr;
	      int   testNum;

	      if((filePtr = fopen(TEST_FILE_NAME,"a"))== NULL) {
		printf("\n\nRank%d: Error in ????????:",Rank);
		printf("\nWrong file = %s for the output data file\n",TEST_FILE_NAME);
		MPI_Abort(MPI_COMM_WORLD,Rank);
	      }
	      fprintf(filePtr,"\n\n\n\n Message from function  twoPartNondiagSD_Calc():");
	      fprintf(filePtr,"\n No more storage for nondiag matrix elements.");
              fprintf(filePtr,"\n head->typeSD = %d  loopSD = %d num_stored = %d",
		                                  head->typeSD, loopSD,num_stored); 
	      MPI_Abort(MPI_COMM_WORLD,Rank);

	  } // End control test

	}  // run through all two-body matr elem for given (kl) 

      } // end l-loop

    } // end k-loop

         /* 
	 ** sort and compress nondiag matrix elements for 
	 ** each |init_SD> after increasing |final_SD>
	 */

    if(number > 1) {
      qsort(matr_ptr, (size_t) number, (size_t)sizeof(STORE_F),
	    (int(*)(const void *, const void *)) id_sort_nondiag);
    }
       // terminate |initSD> block of matr elem by an extra element


    curr_ptr        = matr_ptr + number;
    curr_ptr->value = D_ZERO;
    curr_ptr->final = topBit;
    number++;

    number = id_compress_nondiag_SD_elem(number, matr_ptr,loopSD);

    num_stored += number;
    matr_ptr   += number; // mem position for next <SD'|OP|SD> contribution
    maxMemElem -= number; // free memory storage for next loopSD
 
  } // end loop through all |initSD>



   free(occ); // release temporary memory

   head->numSD    = num_stored;
   head->lastSD   = loopSD - NumProcs;

   if(num_stored) {

     switch(head->typeSD) {

     case 0: id_binary_search_nondiag_asym_elem(head, sd_bas);
             break;
     case 1: id_binary_search_nondiag_sym_elem(head, sd_bas);
             break;
     case 2: id_binary_search_nondiag_nosym_elem(head, sd_bas);
             break;
    } // end switch()
   }

   return; 
 
} // End: function twoPartNondiagSD_Calc()

    /*
    ** The entrance function
    **            threePartNondiagSD_Calc()
    ** calculates nondiagonal identical particle matrix elements 
    **             <SD()'|op(veff(three_part)| SD()> 
    ** from |SD()> = |head->startSD> and store as many as possible
    ** in head->memSize. |SD()> are many-particle basis state.
    ** Pointers to all three-particle matrix elements are stored in 
    ** structure MATR_OP op_veff->id_nondiag_table by the formula
    **  op(k,l,m) =  (num1 - k*(num2 -k))/6 + (l*(num3 - l))/2
    **             - num_orb + m     
    ** where num1 = 3*num_orb^2 -12*num_orb + 11
    **       num2 = 3*num_orb - 6
    **       num3 = 2*num_orb - 3
    ** All <SD'|OP|SD< are stored temporary in memory
    ** The three first elements are reserved for special purpose 
    ** Element two and three are here calculated to 
    ** STORE_F elem2.final = initSD number
    ** STORE_F elem3.final = number of |initSD> groups in 
    ** current mem block
    */

void  threePartNondiagSD_Calc(HEAD *head, SP_BAS *sp_bas,
		        SD_BAS *sd_bas, MATR_OP *op_int)
{
  char      *func = {"i threePartNondiagSD_elem(): "};
  int       num_part, num_part_2,  num_orb, num1, num2, num3,
            num_stored, number,
            k, l, m, orb,*occ, k_phase, kl_phase, klm_phase, phase, count, 
            *k_ptr,*l_ptr,*m_ptr,loopSD, max_mem_elem, klm_elem_lim,
            maxMemElem;
  ULL       *initSD, pos, sd_k, sd_kl, sd_klm, new_sd, topBit,low, high;
  STORE_F   *matr_ptr,*curr_ptr;
  ID_INT    **table,**tab_k,**tab_kl,*hij_ptr;
  int       sumInitMJ;

  topBit     = ULL_ONE << (8 * sizeof(ULL) - 1);
  matr_ptr   = head->memPtr;
  maxMemElem = head->memSize; // available memory storage
  table      = op_int->id_nondiag_table;  // <|H|>
  num_part   = sd_bas->part;
  num_part_2 = num_part - 2;

  num_orb = sd_bas->numm_orb;
  num1    = 3*num_orb*num_orb -12*num_orb + 11;
  num2    = 3*num_orb - 6;
  num3    = 2*num_orb - 3;

  num_stored = 0;
  number     = 0; // initialization

  occ = MALLOC(sd_bas->part, int, func, "occ[]");

  initSD     = sd_bas->SD + head->startSD;

  for(loopSD = head->startSD; loopSD < head->endSD;
             loopSD += NumProcs, initSD += NumProcs) {
    number   = 0; // initialization
    curr_ptr = matr_ptr;
    k_ptr    = occ;
    for(k = 0, orb = 0, pos = ULL_ONE; k < num_part; 
                             k++, orb++, pos <<= 1)  {
      for(;!((*initSD) & pos); orb++, pos <<= 1); 
      *k_ptr = orb;
      k_ptr++;
    }
  
    k_ptr = occ;
    for(k = 0; k < num_part-2; k++, k_ptr++) { 
      k_phase = PHASE(k);
      tab_k    = table + ((num1 - (*k_ptr)*(num2 - (*k_ptr)))*(*k_ptr))/6;
      sd_k     = ((*initSD) ^ (ULL_ONE << (*k_ptr)));
      l_ptr    = k_ptr + 1;
      for(l = k + 1; l <= num_part_2; l++, l_ptr++) {
	kl_phase = k_phase *  PHASE(l - 1);
	tab_kl =  tab_k + ((*l_ptr)*(num3 - (*l_ptr)))/2;
	sd_kl =  (sd_k ^ (ULL_ONE << (*l_ptr)));
	m_ptr = l_ptr +1;
	for(m = l + 1; m < num_part; m++, m_ptr++) {
	  klm_phase = kl_phase*PHASE(m - 2);
	  if((hij_ptr = *(tab_kl + (*m_ptr) - num_orb)) == NULL_PTR) continue; 
	  sd_klm = sd_kl ^ (ULL_ONE << *m_ptr); 

	  // do  {      // run through all ij-pair

	  //****************   NY ekstra kode **********

	  if( hij_ptr->one == ULL_ZERO) continue;

	  for(  ; hij_ptr->one; hij_ptr++) { 

	    for(; sd_klm & hij_ptr->one; hij_ptr++);  

	      if(hij_ptr->one == ULL_ZERO)  break;

	      //  if(!hij_ptr->one) break; // no more ij-triplwts

	      new_sd = sd_klm ^ hij_ptr->one; // contr found, new SD-state



             /* 
	     ** check that |new_sd> has allowed number of particles in all
             ** j-orbits - only active if 0 < allowed number < 2*j + 1
             */

	    if(sp_bas->mask.num > 0) {    /* check masking */
	      for(count = 0; count < sp_bas->mask.num; count++) {
		high = new_sd & sp_bas->mask.list[count];
		phase = 0;
		while(high) {
		  high &= high - 1;
		  phase++;
		}
		if(  (phase > sp_bas->mask.lim[2*count + 1])
		   ||(phase < sp_bas->mask.lim[2*count])) break;
	      }
	      if(count < sp_bas->mask.num) continue; 

	    } //end check masking

	    low   = new_sd & hij_ptr->two; // permutation phase
	    phase = klm_phase;
	    while(low) {
	      low &= low - 1;
	      phase = -phase;
	    }
	    curr_ptr->value = (float)(phase * hij_ptr->val);
	    curr_ptr->final = new_sd;
	    curr_ptr++;
	    number++;

              /* 
	      ** if not more memory is available
	      ** to save another <SD|OP|SD> 
	      ** no mote contribution from interaction
	      */ 
	    
	    if(number >= maxMemElem) {

	           //  terminate the caculation: No more memory

	      FILE      *filePtr;
	      int   testNum;

	      if((filePtr = fopen(TEST_FILE_NAME,"a"))== NULL) {
		printf("\n\nRank%d: Error in ????????:",Rank);
		printf("\nWrong file = %s for the output data file\n",TEST_FILE_NAME);
		MPI_Abort(MPI_COMM_WORLD,Rank);
	      }
	      fprintf(filePtr,"\n\n\n\n Meeasge from function  twoPartNondiagSD_Calc():");
	      fprintf(filePtr,"\n No more storage for nondiag matrix elements.");
              fprintf(filePtr,"\n head->typeSD = %d  loopSD = %d num_stored = %d",
		                                  head->typeSD, loopSD,num_stored); 
	      MPI_Abort(MPI_COMM_WORLD,Rank);

	    } // End control test

	  }  // run through all two-body matr elem for given (kl) 

	} // end m_loop

      } // end l_loop

    } // end k_loop

    
         /*
	 ** sort and compress nondiag matrix elements for 
	 ** each |init_SD> after increasing |final_SD>
	 */
    
    if(number > 1) {
      qsort(matr_ptr, (size_t) number, (size_t)sizeof(STORE_F),
	    (int(*)(const void *, const void *)) id_sort_nondiag);
    }

    // terminate |initSD> block of matr elem by an extra element


    curr_ptr        = matr_ptr + number;
    curr_ptr->value = D_ZERO;
    curr_ptr->final = topBit;
    number++;

    number = id_compress_nondiag_SD_elem(number, matr_ptr,loopSD);


    num_stored += number;
    matr_ptr   += number; // mem position for next <SD'|OP|SD> contribution
    maxMemElem -= number; // free memory storage for next loopSD
 
  } // end loop through all |initSD>

  free(occ); // release temporary memory

  head->numSD  = num_stored;
  head->lastSD = loopSD - NumProcs;

  if(num_stored) {

    switch(head->typeSD) {

    case 0: id_binary_search_nondiag_asym_elem(head, sd_bas);
            break;
    case 1: id_binary_search_nondiag_sym_elem(head, sd_bas);
            break;
    case 2: id_binary_search_nondiag_nosym_elem(head, sd_bas);
            break;
    } // end switch()
  }
  return;
 
} // End: function threePartNondiagSD_elem()

    /*
    ** The function                         
    **        int id_sort_nondiag()                  
    ** is a utility function for the library function qsort() in order to
    ** sort nondiagSD matrix elements of type STORE store[]
    ** after increasing store[].final.
    */

static int id_sort_nondiag(const STORE_F *one, const STORE_F *two)
{
  if(one->final > two->final)       return +1;
  else  if(one->final < two->final) return -1;
  else                              return  0;
 } // End: function id_sort_nondiag() 

      /*
      ** The function 
      **      id_compress_nondiag_SD_elem()
      ** takes num calculated nondiag SD elem which are sorted
      ** after increasing store[....].final and add together all
      ** elements  with same store[...].final.
      ** All elements with store[...]->value < MATRIX_LIMIT is removed 
      ** Note: The function must have num > 1;
      ** The final number of elements is returned.
      */

static int id_compress_nondiag_SD_elem(int num, STORE_F *store, int loopSD)
{
  char      *func = {"id_compress_nondiag_SD_elem(): "};

  int       newNum, restNum;
  STORE_F   *initPtr, *finalPtr;
  ULL       topBit;

  topBit  = ULL_ONE << (8 * sizeof(ULL) - 1);
  newNum  = 0;
  restNum = num - 1;

  finalPtr = store;
  if(num > 1) {
    initPtr = store + 1;
  }
  else if(num == 1) {
    return  num;
  }
  else {
    FILE      *filePtr;

    if((filePtr = fopen(TEST_FILE_NAME,"a"))== NULL) {
      printf("\n\nRank%d: Error in ????????:",Rank);
      printf("\nWrong file = %s for the output data file\n",TEST_FILE_NAME);
      MPI_Abort(MPI_COMM_WORLD,Rank);
    }

    fprintf(filePtr,"\n\n ERROR in function id_compress_nondiag_SD_elem()");
    fprintf(filePtr,"\n num = %d\n\n\n",num);
    fflush(filePtr);
    fclose(filePtr);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }

  while(restNum > 1) {

    if(finalPtr->final == topBit) {
      newNum++;
      finalPtr++;
      finalPtr->final = initPtr->final;
      finalPtr->value = initPtr->value;
      initPtr++;
      restNum--;
      continue;
    }
    
    // Add together  a group of elem

    while(finalPtr->final == initPtr->final) {
      finalPtr->value += initPtr->value;
      initPtr++;
      restNum--;
    }

    if(fabs(finalPtr->value) > MATRIX_LIMIT) {
      newNum++;
      finalPtr++;
    }

    if(restNum == 1)  {
      restNum--;
      break;
    }

    finalPtr->final = initPtr->final;
    finalPtr->value = initPtr->value;
    initPtr++;
    restNum--;

  }  // end while()



  if(restNum == 1) {

    // check next to last 

    if(fabs(finalPtr->value) > MATRIX_LIMIT) {
      newNum++;
      finalPtr++;
    }
  }

  //  Add last topBit

  finalPtr->final = initPtr->final;
  finalPtr->value = initPtr->value;
  newNum++;


  //  Check last elem equal to topBit

  if(finalPtr->final  != topBit) {
    FILE      *filePtr;

    if((filePtr = fopen(TEST_FILE_NAME,"a"))== NULL) {
      printf("\n\nRank%d: Error in ????????:",Rank);
      printf("\nWrong file = %s for the output data file\n",TEST_FILE_NAME);
      MPI_Abort(MPI_COMM_WORLD,Rank);
    }

    fprintf(filePtr,"\n\n ERROR in function id_compress_nondiag_SD_elem()");
    fprintf(filePtr,"\nLast elem in store() not equal to topBit");
    fprintf(filePtr,"\n num = %d  newNum = %d  restNum = %d\n\n\n",
	    num,newNum,restNum);
    fflush(filePtr);
    fclose(filePtr);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }

  return newNum;

}    // END  function new-compress 

    /*
    ** The function 
    **       id_binary_searchnondiag_asym_elem() 
    ** performs a binary search on all time-reversal
    ** asymmetric |final_SD()> and replace them with the
    ** corresponding sequence number.
    */

static void id_binary_search_nondiag_asym_elem(HEAD *head, 
                                      const SD_BAS *sd_bas)
{
  int       num, loopSD,curr_initSD, newNumElem;
  ULL       topBit, low, high, new_tr_sd, init_tr_sd, 
            search;
  double    j_phase;
  STORE_F   *locPtr; 

  newNumElem   = head->numSD;                 // initializatio
  curr_initSD  = head->startSD;
  topBit       = ULL_ONE << (8 * sizeof(ULL) - 1);
  locPtr       = (STORE_F *)head->memPtr;
  j_phase      = sd_bas->J_type ? 1.0 : -1.0;

  for(loopSD = 0; loopSD <  head->numSD; loopSD++, locPtr++) {

       /*
       ** if(locPtr->final == topBit)    
       **  no number search of current elem
       */

   if(locPtr->final == topBit) {
     curr_initSD += NumProcs;
     continue;
   }
    num       = sd_bas->numm_orb; // calculate |new_tr_sd>
    new_tr_sd = ULL_ZERO;
    low       = ULL_ONE;
    high      = (ULL_ONE<<(num - 1));
    while(num--) {
      if(low & locPtr->final) new_tr_sd |= high;
      low  <<= 1;
      high >>= 1;
    }

    //  CASE A: |init_asym>  ---> 1final_asym>

    if((locPtr->final) < new_tr_sd)  { // |new_asymSD>
      low  = curr_initSD;
      high = sd_bas->numSD[0];

      for( ; ;)  { 
        search = (low + high) >> 1;
	if(locPtr->final < sd_bas->SD[search])      high = search - 1;
	else if(locPtr->final > sd_bas->SD[search]) low  = search + 1;
	else                                        break;
      }
      locPtr->final = search << 2;
    } // end new |asymSD> config


    // CASE B:   |init_asym>  ---> 1final_sym>

    else if(locPtr->final == new_tr_sd){ // |new_symSD> 
      if(!sd_bas->J_type) { // Odd total J - no |new_symSD>
	newNumElem--;
	continue;  
      }
      low  = sd_bas->numSD[0];
      high = sd_bas->tot_dimSD;

      for( ; ; )  {
	search = (low + high) >> 1;
	if(locPtr->final < sd_bas->SD[search])      high = search - 1;
	else if(locPtr->final > sd_bas->SD[search]) low  = search + 1;
	else                                        break;
      }
          /*
	  ** Add contribution from the excluded time-reversed 
	  ** |SD'>s when init_tr_sd  < locPtr->final
	  */
      num        =  sd_bas->numm_orb; // calculate |init_tr_sd>
      init_tr_sd = ULL_ZERO;
      low        = ULL_ONE;
      high       = (ULL_ONE<<(num - 1));
      while(num--) {
	if(low & sd_bas->SD[curr_initSD]) init_tr_sd |= high;
	low  <<= 1;
	high >>= 1;
      }    
      locPtr->final = (search << 2)
               + ((init_tr_sd  < locPtr->final) ? 1 : 0);
    }     // end new |sym_SD[]> config

    //  CASE C:   |newSD> = |new_trSD>

    else  { 
	
        /*
	** new_sd is a time reversed |SD> , binary search
	** in asym sd_bas->SD[]. Add contribution from the
	** excluded  time-reversed component.
	*/

      low  = 0;
      high = sd_bas->numSD[0];

      for( ; ; )  {
        search = (low + high) >> 1;
	if(new_tr_sd < sd_bas->SD[search])      high = search - 1;
	else if(new_tr_sd > sd_bas->SD[search]) low  = search + 1;
	else                                           break;
      }
      num        = sd_bas->numm_orb;// calculate |init_tr_sd>
      init_tr_sd = ULL_ZERO;
      low        = ULL_ONE;
      high       = (ULL_ONE<<(num - 1));
      while(num--) {
	if(low & sd_bas->SD[curr_initSD]) init_tr_sd |= high;
	low  <<= 1;
	high >>= 1;
      }
      locPtr->value *= j_phase;
      locPtr->final  =  (search << 2)
                      + ((init_tr_sd < new_tr_sd) ? 3 : 2);
    } // end new contribution to final tr_asym_SD[]
  
 }  // all |final_SD()> transferred to corresponding number


 head->numSD = newNumElem; // update number of elements

 return;

} // End: function id_binary_search_nondiag_asym_elem()

    /*
    ** The function 
    **       id_binary_search_nondiag_sym_elem() 
    ** performs a binary search on all time-reversal symmetric 
    ** |final_SD()> and replace them with the corresponding
    ** sequence number.
    */

static void id_binary_search_nondiag_sym_elem(HEAD *head,
                                         const SD_BAS *sd_bas)
{
  int      loopSD,num,curr_initSD;
  ULL      topBit, low, high, new_tr_sd, search;
  STORE_F  *locPtr;

  curr_initSD = head->startSD;
  topBit      = ULL_ONE << (8 * sizeof(ULL) - 1);
  locPtr      = (STORE_F *)head->memPtr;

  for(loopSD = 0; loopSD < head->numSD; loopSD++,locPtr++) {

    if(locPtr->final == topBit) {
      curr_initSD += NumProcs;
      continue;
    }

    num       = sd_bas->numm_orb; // calculate |new_tr_sd>
    new_tr_sd = ULL_ZERO;
    low       = ULL_ONE;
    high      = (ULL_ONE<<(num - 1));
    while(num--) {
      if(low & locPtr->final) new_tr_sd |= high;
      low  <<= 1;
      high >>= 1;
    }
    if((locPtr->final) < new_tr_sd)  {// |new_asymSD>

      low  = 0;
      high = sd_bas->numSD[0];
      for ( ; ; )  {   //   while(1) {
        search = (low + high) >> 1;
	if(locPtr->final < sd_bas->SD[search])      high = search - 1;
	else if(locPtr->final > sd_bas->SD[search]) low  = search + 1;
	else                                        break;
      }
      locPtr->final = search << 2;
    } // end new |asymSD> config


    // CASE B:  |initSD_sym>   --> |finalSD_syn>

    else if(locPtr->final == new_tr_sd){ // |new_symSD>

      low  = curr_initSD;
      high = sd_bas->tot_dimSD;
      for(  ;  ; )  {      // while(1) {
	search = (low + high) >> 1;
	if(locPtr->final < sd_bas->SD[search])      high = search - 1;
	else if(locPtr->final > sd_bas->SD[search]) low  = search + 1;
	else                                        break;
      }
      locPtr->final = search << 2;

    }     /* end new |sym_SD[]> config */

    //  CASE C:   |initSD_sym>  __> |finalSD_tr_asym>

    else  {  

        /*
	** |new_sd> is a tr_asym_SD[]>, search in 
	** final state sd_bas->SD[0]. Add contribution 
	** from the excluded time-reversed component.
	*/

      low  = 0;
      high = sd_bas->numSD[0];
      for( ; ;)  {         //   while(1) {
        search = (low + high) >> 1;
	if(new_tr_sd < sd_bas->SD[search])      high = search - 1;
	else if(new_tr_sd > sd_bas->SD[search]) low  = search + 1;
	else                                          break;
      }
      locPtr->final =  (search << 2) + 1;

    } // end new contribution to final tr_asym_SD[]
  
  }  // all |final_SD()> transferred to corresponding number
 

}  // End: function id_binary_search_nondiag_sym_elem()



    /*
    ** The function 
    **       id_binary_search_nondiag_nosym_elem() 
    ** performs a binary search on all nosymmetric |final_SD()>
    ** and replace them with the corresponding sequence number.
    ** Then all calculated  <SD()_final|OP|SD()_init> are 
    ** stored on file.
    */

static void id_binary_search_nondiag_nosym_elem(HEAD *head,
                                       const SD_BAS *sd_bas)
{
  int       loopSD, curr_initSD;
  ULL       topBit, low, high, search;
  STORE_F   *locPtr;

  curr_initSD  = head->startSD; // initialization
  topBit       = ULL_ONE << (8 * sizeof(ULL) - 1);
  locPtr       = (STORE_F *)head->memPtr;

  for(loopSD = 0; loopSD < head->numSD; loopSD++, locPtr++) {

    if(locPtr->final == topBit) {
      curr_initSD += NumProcs;
      continue;
    }
    low  = curr_initSD;
    high = sd_bas->numSD[0];
    for( ; ; ) {            // while(1) {
      search = (low + high) >> 1;
      if(locPtr->final < sd_bas->SD[search])      high = search - 1;
      else if(locPtr->final > sd_bas->SD[search]) low  = search + 1;
      else                                     break;
    }
    locPtr->final = search;

  } // all |final_SD()> transferred to corresponding number

} // End: function id_binary_search_nondiag_nosym_elem()
