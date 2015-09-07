/********************  File id-nondiag-calc.c *************/

#include "shell-model.h"

// The file contains two entrance functions 

    /*  
    ** 1. The module entrance function
    **           twoPartNondiagSD_Calc()
    ** calculates all nondiagonal identical particle matrix elements 
    **             <SD()'|OP| SD()> 
    ** where |SD()> are many-particle basis state.
    ** Pointers to all two-particle matrix elements are stored in 
    ** structure id_diag by the formula
    **        op(k,l) = ((num - k) * k)/2 -1 +l
    ** where num = 2 * m_orb - 3.
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
    ** The nondiag matrix elements are saved in memory pointed
    ** to by sd_store->permMemPtr
    */
	    /**** local function declarations ****/

static int id_sort_nondiag(const STORE_F *one, const STORE_F *two);
    /*
    ** is a utility function for the library function qsort() in order to
    ** sort nondiagSD matrix elements of type STORE store[]
    ** after increasing store[].final.
    */

static int id_compress_nondiag_SD_elem(int num, STORE_F *store);
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
    ** calculates all nondiagonal identical particle matrix elements 
    **             <SD()'|OP| SD()> 
    ** where |SD()> are many-particle basis state.
    ** Pointers to all two-particle matrix elements are stored in 
    ** structure id_diag by the formula
    **        op(k,l) = ((num - k) * k)/2 -1 +l
    ** where num = 2 * m_orb - 3.
    */

void twoPartNondiagSD_Calc(HEAD *head, SP_BAS *sp_bas,
		     SD_BAS *sd_bas, MATR_OP *op_int)
{
  char    *func = {"id_nondiag_SD_calc(): "};
  int     num_part, num_part_1, num, num_stored,
          number, k, l, orb,  *occ, kl_phase, phase, count, *k_ptr, *l_ptr,
          loopSD, max_mem_elem, kl_elem_lim;
  ULL     *initSD, pos, sd_k, sd_kl, new_sd, low, high;
  STORE_F *matr_ptr, *curr_ptr;
  ID_INT  **table, **tab_k, *ij_ptr;

  max_mem_elem = head->memSize/sizeof(STORE_F); // storage 
  matr_ptr     = (STORE_F *)head->memPtr;

  kl_elem_lim  =  ((sd_bas->part * (sd_bas->part -1))/2)
                * op_int->maxNondiagElem;

      /*
      ** memory must exceed number of nondiag matrix 
      ** elements from one |initSD>
      */

  if(max_mem_elem < kl_elem_lim) {

    // return without any calculation

    head->numSD   = -1;
    return;
  }

  table          = op_int->id_nondiag_table;  // <|H|>
  num_part       = sd_bas->part;
  num_part_1     = num_part - 1;
  num            = (sd_bas->numm_orb << 1) - 3;
  num_stored     = 0;
  curr_ptr       = matr_ptr;

  occ = MALLOC(sd_bas->part, int, func, "occ[]"); // local memory

  initSD = sd_bas->SD + head->startSD;

  for(loopSD = head->startSD; loopSD < head->endSD; loopSD++, initSD++) {

     /* 
     ** if not enough memory is available to save a new
     ** set of kl_elem_lim nondiag matrix elements save
     ** current elements on file
     ** terminate calculation
     */ 

    if((max_mem_elem - num_stored) < kl_elem_lim) break;

    number   = 0;   // initialization
    curr_ptr = matr_ptr;
    k_ptr    = occ;

    for(k = 0, orb = 0, pos = ULL_ONE; k < num_part; 
                              k++, orb++, pos <<= 1)  {
      for(;!((*initSD) & pos); orb++, pos <<= 1); // particle found
      *(k_ptr++) = orb; //save orbit
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
	do  {  // run through all ij-pair
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
	  (curr_ptr++)->final = new_sd;
	  number++;
	} while((++ij_ptr)->one); // all two-body matr elem for given (kl)
      } // end l-loop
    } // end k-loop

         /* 
	 ** sort and compress nondiag matrix elements for 
	 ** each |init_SD> after increasing |final_SD>
	 */

    if(number > 1) {
      qsort(matr_ptr, (size_t) number, (size_t)sizeof(STORE_F),
	    (int(*)(const void *, const void *)) id_sort_nondiag);
      number = id_compress_nondiag_SD_elem(number, matr_ptr);
    }
    // terminate the block of matrix elements by an extra element

    curr_ptr = matr_ptr + number;
    curr_ptr->value     = D_ZERO;
    (curr_ptr++)->final = ULL_ZERO;
    number++;
    num_stored += number;
    matr_ptr += number;

  } // end loop through all |initSD>

  head->lastSD = loopSD - 1;
  if(num_stored) {
    head->numSD = num_stored;
    switch(head->typeSD) {

    case 0: id_binary_search_nondiag_asym_elem(head, sd_bas);
            break;
    case 1: id_binary_search_nondiag_sym_elem(head, sd_bas);
            break;
    case 2: id_binary_search_nondiag_nosym_elem(head, sd_bas);
            break;
    } // end switch()
  }
  free(occ); // release temporary memory

} // End: function twoPartNondiagSD_Calc()

    /*
    ** The entrance function
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
    ** The nondiag matrix elements are saved in memory pointed
    ** to by sd_store->permMemPtr
    */

void threePartNondiagSD_Calc(HEAD *head, SP_BAS *sp_bas,
				      SD_BAS *sd_bas, MATR_OP *op_veff)
{
  char      *func = {"i threePartNondiagSD_elem(): "};
  int       num_part, num_part_2,  num_orb, num1, num2, num3,
            num_stored, number,
            k, l, m, orb,*occ, k_phase, kl_phase, klm_phase, phase, count, 
            *k_ptr,*l_ptr,*m_ptr,
            loopSD, max_mem_elem, klm_elem_lim;
  ULL       *initSD, pos, sd_k, sd_kl, sd_klm, new_sd, low, high;
  STORE_F   *matr_ptr,*curr_ptr;
  ID_INT    **table,**tab_k,**tab_kl,*hij_ptr;

  max_mem_elem = head->memSize/sizeof(STORE_F); // storage 
  matr_ptr     = (STORE_F *)head->memPtr;

  klm_elem_lim  = (  sd_bas->part 
                   *(sd_bas->part -1) 
		   *(sd_bas->part -2))/6 * op_veff->maxNondiagElem;

     /*
      ** memory must exceed number of nondiag matrix 
      ** elements from one |initSD>
      */

  if(max_mem_elem < klm_elem_lim) {

    // return without any calculation

    head->numSD   = -1;
    return;
  }

  table       = op_veff->id_nondiag_table;  // <|H|>
  num_part    = sd_bas->part;
  num_part_2  = num_part - 2;

  num_orb = sd_bas->numm_orb;
  num1    = 3*num_orb*num_orb -12*num_orb + 11;
  num2    = 3*num_orb - 6;
  num3    = 2*num_orb - 3;

  num_stored = 0;

  occ = MALLOC(sd_bas->part, int, func, "occ[]");

  initSD = sd_bas->SD + head->startSD;

  for(loopSD = head->startSD; loopSD < head->endSD; loopSD++, initSD++) {

     /* 
     ** if not enough memory is available to save a new
     ** set of klm_elem_lim nondiag matrix elements save 
     ** current elements on file
     ** terminate calculation
     */ 

    if((max_mem_elem - num_stored) < klm_elem_lim) break;

    number   = 0; // initialization
    curr_ptr = matr_ptr;
    k_ptr    = occ;

    for(k = 0, orb = 0, pos = ULL_ONE; k < num_part; 
                             k++, orb++, pos <<= 1)  {
      for(;!((*initSD) & pos); orb++, pos <<= 1); 
      *(k_ptr++) = orb; 
    }
    k_ptr = occ;
    for(k = 0; k < num_part_2; k++, k_ptr++) { 
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
	  do  {      // run through all ij-pair
	    for(; sd_klm & hij_ptr->one; hij_ptr++);  
	    if(!hij_ptr->one) break; // no more ij-pairs
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
	    curr_ptr->value     = (float)(phase * hij_ptr->val);
	    (curr_ptr++)->final = new_sd;
	    number++;
	  } while((++hij_ptr)->one); // all VEFF for given (kl)
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
      number = id_compress_nondiag_SD_elem(number, matr_ptr);
    }
    // terminate the block of matrix elements by an extra element

    curr_ptr = matr_ptr + number;
    curr_ptr->value     = D_ZERO;
    (curr_ptr++)->final = ULL_ZERO;
    number++;
    num_stored += number;
    matr_ptr += number;

  } // end loop through all |initSD>

  head->lastSD = loopSD - 1;
  if(num_stored) {
    head->numSD = num_stored;
    switch(head->typeSD) {

    case 0: id_binary_search_nondiag_asym_elem(head, sd_bas);
            break;
    case 1: id_binary_search_nondiag_sym_elem(head, sd_bas);
            break;
    case 2: id_binary_search_nondiag_nosym_elem(head, sd_bas);
            break;
    } // end switch()
  }
  free(occ); // release temporary memory

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
} /* End: function id_sort_nondiag() */

      /*
      ** The function 
      **      id_compress_nondiag_SD_elem()
      ** takes num calculated nondiag SD elem which are sorted
      ** after increasing store[].final and add together all
      ** elements  with same store[].final.
      ** Note: The function must have num > 1;
      ** The final number of elements is returned.
      */

static int id_compress_nondiag_SD_elem(int num, STORE_F *store)
{
   int       k, l;
   STORE_F   *ptr1, *ptr2, *ptr3;

   ptr1  = store;                     /* initialization */
   ptr2  = ptr1 + 1;
   k = num - 1;
   while(--k) {
     if(ptr1->final == ptr2->final) {
       ptr1->value += ptr2->value;
       ptr3 = ptr2 + 1;
       l = k;
       while(--l) {
	 ptr2->final     = ptr3->final;
	 (ptr2++)->value = (ptr3++)->value;
       } /* end l-loop */
       ptr2->final = ptr3->final;
       ptr2->value = ptr3->value;
       ptr2 = ptr1 + 1;
       num--;
     }
     else {
       ptr1++;
       ptr2++;
     }
   } /* end k-loop */

   if(ptr1->final == ptr2->final) {     /* last pair of elements */
     ptr1->value += ptr2->value;
     num--;
   }
   return num;
   
} /* End: function id_compress_nondiag_SD_elem() */

    /*
    ** The function 
    **       id_binary_search_nondiag_asym_elem() 
    ** performs a binary search on all time-reversal
    ** asymmetric |final_SD()> and replace them with the
    ** corresponding sequence number.
    */

static void id_binary_search_nondiag_asym_elem(HEAD *head, 
                                     const SD_BAS *sd_bas)
{
  int       num, loopSD, new_num_elem, curr_initSD;
  ULL       top_bit, low, high, new_tr_sd, init_tr_sd, 
            search;
  double    j_phase;
  STORE_F   *matr_ptr1, *matr_ptr2;
  
  curr_initSD  = head->startSD;
  new_num_elem = head->numSD;
  top_bit      = ULL_ONE << (8 * sizeof(ULL) - 1);
  matr_ptr1    = (STORE_F *)head->memPtr;
  matr_ptr2    = matr_ptr1;
  j_phase      = sd_bas->J_type ? 1.0 : -1.0;

  for(loopSD = 0; loopSD < head->numSD; loopSD++, matr_ptr2++) {

       /*    
       ** change separation elem between groups
       ** of elements with same |init_SD>
       */

    if(matr_ptr2->final == ULL_ZERO) {
      curr_initSD++;
      matr_ptr1->final     = top_bit;
      (matr_ptr1++)->value = D_ZERO;
      continue;
    }
    if(fabs(matr_ptr2->value) < MATRIX_LIMIT) { // remove ZERO elements
      new_num_elem--;
      continue;  
    }
    num       = sd_bas->numm_orb; // calculate |new_tr_sd>
    new_tr_sd = ULL_ZERO;
    low       = ULL_ONE;
    high      = (ULL_ONE<<(num - 1));
    while(num--) {
      if(low & matr_ptr2->final) new_tr_sd |= high;
      low  <<= 1;
      high >>= 1;
    }
    if((matr_ptr2->final) < new_tr_sd)  { // |new_asymSD>
      low  = curr_initSD;
      high = sd_bas->numSD[0];
      while(1) {
        search = (low + high) >> 1;
	if(matr_ptr2->final < sd_bas->SD[search])      high = search - 1;
	else if(matr_ptr2->final > sd_bas->SD[search]) low  = search + 1;
	else                                           break;
      }
      matr_ptr1->value = matr_ptr2->value;
      (matr_ptr1++)->final = search << 2;
    } /* end new |asymSD> config  */

    else if(matr_ptr2->final == new_tr_sd){ // |new_symSD> 
      if(!sd_bas->J_type) { // Odd total J - no |new_symSD>
	new_num_elem--;
	continue;  
      }

      low  = sd_bas->numSD[0];
      high = sd_bas->tot_dimSD;
      while(1) {
	search = (low + high) >> 1;
	if(matr_ptr2->final < sd_bas->SD[search])      high = search - 1;
	else if(matr_ptr2->final > sd_bas->SD[search]) low  = search + 1;
	else                                           break;
      }
          /*
	  ** Add contribution from the excluded time-reversed 
	  ** |SD'>s when init_tr_sd  < matr_ptr2->final
	  */

      num        = sd_bas->numm_orb; // calculate |init_tr_sd>
      init_tr_sd = ULL_ZERO;
      low        = ULL_ONE;
      high       = (ULL_ONE<<(num - 1));
      while(num--) {
	if(low & sd_bas->SD[curr_initSD]) init_tr_sd |= high;
	low  <<= 1;
	high >>= 1;
      }    
      matr_ptr1->value = matr_ptr2->value;
      (matr_ptr1++)->final = (search << 2)
               + ((init_tr_sd  < matr_ptr2->final) ? 1 : 0);
    }     /* end new |sym_SD[]> config */

    else  {      /* |newSD> = |new_trSD> */
	
        /*
	** new_sd is a time reversed |SD> , binary search
	** in asym sd_bas->SD[]. Add contribution from the
	** excluded  time-reversed component.
	*/

      low  = 0;
      high = sd_bas->numSD[0];
      while(1) {
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
      matr_ptr1->value = j_phase * matr_ptr2->value;
      (matr_ptr1++)->final =  (search << 2)
            + ((init_tr_sd < new_tr_sd) ? 3 : 2);
    } // end new contribution to final tr_asym_SD[]
  
  }  // all |final_SD()> transferred to corresponding number
 
  head->numSD = new_num_elem; // update number of elements

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
  int      loop, new_num_elem, num, curr_initSD;
  ULL      top_bit, low, high, new_tr_sd, search;
  double   j_phase;
  STORE_F  *matr_ptr1, *matr_ptr2;

  curr_initSD  = head->startSD;  // initialization
  new_num_elem = head->numSD;
  top_bit      = ULL_ONE << (8 * sizeof(ULL) - 1);
  matr_ptr1    = (STORE_F *)head->memPtr;
  matr_ptr2    = matr_ptr1;
  j_phase      = sd_bas->J_type ? 1.0  : -1.0;
  for(loop = 0; loop < head->numSD; loop++, matr_ptr2++) {

       /*    
       ** change separation elem between groups
       ** of elements with same |init_SD>
       */

    if(matr_ptr2->final == ULL_ZERO) {
      curr_initSD++;
      matr_ptr1->final     = top_bit;
      (matr_ptr1++)->value = D_ZERO;
      continue;
    }
    if(fabs(matr_ptr2->value) < MATRIX_LIMIT)  {
      new_num_elem--;
      continue;  
    }
    num       = sd_bas->numm_orb; // calculate |new_tr_sd>
    new_tr_sd = ULL_ZERO;
    low       = ULL_ONE;
    high      = (ULL_ONE<<(num - 1));
    while(num--) {
      if(low & matr_ptr2->final) new_tr_sd |= high;
      low  <<= 1;
      high >>= 1;
    }
    if((matr_ptr2->final) < new_tr_sd)  {// |new_asymSD>
      low  = 0;
      high = sd_bas->numSD[0];
      while(1) {
        search = (low + high) >> 1;
	if(matr_ptr2->final < sd_bas->SD[search])      high = search - 1;
	else if(matr_ptr2->final > sd_bas->SD[search]) low  = search + 1;
	else                                          break;
      }
      matr_ptr1->value = matr_ptr2->value;
      (matr_ptr1++)->final = search << 2;
    } // end new |asymSD> config

    else if(matr_ptr2->final == new_tr_sd){ // |new_symSD>
      low  = curr_initSD;
      high = sd_bas->tot_dimSD;
      while(1) {
	search = (low + high) >> 1;
	if(matr_ptr2->final < sd_bas->SD[search])      high = search - 1;
	else if(matr_ptr2->final > sd_bas->SD[search]) low  = search + 1;
	else                                           break;
      }
      matr_ptr1->value = matr_ptr2->value;
      (matr_ptr1++)->final = search << 2;
    }     /* end new |sym_SD[]> config */

    else  {  // |newSD> = |new_trSD>

        /*
	** |new_sd> is a tr_asym_SD[]>, search in 
	** final state sd_bas->SD[0]. Add contribution 
	** from the excluded time-reversed component.
	*/
      low  = 0;
      high = sd_bas->numSD[0];
      while(1) {
        search = (low + high) >> 1;
	if(new_tr_sd < sd_bas->SD[search])      high = search - 1;
	else if(new_tr_sd > sd_bas->SD[search]) low  = search + 1;
	else                                          break;
      }
      matr_ptr1->value = matr_ptr2->value;
      (matr_ptr1++)->final =  (search << 2) + 1;
    } // end new contribution to final tr_asym_SD[]
  
  }  // all |final_SD()> transferred to corresponding number
 
  head->numSD = new_num_elem;  // update number of elements

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
  int
              loop, new_num_elem, curr_initSD;
  ULL
              top_bit, low, high, search;
  STORE_F
              *matr_ptr1, *matr_ptr2;

  curr_initSD  = head->startSD; // initialization
  new_num_elem = head->numSD;
  top_bit      = ULL_ONE << (8 * sizeof(ULL) - 1);
  matr_ptr1    = (STORE_F *)head->memPtr;
  matr_ptr2    = matr_ptr1;
  search       = curr_initSD;

  for(loop = 0; loop < head->numSD; loop++, matr_ptr2++) {

       /*    
       ** change separation elem between groups
       ** of elements with same |init_SD>
       */

    if(matr_ptr2->final ==  ULL_ZERO) {
      curr_initSD++;
      search               = curr_initSD;
      matr_ptr1->final     = top_bit;
      (matr_ptr1++)->value = D_ZERO;
      continue;
    }
    if(fabs(matr_ptr2->value) < MATRIX_LIMIT) { // remove ZERO elements
     new_num_elem--;
      continue;  
    }
    low  = search + 1;
    high = sd_bas->numSD[0];
    while(1) {
      search = (low + high) >> 1;
      if(matr_ptr2->final < sd_bas->SD[search])      high = search - 1;
      else if(matr_ptr2->final > sd_bas->SD[search]) low  = search + 1;
      else                                          break;
    }
    matr_ptr1->value = matr_ptr2->value;
    (matr_ptr1++)->final = search;
  } // all |final_SD()> transferred to corresponding number

  head->numSD = new_num_elem; // update number of elements

} // End: function id_binary_search_nondiag_nosym_elem()
