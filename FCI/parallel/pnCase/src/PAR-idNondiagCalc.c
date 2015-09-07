
//**************** File: idNondiag_calc.c *************/

#include "PAR-pnShellModel.h"

    /*
    ** The entrance function
    **           idTwoPartNondiagSD_Calc()
    ** calculates nondiagonal identical particle matrix elements 
    **             <SD()'|OP(typeCM)|SD()>
    ** typeCM = 0 - VEFF_INT and ANG_INT: = 1 - CN_INT
    ** from |SD()> = |head->startSD> and store as many as possible
    ** in head->memPtr.
    ** |SD()> are many-particle basis state.
    ** Pointers to all two-particle matrix elements are stored in 
    ** structure id_diag by the formula
    **        op(k,l) = ((num - k) * k)/2 -1 +l
    ** where num = 2 * m_orb - 3.
    ** The function returns endOfCalc = YES (if all elements 
    ** are calculated), otherwise NO
    */

/***************   local function declarations*********/

static int id_sort_nondiag(const STORE_F *one, const STORE_F *two);
      /*
      **        int id_sort_nondiag()                  
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

static void id_binary_search_nondiag_elem(int type, 
                       GROUP *group, int numStored,
                         STORE_F *tempMemNondiagSD);
      /*
      ** performs a binary search on all nosymmetric |final_SD()>
      ** and replace them with the corresponding sequence number.
      ** Then all calculated  <SD()_final|OP|SD()_init> are 
      ** stored on file.
      */

static int occLim(int type,GR_BAS *grBas,ULL new_sd);
    /* 
    ** calculates sumN for |new_SD>.
    ** Return FALSE if sumN > N_excitation
    ** otherwise TRUE
    */
/**************  end local function declarations ********/

    /*
    ** The entrance function
    **           idTwoPartNondiagSD_Calc()
    ** calculates nondiagonal identical particle matrix elements 
    **             <SD()'|OP(typeCM)|SD()>
    ** typeCM = 0 - VEFF_INT and ANG_INT: = 1 - CN_INT
    ** from |SD()> = |head->startSD> and store as many as possible
    ** in head->memPtr.
    ** |SD()> are many-particle basis state.
    ** Pointers to all two-particle matrix elements are stored in 
    ** structure id_diag by the formula
    **        op(k,l) = ((num - k) * k)/2 -1 +l
    ** where num = 2 * m_orb - 3.
    ** The function returns endOfCalc = YES (if all elements 
    ** are calculated), otherwise NO
    */

int idTwoPartNondiagSD_Calc(int type, int typeCM,GR_BAS *grBas,
                             GROUP *group,STORE_F *idMemBlockPtr)
{
  int        num_part, num_part_1, num, numStored,loopSD_start,numSD,
             k, l,*occ,kl_phase, phase,count,*k_ptr, *l_ptr, loopSD,
             loopStep,occStep, number;
  ULL        *initSD, sd_k, sd_kl, new_sd, low, high;
  STORE_F    *currPtr, *startLoopPtr;
  ID_INT     **table, **tab_k, *ij_ptr;

  table         = (type == PROTON) ? grBas->op_pp.id_nondiag_table
                                  : grBas->op_nn.id_nondiag_table;
  num_part      = grBas->spBas[type]->part;
  num_part_1    = num_part - 1;
  num           = (grBas->spBas[type]->numm_orb << 1) - 3;

  occ = group->occ[type]; 
  
  if(type == PROTON) {
    loopSD_start  = 0;
    numSD         = group->numSD[PROTON] - 1;
    initSD        = group->SD[type];
    loopStep      = 1;
  }
  else {  // NEUTRON
    loopSD_start  = Rank;
    numSD         = group->numSD[NEUTRON] - 1;
    initSD        = group->SD[type] + Rank;
    loopStep      = NumProcs;
  }
  numStored     = 0;  // initalization
  currPtr       = idMemBlockPtr;

  for(loopSD = loopSD_start, occStep = 0; loopSD < numSD; 
             loopSD += loopStep,occStep++,initSD += loopStep) {

    number       = 0; // num elem in current loopSD
    startLoopPtr = currPtr;

    k_ptr    = &occ[occStep*(num_part + 1)];

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

	  // check oscillator exitation limit M_excitation

	  if(occLim(type,grBas,new_sd) == FALSE) continue;

             /* 
	     ** check that |new_sd> has allowed number of 
	     ** particles in all j-orbits - only active if
	     **  0 < allowed number < 2*j + 1
             */
	  
	  if(grBas->spBas[type]->mask.num > 0) {    /* check masking */
	    for(count = 0; count < grBas->spBas[type]->mask.num; count++) {
	      high = new_sd & grBas->spBas[type]->mask.list[count];
	      phase = 0;
	      while(high) {
		high &= high - 1;
		phase++;
	      }
	     if(  (phase > grBas->spBas[type]->mask.lim[2*count + 1])
		  ||(phase < grBas->spBas[type]->mask.lim[2*count])) {
               break; // |new_sd> not accepted
	     }
	    }
            if(count < grBas->spBas[type]->mask.num) continue; 

	  } //end check masking

	  low   = new_sd & ij_ptr->two; // permutation phase
	  phase = kl_phase;
	  while(low) {
	    low &= low - 1;
	    phase = -phase;
	  }

	  {   // Has |new_sd> been calc allready ?
	    int kk;
	    for(kk= 0; kk < number;kk++) {
	      if(new_sd == startLoopPtr[kk].final) {
		startLoopPtr[kk].value 
                    += (float)(phase * ij_ptr->val[typeCM]);
		break;  // return to new two-part ij_ptr matr elem
	      }
	    }
	    if(kk == number) {// |new_sd>  is NEW
	      currPtr->value = (float)(phase * ij_ptr->val[typeCM]);
	      currPtr->final = new_sd;
	      currPtr++;
	      number++;
	    }
	  }

	} // end all two-body matr elem for given (kl)  
      }   // end l-loop
    }     // end k-loop

         /* 
	 ** sort and compress nondiag matrix elements for 
	 ** each |init_SD> after increasing |final_SD>
	 */
  
    if(number > 1) {
      qsort(startLoopPtr, (size_t) number, (size_t)sizeof(STORE_F),
	    (int(*)(const void *, const void *)) id_sort_nondiag);
    }
    // terminate |initSD> block of matr elem by an extra element

    currPtr        = startLoopPtr + number;
    currPtr->value = D_ZERO;
    currPtr->final = ULL_ZERO;
    currPtr++;
    number++;

    numStored += number;

  } // end loop through all |initSD>

  id_binary_search_nondiag_elem(type,group,numStored,idMemBlockPtr);

  return  numStored;
    
} // End: function idTwoPartNondiagSD_Calc()

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
      ** after increasing store[].final and add together all
      ** elements  with same store[].final
      ** Then all element store{].value < MATRIX_LIMIT
      ** are removed.  .
      ** Note: The function must have num > 1;
      ** The final number of elements is returned.
      */

static int id_compress_nondiag_SD_elem(int num, STORE_F *store)
{
  int  lim, k, l;
   
   lim = num - 1;
   k   = 0;
   do {
     for( ; k < lim; k++) {
       if(store[k].final == store[k + 1].final) {
	 store[k].value += store[k + 1].value;
	 for(l = k + 1; l < lim; l++) {
	   store[l].value = store[l + 1].value;
	   store[l].final = store[l + 1].final;
	 } // end loop l
	 lim--;
	 break;
       }  //end if-test
     }  // end loop k
   } while(k < lim);

   lim++;
   k = 0;
   do {
     for( ; k < lim; k++) {
       if(fabs(store[k].value) > MATRIX_LIMIT) continue;
       store[k].value = store[k + 1].value;
       store[k].final = store[k + 1].final;
       for(l = k + 1; l < lim; l++) {
	 store[l].value = store[l + 1].value;
	 store[l].final = store[l + 1].final;
       } // end loop l
       lim--;
       break;
     }  // end loop k
   } while(k < lim);

   return lim;

} // End: function id_compress_nondiag_SD_elem()

    /*
    ** The function 
    **       id_binary_search_nondiag_elem() 
    ** performs a binary search on all nosymmetric |final_SD()>
    ** and replace them with the corresponding sequence number.
    ** Then all calculated  <SD()_final|OP|SD()_init> are 
    ** stored on file.
    */

static void id_binary_search_nondiag_elem(int type, 
                       GROUP *group, int numStored,
                         STORE_F *tempMemNondiagSD)
{
  int     loopSD, numSD;
  ULL     top_bit, low, high, search;
  STORE_F *sdPtr;

  numSD   = (type == PROTON) ? group->numSD[0] 
                             : group->numSD[1];
  top_bit = ULL_ONE << (8 * sizeof(ULL) - 1);
  sdPtr   = tempMemNondiagSD;

  for(loopSD = 0;loopSD <  numStored;loopSD++) {

       /*    
       ** change separation elem between groups
       ** of elements with same |init_SD>
       */

    if(sdPtr[loopSD].final ==  ULL_ZERO) {
      sdPtr[loopSD].final = top_bit;
      sdPtr[loopSD].value = D_ZERO;
      continue;
    }

    low  = 0;    // low  = search + 1;
    high = numSD;
    for( ; ;) {            // while(1) {
      search = (low + high) >> 1;
      if(sdPtr[loopSD].final < group->SD[type][search])      high = search - 1;
      else if(sdPtr[loopSD].final > group->SD[type][search]) low  = search + 1;
      else                                                break;
    }
    sdPtr[loopSD].final = search;
  } // all |final_SD()> transferred to corresponding number

}  // End: function id_binary_search_nondiag_nosym_elem()

    /* 
    ** The function
    **       occLim()
    ** calculates sumN for |new_SD>.
    ** Return FALSE if sumN > N_excitation
    ** otherwise TRUE
    */

static int occLim(int type,GR_BAS *grBas,ULL new_sd)
  {
    int    num_part,par, orb, sumN, maxOscN;
    ULL   pos;
    MBAS  *mBasPtr;

    if(type == PROTON) {
      maxOscN = grBas->maxOscN_Z;
    }
    else {
      maxOscN = grBas->maxOscN_N;;
    }

    num_part     = grBas->spBas[type]->part; // initialization
    mBasPtr      = grBas->spBas[type]->mbas;

    sumN     = 0;

    for(par = 0, orb = 0, pos = ULL_ONE; par < num_part;
                                 par++, orb++, pos <<= 1)  {
      for(;!(new_sd & pos); orb++, pos <<= 1);
      sumN += mBasPtr[orb].N; 
      
    } // end loop through all |SD> in sd[]
    
    if(sumN > maxOscN) {
      return FALSE;
    }
    else  {
      return TRUE;
    }
  
  } // End: function ooccLim()
