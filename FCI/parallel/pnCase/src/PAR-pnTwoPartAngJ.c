
/*******************  The module pn_two_part_angJ.c  ******************/

#include  "PAR-pnShellModel.h"

     /*
     ** The entrance function 
     **         pn_two_part_angJ()
     ** reads from file angular momentum coupled two-particle
     ** identical particle matrix elements, transfer the elements
     ** into m-scheme representation and store diagonal matrix
     ** elements in 
     **         MATR_OP op_int->id_diag[] 
     ** and pointers to the nondiag matrix elements in 
     **         MATR_OP op_int->id_nondiag_table[]
     */


          /******  local function declarations  ********/

static void pn_m_pot_diag(int num_mZ, MBAS *mbasZ,
	       int num_mN, MBAS *mbasN, double *diag_ang);

     /*
     ** calculates all diagonal m-scheme proton-neutron
     ** two-particle matrix elements
     **       < p(k), n(l) |J**2 | p(k), l(n) >
     ** The matrix elements are  stored  in  
     **          diag_ang[k * num_mN + l].diag_val
     */


static int pn_max_nondiag_elem(int parZ, SP_BAS *spBasZ, 
                                         SP_BAS *spBasN);
     /*
     ** calculates and returns the maximum  possible number of 
     ** proton-neutron  non-diagonal two-body matrix elements.
     */

static int pn_m_pot_nondiag(int parZ, SP_BAS *spBasZ, SP_BAS *spBasN,
			    MATR_PN *op_int);
     /*
     ** calculates all nondiagonal m-scheme proton-neutron two-particle
     ** matrix elements < p(i), n(j) | J**2 | p(k), l(n)>
     ** The nondiagonal matrix elements  are  stored in a list
     ** pn_nondiag_elem[] pointed to by  
     ** table[((k * num_mN+ l].start_nondiag
     ** The function returns number of calculated matrix elements.
     */


static double  pn_m_scheme_nondiag_ang(int m_orb_a, int m_orb_b,
	       int m_orb_c, int m_orb_d, MBAS *mbasZ, MBAS *mbasN);
     /*
     ** calculates and return  non-diagonal two-particle matrix
     ** element in m-scheme of angular momentum  operator                                    
     **             <i,j| J**2 |k,l>                
     ** for proton-neutron case.
     */

static void  addPN_1to2SinglePartTerms(int type, int part_num, 
	  int num_mZ, int num_mN, MBAS *mbas, MATR_PN *op_ang);
     /*
     ** adds single-particle contribution for protons(neutrons) when type = 0(1)
     ** to two-particle matrix elements  <k(p).l(n) |J**2 | k(p).l(n)>
     */

          /******  end local function declarations *****/


     /*
     ** The entrance function 
     **         pn_two_part_angJ()
     ** calculates m-scheme two-particle matrix
     ** elements of J**2 and store the result in 
     ** MATR_PN op_ang->pn_diag[] and 
     ** MATR_PN op_ang->pn_nondiag_table[]
     */

void pn_two_part_angJ(int parZ, SP_BAS *spBasZ,SP_BAS *spBasN, 
                                               MATR_PN *op_ang)
{
  char     *func = {" pn_two_part_int(): "};
  int      num_kl_elem, num_kl_tables = 0, k, numCalc;

  // if possible store data in MATR_ID *op_int for VEFF

  if(op_ang->numNondiagElem == 0) {

    // num del_m and del_p kl-tables

    num_kl_tables = (spBasZ->mbas[0].m + 1) * (parZ ? 1 : 2);

    op_ang->kl_list = MALLOC(num_kl_tables, PN_INT **, func,
			               "op_ang->kl_list[]");

    num_kl_elem = spBasZ->numm_orb * spBasN->numm_orb;

    for(k = 0; k < num_kl_tables; k++) {
      op_ang->kl_list[k]
	= MALLOC(num_kl_elem, PN_INT *, func,"pn_ang_kl_table[]");
    }
    op_ang->numNondiagElem  = pn_max_nondiag_elem(parZ, spBasZ, spBasN);

    op_ang->nondiag_bas = MALLOC(op_ang->numNondiagElem, PN_INT, func,
				              "op_ang->nondiag_bas[]");
  }
    // calculate and store diag m-scheme two-particle matr elem

  pn_m_pot_diag(spBasZ->numm_orb,spBasZ->mbas,spBasN->numm_orb,
                                    spBasN->mbas,op_ang->pn_diag);

  // calculate and store nondiag m-scheme two-particle matr elem

  numCalc = pn_m_pot_nondiag(parZ, spBasZ, spBasN,op_ang);

  if(numCalc > op_ang->numNondiagElem) {
    printf("\n\nRank%d: Error in pn_two_part_angJ(): ",Rank);
    printf("\n pn_matrElem: mem reserve = %d - calc = %d",
	                   op_ang->numNondiagElem, numCalc);
    printf("\nThis may not occur!!!!!\n");
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }
  if(numCalc == 0) {
    free(op_ang->nondiag_bas);
    op_ang->nondiag_bas      = (PN_INT *)NULL_PTR;
    op_ang->numNondiagElem= 0;
    for(k = num_kl_tables - 1; k >= 0; k--) {
      free(op_ang->kl_list[k]);
    }
    free(op_ang->kl_list);
  }
  else if(numCalc < op_ang->numNondiagElem) { // reduce size of nondiag[] to new number

/***********************************
    op_ang->numNondiagElem = numCalc;
    nondiag_bas = op_ang->nondiag_bas;
    op_ang->nondiag_bas = REALLOC(nondiag_bas,op_ang->numNondiagElem,
                                    PN_INT, func,"new nondiag_bas[]");
    if(op_ang->nondiag_bas != nondiag_bas) { 
      printf("\n\nRank%d: Error in pn_two_part_angJ(): ");
      printf("\nREALLOC() does not work properly - move mem loc");
      printf("\nof nondiag matrElem to a new place\n");
      MPI_Abort(MPI_COMM_WORLD,Rank);
    }
************************************/



  }
       /*
       ** If number of proton and/or neutron is less
       ** than two the single-particle contributions
       ** to effective interation are added here 
       */
  
  if(spBasZ->part == 1) {
    addPN_1to2SinglePartTerms(PROTON, spBasN->part, spBasZ->numm_orb,
			     spBasN->numm_orb, spBasZ->mbas, op_ang);
  }
  if(spBasN->part == 1) { 
    addPN_1to2SinglePartTerms(NEUTRON, spBasZ->part, spBasZ->numm_orb,
                             spBasN->numm_orb,spBasN->mbas, op_ang);
  }
  return;

} // End: function pn_two_part_angJ()

     /*
     ** The function                                          
     **              pn_m_pot_diag(...)                    
     ** calculates all diagonal m-scheme proton-neutron
     ** two-particle matrix elements
     **       < p(k), n(l) |J**2 | p(k), l(n) >
     ** The matrix elements are  stored  in  
     **          diag_ang[k * num_mN + l].diag_val
     */

static void pn_m_pot_diag(int num_mZ, MBAS *mbasZ,
                          int num_mN, MBAS *mbasN, double *diag_ang)
{
  int   k, l;

  for(k = 0; k < num_mZ; k++) { // proton right
    for(l = 0; l < num_mN; l++, diag_ang++) {  // neutron right
      *diag_ang = 0.5 * ((double) (mbasZ[k].m * mbasN[l].m));
    } // end index l, neutron right
  } // end index k, proton right
   
} // End: function pn_pn_m_pot_diag()

     /*
     ** The function
     **        pn_max_nondiag_elem()
     ** calculates and returns the maximum  possible number of 
     ** proton-neutron  non-diagonal two-body matrix elements.
     */

static int pn_max_nondiag_elem(int parZ, SP_BAS *spBasZ, SP_BAS *spBasN)
{
  int   total_num, i, j, k, l, m_k, m_l,par_k, par_l,
        max_del_mp, del_m, del_par, par_p, limit, count, number;
    
  // Calculate num of non-diag matr el<em with i(p) >= k(p) and j(n) >= l(n)

  max_del_mp = spBasZ->mbas[0].m;

  limit = parZ ? 1 : 2; // number of parity change in matr.elem

  total_num = 0;  // initialization

  for(k = 0; k < spBasZ->numm_orb; k++) { // right proton particle
    m_k   = spBasZ->mbas[k].m;
    par_k = spBasZ->mbas[k].par; 
    for(l = 0; l < spBasN->numm_orb; l++) { // right neutron particle
      m_l   = spBasN->mbas[l].m;
      par_l = spBasN->mbas[l].par;
      
      number = 0; // counter for each (kl) - group
      
      for(del_m = 0; del_m <= max_del_mp; del_m++)  {
         
	for(del_par = 0,par_p = +1;del_par < limit; del_par++, par_p = -1) {

	  count = 0; // element counter for each block with given del_par
         
	  total_num++; // header element
            
	  for(i = ((del_m == 0) && (del_par == 1)) ? 0 : k; 
                                    i < spBasZ->numm_orb; i++)   {
	    for( ; !(  (del_m == ((m_k - spBasZ->mbas[i].m) / 2))
                     &&(par_p == (par_k * spBasZ->mbas[i].par)))      
                     && ( i < spBasZ->numm_orb); i++);
	    if(i ==  spBasZ->numm_orb) break;
	    for(j = (i == k) ? l + 1: 0; j < spBasN->numm_orb; j++)   {
	      for( ;   !(  (del_m == ((spBasN->mbas[j].m - m_l) / 2)) 
		     &&(par_p == (spBasN->mbas[j].par * par_l)))
                     &&(j < spBasN->numm_orb); j++);
	      if(j == spBasN->numm_orb) break; 
                  
	      count++;
	      number++;                 
	      total_num++; // matrix element
               
	    } // end loop through j-neutrons
	  } // end loop through i-protons
	  if(!count)  {
	    total_num--;// no matrix element in the present group
	  }
	} // end del_par loop

      } // end del _m for the two parity blocks of matrix elements
         
      if(number) total_num++; // extra bottom element

    } // end of index l, neutron particle right
  } // end of  index k, proton particle right

  return total_num;

} // End: function pn_max_nondiag_elem()

     /*
     ** The function                                          
     **              pn_m_pot_nondiag(,...)                    
     ** calculates all nondiagonal m-scheme proton-neutron two-particle
     ** matrix elements < p(i), n(j) | J**2 | p(k), l(n)>
     ** The nondiagonal matrix elements  are  stored in a list
     ** pn_nondiag_elem[] pointed to by  
     ** table[((k * num_mN+ l].start_nondiag
     ** The function returns number of calculated matrix elements.
     */

static int pn_m_pot_nondiag(int parZ, SP_BAS *spBasZ, SP_BAS *spBasN,
			                             MATR_PN *op_ang)
{
  int       num_mZ, num_mN, i, j, k, l, m_k, m_l,par_k, par_l, 
            del_m, del_par, num_m_par, num_parity, par_p,
            number, tot_num;
  MBAS      *mbasZ, *mbasN;
  PN_INT    **kl_table, *nondiag;

  num_mZ     = spBasZ->numm_orb;
  num_mN     = spBasN->numm_orb;
  num_parity = parZ ? 1 : 2;  // number of parity change in matr.elem
  mbasZ      = spBasZ->mbas;  
  mbasN      = spBasN->mbas;

  nondiag   = op_ang->nondiag_bas;

  tot_num = 0;

  num_m_par = 0;
  for(del_m = 0; del_m <=  spBasZ->mbas[0].m; del_m++){
 
   for(del_par = 0, par_p = +1; del_par < num_parity; 
                       del_par++, par_p = -1) {

       /* 
       ** for each del_m and del_p: calculate non-diagonal
       **  matrix elements with i(p) >= k(p) and j(n) >= l(n),
       */
    
     kl_table =  op_ang->kl_list[num_m_par];
     num_m_par++;
     
      for(k = 0; k < num_mZ; k++) { // right proton particle
	m_k   = mbasZ[k].m;
	par_k = mbasZ[k].par;
	for(l = 0; l < num_mN; l++, kl_table++) {  // right neutron particle
	  m_l   = mbasN[l].m;
	  par_l = mbasN[l].par;

	  *kl_table = nondiag;
	  number = 0;              /* counter for each (kl) - group */


	  for(i = ((del_m == 0) && (del_par == 1)) ? 0 : k; i < num_mZ; i++)   {
	    for( ;  !(  (del_m == ((m_k   - mbasZ[i].m) / 2))
                       &&(par_p == (par_k * mbasZ[i].par)))      
                       && ( i < num_mZ); i++);
	    if(i == num_mZ) break;
	    for(j = (i == k) ? l + 1: 0; j < num_mN; j++)   {
	      for( ;   !(  (del_m == ((mbasN[j].m - m_l) / 2)) 
                          &&(par_p == (mbasN[j].par * par_l)))
                          && ( j < num_mN); j++);
	      if(j == num_mN) break; 
              nondiag->val[0] = pn_m_scheme_nondiag_ang(i,j,k,l,mbasZ,mbasN);

	      if(fabs(nondiag->val[0]) > MATRIX_LIMIT) {
		nondiag->one[0] =  ULL_ONE<<i;
		nondiag->one[1] =  ULL_ONE<<j;
		nondiag->two[0]
                            = (i == k) ? 0 : (ULL_ONE<<MAX(i,k)) - (ULL_ONE<<(MIN(i,k) + 1));
		nondiag->two[1] 
                              = (j == l) ? 0 : (ULL_ONE<<MAX(j,l)) - (ULL_ONE<<(MIN(j,l)+1));
		number++;
		tot_num++;
		nondiag++;
	      }
	    }  // end loop through j-neutrons
	  } // end loop through i-protons
	  if(number)   { 
	    nondiag->one[0]     = 0;  // values to terminate a block of matr elem
	    (nondiag++)->one[1] = 0;
	    tot_num++;
	  }     
	  else   {  // no matrix elements in the present block
	    *kl_table = (PN_INT *)NULL_PTR;
	  }
	} // end of index l, neutron particle right
      } // end of  index k, proton particle right

    } // end del_par loop
  } // end del _m for the two parity blocks of matr elem

  return tot_num; // return total number of matr elem

} // End: function pn_m_pot_nondiag()

     /*
     ** The function                                
     **       pn_m_scheme_nondiag_ang()            
     ** calculates and return  non-diagonal two-particle matrix
     ** element in m-scheme of angular momentum  operator                                    
     **             <i,j| J**2 |k,l>                
     ** for proton-neutron case.
     */

static double  pn_m_scheme_nondiag_ang(int m_orb_a, int m_orb_b,
                        int m_orb_c, int m_orb_d, MBAS *mbasZ, MBAS *mbasN)
{
   int       j_a, m_a, j_b, m_b, j_c, m_c, j_d, m_d, value;

   if(  (mbasZ[m_orb_a].orb != mbasZ[m_orb_c].orb) 
     || (mbasN[m_orb_b].orb != mbasN[m_orb_d].orb)) {
     return D_ZERO; // J**2 matrix element is ZERO
   }

   j_a = mbasZ[m_orb_a].j; // m_orb_a and m_orb_c are particle one
   m_a = mbasZ[m_orb_a].m;
   j_c = mbasZ[m_orb_c].j;
   m_c = mbasZ[m_orb_c].m;

   j_b = mbasN[m_orb_b].j; // m_orb_b and m_orb_d are particle one
   m_b = mbasN[m_orb_b].m;
   j_d = mbasN[m_orb_d].j;
   m_d = mbasN[m_orb_d].m;

   value = 0;

   if((m_a == (m_c - 2)) && (m_b == (m_d + 2)))  { // m_a=m_c-1 m_b=m_d+1
      value = (j_c + m_c) * (j_c - m_c + 2) * (j_d - m_d) * (j_d + m_d + 2);
   }
   else if((m_a == (m_c + 2)) && (m_b == (m_d - 2))) { // m_a=m_c+1 m_b=m_d-1
      value = (j_c - m_c) * (j_c + m_c + 2) * (j_d + m_d) * (j_d - m_d + 2);
   }
   if(value)   {
      return (  0.25 * sqrt((double) value)
		* mbasZ[m_orb_a].phase * mbasZ[m_orb_c].phase // Whitehead time-reversal phase
              * mbasN[m_orb_b].phase * mbasN[m_orb_d].phase);
   }  
   else       return D_ZERO;

} // End: function pn_pn_m_scheme_nondiag_ang_mom()

     /*
     ** The function
     **                  addPN_1to2SinglePartTerms()
     ** adds single-particle contribution for protons(neutrons) when type = 0(1)
     ** to two-particle matrix elements  <k(p).l(n) |J**2 | k(p).l(n)>
     */

static void  addPN_1to2SinglePartTerms(int type, int part_num, 
              int num_mZ, int num_mN, MBAS *mbas, MATR_PN *op_ang)
{
   int       k, l;
   double    *pn_diag, ang_mom, num;

   pn_diag = op_ang->pn_diag;

   num = 1.0/((double)(part_num));
   for(k = 0; k < num_mZ; k++) {
      if(type == 0)  {
         ang_mom = (0.25 * mbas[k].j * (mbas[k].j + 2.0))* num;
      }
      for(l = 0; l < num_mN; l++, pn_diag++) {
         if(type == 1)  {
            ang_mom = (0.25 * mbas[l].j * (mbas[l].j + 2.0)) * num;
         }
	 *pn_diag += ang_mom;
      } // end l-loop
   } // end k-loop

} // End: function   addPN_1to2SinglePartTerms()
