
/******************  The module PAR-id-two-part-angJ.c  ******************/

#include "PAR-pnShellModel.h"

     /*
     ** The module entrance function 
     **    id_two_part_angJ()
     ** calculates two-particle matrix elements of J**2 
     ** in m-scheme representation and store the result 
     ** in tables:ID_TABLE id_table[]
     */


       /****   local function declarations **********/

static void id_m_pot_diag(SP_BAS *spBas, MATR_ID *op_ang);
     /*
     ** The function                                          
     **              id_m_pot_diag(...)                    
     ** calculates all diagonal m-scheme J**2 matrix elements
     ** <k.l| J**2 |k.l>  for identical particles 
     ** The diagonal matrix elements  for k < l are  stored  in vector: 
     **         table[((2*num_m - k - 3) * k)/2 + l].diag_val
     */

static void id_m_pot_diag(SP_BAS *spBas, MATR_ID *op_ang);
     /*
     ** calculates all diagonal m-scheme J**2 matrix elements
     ** <k.l| J**2 |k.l>  for identical particles 
     ** The diagonal matrix elements  for k < l are  stored  in vector: 
     **         table[((2*num_m - k - 3) * k)/2 + l].diag_val
     */

static void id_m_pot_diag(SP_BAS *spBas, MATR_ID *op_ang);
     /*
     ** calculates all diagonal m-scheme J**2 matrix elements
     ** <k.l| J**2 |k.l>  for identical particles 
     ** The diagonal matrix elements  for k < l are  stored  in vector: 
     **         table[((2*num_m - k - 3) * k)/2 + l].diag_val
     */

static double  id_m_scheme_diag_ang(int k, int l, MBAS *mbas);
     /*
     ** calculates the diagonal two-particle matrix element 
     ** in m-scheme of 
     **      <k,l| J^2 |k,l> for for identical particles.
     */

static void addID_1to2SinglePartTerms(SP_BAS *spBas,MATR_ID *op_ang);
     /*
     ** adds contributions from single-particle terms of J**2
     ** to the two-particle matrix elements  <k.l |J**2| k.l>
     */

static int id_m_pot_nondiag(SP_BAS *spBas,MATR_ID *op_ang);
     /*
     ** calculates all nondiagonal m-scheme two-particle matrix elements
     ** <i,j| J**2 |k.l>  for identical particles
     ** The nondiagonal matrix elements  are  stored in a table pointed to by  
     ** table[((2*num_m - k - 3) * k)/2 + (l - 1)]
     ** The function returns number of calculated matrix elements.
     */

static double id_m_scheme_nondiag_ang(int i, int j, int k, int l,
					                     MBAS *mbas);
     /*
     ** calculates and returns non-diagonal two-particle 
     ** matrix element in m-scheme 
     **             <i,j| J**2 |k,l>
     ** for the identical particle case.
     ** It is assumed that i < j and k < l.         
     */

                /**** End: function declarations ****/


               /**** The function definitions  ****/ 

     /*
     ** The module entrance function 
     **    id_two_part_angJ()
     ** calculates two-particle matrix elements of J**2 
     ** in m-scheme representation and store the result 
     ** in tables:ID_TABLE id_table[]
     */

void id_two_part_angJ(SP_BAS *spBas,MATR_ID *op_ang)
{
  int      num;

  // if possible store data in MATR_ID *op_int for VEFF

  //  num_elem = op_ang->num_nondiag_elem;
  //if(num_elem == 0) {


 /******************  Fjernes ********************  
    num_elem = (spBas->numm_orb * (spBas->numm_orb - 1)) >> 1;
    op_ang->id_nondiag_table
         = MALLOC(num_elem, ID_INT *, func,"id_kl_ang_nondiag_table[]");
    num_elem            = max_nondiag_id_elem(spBas->numm_orb, spBas->mbas);
    id_nondiag_elem_bas = MALLOC(num_elem, ID_INT, func,"id_nondiag_ang_elements[]");
    op_ang->id_nondiag_bas = id_nondiag_elem_bas;

    //}
 *******************  Til HIT  **********/

  // calculate and store diag m-scheme two-particle matrix elements

  id_m_pot_diag(spBas,op_ang);

        /* 
	** add single-particle J**2 energies to the matrix
	** elements for the current particle number
	*/

  addID_1to2SinglePartTerms(spBas,op_ang);

  // calculate and store nondiag m-scheme J**2 matrix elements

  num = id_m_pot_nondiag(spBas,op_ang); 
  op_ang->num_nondiag_elem = num;  //reserved for later use



 /****************   Fjernes  *****************

  if(num > num_elem) {
    printf("\n\nRank%d: Error in function id_two_part_angJ(): ",Rank);
    printf(" = %d - calculated = %d",num_elem, num); 
    printf("\n This may not occur!!!!!\n");
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }
  if(num == 0) {
    free( op_ang->id_nondiag_bas);
    op_ang->id_nondiag_bas   = NULL_PTR;
    free(op_ang->id_nondiag_table);
    op_ang->id_nondiag_table = NULL_PTR;
    op_ang->num_nondiag_elem = 0;
  }
  else if(num < num_elem) {  // reduce size of nondiag[] to new number
    op_ang->id_nondiag_bas = REALLOC(id_nondiag_elem_bas, num, ID_INT,
                                            func,"new id_nondiag_bas[]");
    if(op_ang->id_nondiag_bas != id_nondiag_elem_bas) { 
      printf("\n\nRank%d: Error in function  id_two_part_angJ():",Rank);
      printf("\nfunction REALLOC() does not work properly - move mem loc");
      printf("\nof id_nondiag matrix elements to new place\n");
      MPI_Abort(MPI_COMM_WORLD,Rank);
    }
  }

******************  til HIT  *******/

} // End: function id_two_part_angJ()

     /*
     ** The function                                          
     **              id_m_pot_diag(...)                    
     ** calculates all diagonal m-scheme J**2 matrix elements
     ** <k.l| J**2 |k.l>  for identical particles 
     ** The diagonal matrix elements  for k < l are  stored  in vector: 
     **         table[((2*num_m - k - 3) * k)/2 + l].diag_val
     */

static void id_m_pot_diag(SP_BAS *spBas, MATR_ID *op_ang)
{
  int     k, l, limit;
  double *id_diag;
      
  limit   = spBas->numm_orb - 1; // initialization
  id_diag = op_ang->id_diag;

  for(k = 0; k < limit; k++) {
    for(l = k+1; l <= limit; l++, id_diag++) {
      *id_diag = id_m_scheme_diag_ang(k,l,spBas->mbas);
    } // end l index
  } // end k index 

} // End: function id_m_pot_diag()

     /*
     ** The function                                
     **          id_m_scheme_diag_ang()            
     ** calculates the diagonal two-particle matrix element 
     ** in m-scheme of 
     **      <k,l| J^2 |k,l> for for identical particles.
     */

static double  id_m_scheme_diag_ang(int k, int l, MBAS *mbas)
{
  int      n_k, n_l, l_k, l_l, j_k, j_l, m_k, m_l;
  double   value;

  n_k = mbas[k].osc;
  n_l = mbas[l].osc;
  j_k = mbas[k].j;
  j_l = mbas[l].j;
  l_k = mbas[k].l;
  l_l = mbas[l].l;
  m_k = mbas[k].m;
  m_l = mbas[l].m;

  value = 2 * m_k * m_l;  

  if((n_k == n_l) && (l_k == l_l) && (j_k == j_l))  {
    if(m_k == (m_l + 2))  {
         value -= (double) ((j_k + m_k) * (j_k - m_k + 2));
    }
    else if(m_k == (m_l - 2)) {
      value -= (double) ((j_k - m_k) * (j_k + m_k + 2));
    }
  }
  return (0.25 *  value);

} // End: function id_m_scheme_diag_ang()

     /*
     ** The function                                          
     **             id_m_pot_nondiag(,...)                    
     ** calculates all nondiagonal m-scheme two-particle matrix elements
     ** <i,j| J**2 |k.l>  for identical particles
     ** The nondiagonal matrix elements  are  stored in a table pointed to by  
     ** table[((2*num_m - k - 3) * k)/2 + (l - 1)]
     ** The function returns number of calculated matrix elements.
     */

static int id_m_pot_nondiag(SP_BAS *spBas,MATR_ID *op_ang)
{
  int      i, j, k, l, m_kl, parity_kl, limit, count, number;
  ID_INT   *nondiag, **nondiag_table = NULL;

  nondiag_table = op_ang->id_nondiag_table;
  op_ang->maxNondiagElem = 0;

  limit   = spBas->numm_orb - 1; // initialization
  nondiag = op_ang->id_nondiag_bas;
  number  = 0;

  for(k = 0; k < limit; k++) { // first right */ 
    for(l = k+1; l <= limit; l++, nondiag_table++)  { // second right
      m_kl = spBas->mbas[k].m + spBas->mbas[l].m;
      parity_kl = spBas->mbas[k].par * spBas->mbas[l].par;

      *nondiag_table = nondiag; // store pointer

      for(count = 0,j = l; j <= limit; j++)  { // second lef
	for(i = (j == l) ? k + 1: 0; i < j; i++){  // first left
	  if(   (m_kl == spBas->mbas[i].m + spBas->mbas[j].m) // test mj and parity
		&& (parity_kl == spBas->mbas[i].par * spBas->mbas[j].par))  {

   nondiag->val[0] = id_m_scheme_nondiag_ang(i,j,k,l, spBas->mbas);

	    if(fabs(nondiag->val[0]) > MATRIX_LIMIT) {
	      nondiag->one    = ULL_ONE<<j ^ ULL_ONE<<i;
	      nondiag->two    = (ULL_ONE<<j) - (ULL_ONE<<(i+1));
	      count++;
	      number++;  
	      nondiag++;
	    }
	  } // end mj and parity test
		     	 	    
	} // end of index i, second particle left */
      } /* end of index j, first particle left */
      if(!count)  {                  /* no matrix elements for present (k,l) */
	*nondiag_table = (ID_INT *)NULL_PTR;
      }
      else   {                      /* count matrix elements for present (k,l) */
	op_ang->maxNondiagElem = MAX(op_ang->maxNondiagElem, count);
	number++; 
	nondiag->one = ULL_ZERO;    /* terminate each (ij)-set with a ZERO */
	nondiag++;
      }
    } /* end of index l, second right */
  } /* end of  index k, first right */

  return number;     /* return the total number of elements in nondiag[] */

} /* End: function id_m_pot_nondiag() */


     /*
     ** The function                                
     **       id_m_scheme_nondiag_ang()            
     ** calculates and returns non-diagonal two-particle 
     ** matrix element in m-scheme 
     **             <i,j| J**2 |k,l>
     ** for the identical particle case.
     ** It is assumed that i < j and k < l.         
     */

static double id_m_scheme_nondiag_ang(int i, int j, int k, int l,
                                                             MBAS *mbas)
{
  int  n_i, n_j, l_i, l_j, j_i, j_j, m_i, m_j,
        n_k, n_l, l_k, l_l, j_k, j_l, m_k, m_l,
        phase, value;

  n_i = mbas[i].osc;  // i and k are particle one 
  n_k = mbas[k].osc;
  l_i = mbas[i].l;
  l_k = mbas[k].l;
  j_i = mbas[i].j;
  j_k = mbas[k].j;
  m_i = mbas[i].m;
  m_k = mbas[k].m;

  n_j = mbas[j].osc;  // j and l are particle two
  n_l = mbas[l].osc;
  l_j = mbas[j].l;
  l_l = mbas[l].l;
  j_j = mbas[j].j;
  j_l = mbas[l].j;
  m_j = mbas[j].m;
  m_l = mbas[l].m;

  value = 0;

  // "Whitehead" phase is included

  phase =   mbas[i].phase * mbas[j].phase
          * mbas[k].phase * mbas[l].phase;

  if(   ((n_i == n_k) && (l_i == l_k) && (j_i == j_k)) 
     && ((n_j == n_l) && (l_j == l_l) && (j_j == j_l)))  {

    // (n_i,j_i ) = ( n_k,j_k) & (n_j,j_j) = (n_l,j_l) 

             /* m_i = m_k - 1  m_j = m_l + 1 */
    if((m_i == (m_k - 2)) && (m_j == (m_l + 2)))  {
      value = (j_k + m_k)*(j_k - m_k + 2)*(j_l - m_l)*(j_l + m_l + 2);
    }
             /* m_i = m_k + 1  m_j = m_l - 1 */
    else if((m_i == (m_k + 2)) && (m_j == (m_l - 2))) { 
         value = (j_k - m_k)*(j_k + m_k + 2)*(j_l + m_l)*(j_l - m_l + 2);
    }
  }

  if(   ((n_i == n_l) && (l_i == l_l) && (j_i == j_l)) 
      && ((n_j == n_k) && (l_j == l_k) && (j_j == j_k)))  {

    // (n_i,j_i ) = ( n_l,j_l) & (n_j,j_j) = (n_k,j_k)

                  /* m_i = m_l + 1  m_j = m_k - 1 */

      if((m_i == (m_l + 2)) && (m_j == (m_k - 2)))  { 
         value = (j_k + m_k)*(j_k - m_k + 2)*(j_l - m_l)*(j_l + m_l + 2);
         phase *= -1;
      } 

                  /* m_i = m_l - 1  m_j = m_k + 1 */
      if((m_i == (m_l - 2)) && (m_j == (m_k + 2)))  { 
	value = (j_k - m_k)*(j_k + m_k + 2)*(j_l + m_l)*(j_l - m_l + 2);
	phase *= -1;
      }
  }
  if(abs(value))      return (0.25 * phase * sqrt((double) value));
  else                return 0.0;

} // End: function id_m_scheme_nondiag_ang() 

     /*
     ** The function
     **        addID_1to2SinglePartTerms()
     ** adds contributions from single-particle terms of J**2
     ** to the two-particle matrix elements  <k.l |J**2| k.l>
     */

static void addID_1to2SinglePartTerms(SP_BAS *spBas,MATR_ID *op_ang)
{
   int         k, l, limit;
   double      *diag, num;

   limit  = spBas->numm_orb - 1;  // max. k - values
   num    = 1.0 /((double)(spBas->part - 1));
   diag   = op_ang->id_diag;
   for(k = 0; k < limit; k++) { 
      for(l = k+1; l <= limit; l++, diag++) {   
	*diag += 0.25*(double)(spBas->mbas[k].j * (spBas->mbas[k].j + 2)
                             + spBas->mbas[l].j*(spBas->mbas[l].j + 2)) * num;
      } // end l-loop
   } // end k-loop

} // End: function  addID_1to2SinglePartTerms()
