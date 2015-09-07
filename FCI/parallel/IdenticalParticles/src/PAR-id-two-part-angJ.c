   /******************  The module PAR-id-two-part-angJ.c  ******************/

#include "PAR-shell-model.h"

     /*
     ** The mdule entrance function 
     **         id_ang_mom_interaction()
     ** calculates the two-particle m-scheme matrix 
     ** elements of J**2 and store the
     ** result in table: op_ang->id_nondiag_table[]
     */

//***** local function declarations  ********** 

static int max_nondiag_id_elem(SP_BAS *sp_bas);
      /*
      ** calculates and returns the maximum  possible number non-diagonal
      ** two-body matrix elements for identical particles.
      */

static void id_m_pot_diag(SP_BAS *sp_bas, MATR_OP *op_ang);
      /*
      ** calculates all diagonal m-scheme two-particle matrix elements
      ** <k.l| J**2 |k.l>  for identical particles 
      ** The diagonal matrix elements  for k < l are  stored  in vector: 
      **         table[((2*num_m - k - 3) * k)/2 + l].diag_val
      */

static double  id_m_scheme_diag_ang_mom(int k, int l, MBAS *mbas);
      /*
      **          id_m_scheme_diag_ang_mom()            
      ** calculates the diagonal two-particle matrix  element in m-scheme of the
      ** angular momentum operator <k,l| J^2 |k,l> for for identical particles.
      */
static int id_m_pot_nondiag(SP_BAS *sp_bas, MATR_OP *op_ang, ID_INT *nondiag_bas);
      /*
      ** calculates all nondiagonal m-scheme two-particle matrix elements
      ** <i,j|J**2|k.l>  for identical particles
      ** The nondiagonal matrix elements  are  stored in a table pointed to by  
      ** table[((2*num_m - k - 3) * k)/2 + (l - 1)]
      ** The function returns number of calculated matrix elements.
      */

static double id_m_scheme_nondiag_ang_mom(int i, int j, int k, int l, MBAS *mbas);
      /*
      ** calculates and returns non-diagonal two-particle matrix element
      ** in m-scheme of the angular momentum  operator
      **             <i,j| J**2 |k,l>
      ** for the identical particle case. It is assumed that i < j and k < l.         
      */

 //***** End local function declarations  ********** 


               /**** The function definitions  ****/ 

     /*
     ** The entrance function 
     **         id_ang_mom_interaction()
     ** calculates the two-particle m-scheme matrix elements 
     ** of J**2 and store the result in table:
     ** op_ang->id_nondiag_table[]
     */

void id_ang_mom_interaction(SP_BAS *sp_bas, MATR_OP *op_ang)
{
   char     *func = {"id_ang_mom_interaction(): "};
   int      num_elem, num;
   ID_INT   *id_nondiag_elem_bas;

           /* memory for  kl-table */

   num_elem = sp_bas->numm_orb * (sp_bas->numm_orb - 1) >> 1;
   op_ang->id_diag = MALLOC(num_elem, double, func, "id_kl_ang_table[]");  
   op_ang->id_nondiag_table
            = MALLOC(num_elem, ID_INT *, func,"kl_ang_nondiag_table[]");

   num_elem            = max_nondiag_id_elem(sp_bas);
   id_nondiag_elem_bas = MALLOC(num_elem,ID_INT,func, "id_effective_nondiag_ang_elements[]");
   op_ang->id_nondiag_bas = id_nondiag_elem_bas;  

             /* calculate and store diag m-scheme of J**2 */

   id_m_pot_diag(sp_bas, op_ang);

          /* calculate and store nondiag m-scheme J**2 matrix elements */

   num = id_m_pot_nondiag(sp_bas, op_ang, id_nondiag_elem_bas); 

   if(num > num_elem) {
      printf("\n\nError in function id_effective_interaction(): ");
      printf(" = %d - calculated = %d",num_elem, num); 
      printf("\n This may not occur!!!!!\n");
       MPI_Abort(MPI_COMM_WORLD,Rank);
   }

/************  NO REALLOC   ****************

   if(num < num_elem) {           // reduce size of nondiag[] to new number
      op_ang->id_nondiag_bas = REALLOC(id_nondiag_elem_bas, num, ID_INT, func,
                               "new id_nondiag_J**2_elements[]");
      if( op_ang->id_nondiag_bas != id_nondiag_elem_bas) { 
	printf("\n\nError in function  id_ang_mom_interaction():");
	printf("\nFunction REALLOC() does not work properly - move memory location");
	printf("\nof id_nondiag matrix elements to a new place\n");
	 MPI_Abort(MPI_COMM_WORLD,Rank);
      }
   }

**************  END NO REALLOC *********/

} // End: function id_ang_mom_interaction()

     /*
     ** The function
     **        max_nondiag_id_elem()
     ** calculates and returns the maximum  possible number non-diagonal
     ** two-body matrix elements for identical particles.
     */

static int max_nondiag_id_elem(SP_BAS *sp_bas)
{
   int
        i, j, k, l, limit, m_kl, parity_kl, count;

   count = 0;                     /* initialization */
   limit = sp_bas->numm_orb - 1;

   for(k = 0; k < limit; k++)  {         /* first index right  */
      for(l = k+1; l <= limit; l++) {     /* second index right */
         m_kl = sp_bas->mbas[k].m + sp_bas->mbas[l].m;        /* m value and parity */
         parity_kl = sp_bas->mbas[k].par * sp_bas->mbas[l].par;
         for(j = l; j <= limit; j++) {                   /* second index left */
            for(i = (j == l) ? k + 1: 0; i < j; i++){     /* first index left */   

                        /* test m value and parity */

               if(   (m_kl == (sp_bas->mbas[i].m + sp_bas->mbas[j].m))
                  && (parity_kl == (sp_bas->mbas[i].par * sp_bas->mbas[j].par))) count++;
            } /* end of index i */
         } /* end of index j */
         count++; /* add space for an extra ZERO */
      } /* end of index l */
   } /* end of index k */
  return count;

} /* End: function max_id_nondiag_elem() */

     /*
     ** The function                                          
     **              id_m_pot_diag(...)                    
     ** calculates all diagonal m-scheme two-particle matrix elements
     ** <k.l| J**2 |k.l>  for identical particles 
     ** The diagonal matrix elements  for k < l are  stored  in vector: 
     **         table[((2*num_m - k - 3) * k)/2 + l].diag_val
     */

static void id_m_pot_diag(SP_BAS *sp_bas, MATR_OP *op_ang)
{
   int
             k, l, limit, number;
   double    *id_diag;

   limit   = sp_bas->numm_orb - 1;          /* initialization */
   id_diag = op_ang->id_diag; 

   number = 0;
   for(k = 0; k < limit; k++) {
      for(l = k+1; l <= limit; l++, id_diag++) {
	*id_diag = id_m_scheme_diag_ang_mom(k,l, sp_bas->mbas);
	number++;
      }   /* end l index */
   }  /* end k index */


} /* End: function id_m_pot_diag() */

     /*
     ** The function                                
     **          id_m_scheme_diag_ang_mom()            
     ** calculates the diagonal two-particle matrix  element in m-scheme of the
     ** angular momentum operator <k,l| J^2 |k,l> for for identical particles.
     */

static double  id_m_scheme_diag_ang_mom(int k, int l, MBAS *mbas)
{
   int      n_k, n_l, l_k, l_l, j_k, j_l, m_k, m_l;
   double                      value;

   n_k = mbas[k].osc;
   n_l = mbas[l].osc;
   l_k = mbas[k].l;
   l_l = mbas[l].l;
   j_k = mbas[k].j;
   j_l = mbas[l].j;
   m_k = mbas[k].m;
   m_l = mbas[l].m;

   value = (double)(2 * m_k * m_l);  

   if((n_k == n_l) && (l_k == l_l) && (j_k == j_l))  {
      if(m_k == (m_l + 2))  {
         value -= (double) ((j_k + m_k) * (j_k - m_k + 2));
      }
      else if(m_k == (m_l - 2)) {
         value -= (double) ((j_k - m_k) * (j_k + m_k + 2));
      }
   }
   return (0.25 *  value);

} /* End: function id_m_scheme_diag_ang_mom() */

     /*
     ** The function                                          
     **             id_m_pot_nondiag(,...)                    
     ** calculates all nondiagonal m-scheme two-particle matrix elements
     ** <i,j|J**2|k.l>  for identical particles
     ** The nondiagonal matrix elements  are  stored in a table pointed to by  
     ** table[((2*num_m - k - 3) * k)/2 + (l - 1)]
     ** The function returns number of calculated matrix elements.
     */

static int id_m_pot_nondiag(SP_BAS *sp_bas, MATR_OP *op_ang, ID_INT *nondiag_bas)
{
  int
            i, j, k, l, m_kl, parity_kl, limit, count, number;
  
  ID_INT
            *nondiag, **nondiag_table;

  nondiag_table = op_ang->id_nondiag_table;
  op_ang->maxNondiagElem = 0;
  limit   =  sp_bas->numm_orb - 1;                           /* initialization */
  nondiag = nondiag_bas;
  number  = 0;
  for(k = 0; k < limit; k++) {                                              /* first right */ 
    for(l = k+1; l <= limit; l++, nondiag_table++)  {                       /* second right */
      m_kl = sp_bas->mbas[k].m + sp_bas->mbas[l].m;
      parity_kl = sp_bas->mbas[k].par * sp_bas->mbas[l].par;

      *nondiag_table = nondiag;                /* store pointer */

      for(count = 0,j = l; j <= limit; j++)  {                    /* second left*/
	for(i = (j == l) ? k + 1: 0; i < j; i++){      /* first left  */
	  if(   (m_kl == sp_bas->mbas[i].m + sp_bas->mbas[j].m) /* test mj and parity */
		&& (parity_kl == sp_bas->mbas[i].par * sp_bas->mbas[j].par))  {

	    nondiag->val = id_m_scheme_nondiag_ang_mom(i,j,k,l, sp_bas->mbas);

	    if(fabs(nondiag->val) > MATRIX_LIMIT) {
	      nondiag->one    = ULL_ONE<<j ^ ULL_ONE<<i;
	      nondiag->two    = (ULL_ONE<<j) - (ULL_ONE<<(i+1));
	      count++;
	      number++;  
	      nondiag++;
	    }
	  } /* end mj and parity test */
		     	 	    
	}  /* end of index i, second particle left */
      } /* end of index j, first particle left */
      if(!count)  {                  /* no matrix elements for present (k,l) */
	*nondiag_table = NULL_PTR;
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
     **       id_m_scheme_nondiag_ang_mom()            
     ** calculates and returns non-diagonal two-particle matrix element
     ** in m-scheme of the angular momentum  operator
     **             <i,j| J**2 |k,l>
     ** for the identical particle case. It is assumed that i < j and k < l.         
     */

static double id_m_scheme_nondiag_ang_mom(int i, int j, int k, int l, MBAS *mbas)
{
   int      n_i, n_j, l_i, l_j, j_i, j_j, m_i, m_j,
            n_k, n_l, l_k, l_l, j_k, j_l, m_k, m_l,
            phase, value;

   n_i = mbas[i].osc;  /* i and k are particle one */
   n_k = mbas[k].osc;
   l_i = mbas[i].l;
   l_k = mbas[k].l;
   j_i = mbas[i].j;
   j_k = mbas[k].j;
   m_i = mbas[i].m;
   m_k = mbas[k].m;

   n_j = mbas[j].osc;  /* j and l are particle two */
   n_l = mbas[l].osc;
   l_j = mbas[j].l;
   l_l = mbas[l].l;
   j_j = mbas[j].j;
   j_l = mbas[l].j;
   m_j = mbas[j].m;
   m_l = mbas[l].m;

   value = 0;

   phase =   mbas[i].phase * mbas[j].phase /* "Whitehead" phase is included */
           * mbas[k].phase * mbas[l].phase;

   if(  ((n_i == n_k) && (l_i == l_k) && (j_i == j_k)) 
     && ((n_j == n_l) && (l_j == l_l) && (j_j == j_l)))  {

      /* (n_i,l_i,j_i ) = ( n_k,l_k,j_k) & (n_j,l_j,j_j) = (n_l,l_l,j_l) */

      if((m_i == (m_k - 2)) && (m_j == (m_l + 2)))  { /* m_i = m_k - 1  m_j = m_l + 1 */
         value = (j_k + m_k)*(j_k - m_k + 2)*(j_l - m_l)*(j_l + m_l + 2);
      }
      else if((m_i == (m_k + 2)) && (m_j == (m_l - 2))) { /* m_i = m_k + 1  m_j = m_l - 1 */
         value = (j_k - m_k)*(j_k + m_k + 2)*(j_l + m_l)*(j_l - m_l + 2);
      }
   }

   if(   ((n_i == n_l) && (l_i == l_l) && (j_i == j_l))
      && ((n_j == n_k) && (l_j == l_k) && (j_j == j_k)))  {

      /** (n_i,l_i,j_i ) = ( n_l,l_l,j_l) & (n_j,l_j,j_j) = (n_k,l_k,j_k) **/

      if((m_i == (m_l + 2)) && (m_j == (m_k - 2)))  { /* m_i = m_l + 1  m_j = m_k - 1 */
         value = (j_k + m_k)*(j_k - m_k + 2)*(j_l - m_l)*(j_l + m_l + 2);
         phase *= -1;
      }
      if((m_i == (m_l - 2)) && (m_j == (m_k + 2)))  { /* m_i = m_l - 1  m_j = m_k + 1 */
         value = (j_k - m_k)*(j_k + m_k + 2)*(j_l + m_l)*(j_l - m_l + 2);
         phase *= -1;
      }
   }
   if(abs(value))      return (0.25 * phase * sqrt((double)value));
   else                return 0.0;

} /* End: function id_m_scheme_nondiag_ang_mom() */
