
/******************  The module PAR-id-two-part-intJ.c  ******************/

#include "PAR-shell-model.h"

      /*
      ** The module entrance function 
      **    id_two_part_intJ()
      ** reads from file angular momentum coupled two-particle
      ** identical particle matrix elements, transfer the elements 
      ** into m-scheme representation and store the result in tables:
      ** ID_TABLE id_table[]
      */

// data definitions for two-particle matrix elements

     /* Storage of the spherical j-j coupled effective
     ** matrix elements for identical particles
     **         <n_a,n_b | V | n_c, n_d>
     */
 
/*******************   local data structure  ***********/

   typedef struct  {
          int           J;  // twice total J
          ULL      config;  // = n_a << 3*length + n_b << 2*length
			    // + n_c<<length + n_d,
	                    // with length = sizeof(ULL) 
          double value[2];  // value of veff element[0] and 
                            // center_of_mass element [1]
   } TWO_PART;

   typedef struct  {
          int            num,  // total number of matrix elements
                      *start;  // start position of given J group
          TWO_PART *two_part;  // pointer to the matrix elements
   } J_POT;

                  /*
                  ** local function declarations in 
		  **      file lanc-input.c
                  */
static void read_id_two_part_matr_elem(char *filename, SP_BAS *sp_bas,
				                         J_POT *j_pot);
     /*
     ** reads angular momentum coupled two-particle matrix elements
     ** for identical particles. Each matrix element
     ** is given an identity number. Then they are  sorted into groups of  total_J 
     ** in a monotonic increasing sequence and finally stored in structure J_POT j_pot[].
     */

static ORB_ID id_matr_elem_identity(SP_BAS *sp_bas, int *int_data);
      /*
      ** converts particle matrix elements identification from
      ** (n,l,j) form into an orbital number specified according
      ** to sp_bas->jbas[]. The function returns four orbital
      ** numbers in structure ORB_ID orb;  
      */

static double id_orbit_interchange(SP_BAS *sp_bas, int tot_J, ORB_ID *orb, int length);
      /*
      ** reorder in case of identical particles orb_a, orb_b
      ** orb_c and orb_d such that    
      **        orb_a <= orb_b , orb_c <= orb_d                   
      ** if(orb_a == orb_b) a normalization factor of sqrt2 is produced
      ** to obtain unnormalized two-particle matrix elements. A similar
      ** factor is added if(orb_c == orb_d). Calculate and store in
      ** structure orb_id left and right identity number.
      ** The function returns the corresponding phase_norm factor.               
      */

static void sort_matrix_elements(int num_J, J_POT *j_pot);
      /*
      ** sorts all matrix elements into groups of same total J and
      ** for each group sort after increasing configuration number
      */

static int total_J_comp(const TWO_PART *one, const TWO_PART *two);
      /*
      ** is a utility function for the library function qsort() in order to sort the
      ** two-body matrix elements into groups of increasing value of the total J-value.
      */

static int config_comp(const TWO_PART *one, const TWO_PART *two);
      /*
      ** is a utility function for the library function qsort() in order to sort the two-body
      ** matrix elements for fixed total J into groups of increasing conf_ident number.
      */

static void check_id_matr_elem(SP_BAS *sp_bas, int num_J, J_POT j_pot);
      /*
      ** checks current identical two-particle matrix elements to
      ** see if there are any missing compared to the available
      ** single-particle orbits in the model.
      */

static int max_nondiag_id_elem(SP_BAS *sp_bas);
      /*
      ** calculates and returns the maximum  possible number non-diagonal
      ** two-body matrix elements for identical particles.
      */

static void id_m_pot_diag(SP_BAS *sp_bas, MATR_OP *op_int, J_POT j_pot);
      /*
      ** calculates all diagonal m-scheme two-particle matrix elements
      ** <k.l| V |k.l>  for identical particles of an interaction
      ** determined by the parameter 
      **   interaction = 0:  effective interaction 
      **               = 1:  angular momentum interaction
      ** The diagonal matrix elements  for k < l are  stored  in vector: 
      **         table[((2*num_m - k - 3) * k)/2 + l].diag_val
      */

static double  id_m_scheme_diag_int(int m_orb_a, int m_orb_b, MBAS *mbas,
                                                                    J_POT j_pot);
      /*
      ** calculates and return for identical particles diagonal effective interaction
      ** in m-scheme
      **  <m_orb_a, m_orb_b | Int | m_orb_a, m_orb_b > =                        
      **    SUM(jtot) (  
      **       clebsch(j_a, m_a, j_b, m_b, jtot) * clebsch(j_a, m_a, j_b, m_b, jtot)         
      **     * <j_orb_a, j_orb_b, jtot| Int | j_orb_a, j_orb_b, jtot >
      */

static int id_m_pot_nondiag(SP_BAS *sp_bas, MATR_OP *op_int,  J_POT j_pot, ID_INT *nondiag_bas);
      /*
      ** calculates all nondiagonal m-scheme two-particle matrix elements
      ** <i,j| V |k.l>  for identical particles of an interaction
      ** determined by the parameter 
      **   interaction = 0:  effective interaction 
      **               = 1:  angular momentum interaction
      ** The nondiagonal matrix elements  are  stored in a table pointed to by  
      ** table[((2*num_m - k - 3) * k)/2 + l].start_nondiag
      ** The function returns number of calculated matrix elements.
      */
  
static double  id_m_scheme_nondiag_int(int m_orb_a, int m_orb_b, int m_orb_c,
                                               int m_orb_d, MBAS *mbas, J_POT j_pot);
      /*
      ** calculates and return effective interaction in m-scheme (a <= b and c <= d) 
      **   <m_orb_a, m_orb_b | Int | m_orb_c, m_orb_d > =                        
      **      SUM(jtot) (  clebsch(j_a, m_a, j_b, m_b, jtot)         
      **                  * clebsch(j_c, m_c, j_d, m_d, jtot)         
      **                  * <j_orb_a, j_orb_b, jtot| Int | j_orb_c, j_orb_d, jtot >
      */

                /**** End: function declarations ****/


               /**** The function definitions  ****/ 

     /*
     ** The entrance function 
     **         id_two_part_intJ()
     ** reads from file angular momentum coupled two-particle
     ** identical particle matrix elements, transfer the elements 
     ** into m-scheme representation and store the result in table:
     ** op_int->id_nondiag_table[]
     */

  void id_two_part_intJ(char *filename, SP_BAS *sp_bas, MATR_OP *op_int)
{
   char     *func = {"id_two_part_intJ(): "};
   int      num_elem, num;
   J_POT    j_pot;
   ID_INT   *id_nondiag_elem_bas;

   // memory for  kl-table

   num_elem = (sp_bas->numm_orb * (sp_bas->numm_orb - 1)) >> 1;
   op_int->id_diag
            = MALLOC(num_elem, double, func, "id_kl_int_diag[]");  
   op_int->id_nondiag_table
            = MALLOC(num_elem, ID_INT *, func,"id_kl_int_nondiag_table[]");

   num_elem       = max_nondiag_id_elem(sp_bas);
   id_nondiag_elem_bas = MALLOC(num_elem, ID_INT, func,"id_nondiag_int_elements[]");
   op_int->id_nondiag_bas = id_nondiag_elem_bas;
 
   // read and store two-particle matrix elements in J-scheme

   read_id_two_part_matr_elem(filename, sp_bas, &j_pot);

   // calculate and store diag m-scheme two-particle matrix elements

   id_m_pot_diag(sp_bas, op_int, j_pot);

   // calculate and store nondiag m-scheme effective two-particle matrix elements

   num = id_m_pot_nondiag(sp_bas,op_int, j_pot, id_nondiag_elem_bas); 

   if(num > num_elem) {
     printf("\n\nError in function id_effective_interaction(): ");
     printf(" = %d - calculated = %d",num_elem, num); 
     printf("\n This may not occur!!!!!\n");
   MPI_Abort(MPI_COMM_WORLD,Rank);
   }
   if(num == 0) {
     free( op_int->id_nondiag_bas);
     op_int->id_nondiag_bas   = NULL_PTR;
     free(op_int->id_nondiag_table);
     op_int->id_nondiag_table = NULL_PTR;
     op_int->num_nondiag_elem = 0;
   }
   else if(num < num_elem) {           /* reduce size of nondiag[] to new number */
      op_int->id_nondiag_bas = REALLOC(id_nondiag_elem_bas, num, ID_INT, func,
                                                     "new id_nondiag_bas[]");
      if(op_int->id_nondiag_bas != id_nondiag_elem_bas) { 
	printf("\n\nError in function id_effective_interaction():");
	printf("\nFunction REALLOC() does not work properly - move memory location");
	printf("\nof id_nondiag matrix elements to a new place\n");
        MPI_Abort(MPI_COMM_WORLD,Rank);
      }
   }
   free(j_pot.start);        /* release temporary memory in pp_pot */
   free(j_pot.two_part);  

} // End: function id_two_part_intJ()

     /*
     ** The function
     **         read_id_two_part_matr_elem()
     ** reads angular momentum coupled two-particle matrix elements
     ** for identical particles. Each matrix element
     ** is given an identity number. Then they are  sorted into groups of  total_J 
     ** in a monotonic increasing sequence and finally stored in structure J_POT j_pot[].
     */

static void read_id_two_part_matr_elem(char *filename, SP_BAS *sp_bas, J_POT *j_pot)
{
   char      ch, *func = {"read_id_two_part_matr_elem(): "};
   int       loop, int_data[14], num_elem, new_num_elem, num_J;
   ULL        length, max_orb, temp;
   double    dbl_val;
   ORB_ID    orb_ident;
   TWO_PART  *two_part, *ptr;
   FILE      *file_ptr;

   length  = sizeof(ULL) << 1;       /* number of bits for indentity number */
   max_orb = (ULL_ONE << length);     /* maximum number of spherical orbits */

   if(sp_bas->numj_orb > max_orb)   {
      printf("\n\nError in function read_id_two_part_matr_elem():");
      printf("\nToo many spherical particle single-particle orbits"); 
      printf("\nsp_bas->numj_orb = %d  -- must not exceed %d!",
                                 sp_bas->numj_orb, (int)max_orb);
      printf("\n since an integer is only %d bytes\n", (int)sizeof(int));
    MPI_Abort(MPI_COMM_WORLD,Rank);
   }
   if( (file_ptr = fopen(filename,"r")) == NULL) {              /* open input file */
      printf("\n\nError in function read_id_two_part_matr_elem();");
      printf("\nWrong file = %s for input of eff.interaction.\n",filename);
    MPI_Abort(MPI_COMM_WORLD,Rank);
   }


   /*****************************************************
   if((fscanf(file_ptr,"%d%c", &num_elem, &ch) != 2) | (ch !=  '\n')) { 
      printf("\n\n Error in function read_two_part_matr_elem(): ");
      printf("\nFirst number in file %s", filename);
      printf(" - number of matrix element - is not correct\n");
    MPI_Abort(MPI_COMM_WORLD,Rank);
   } 
   *************************************/

   fscanf(file_ptr,"%d",&num_elem);

            /* memory to store two-particle matrix elements */

   two_part = MALLOC(num_elem, TWO_PART, func,"j_pot.two_part[]");
   
   for(loop = 0, new_num_elem = 0, ptr = two_part; loop < num_elem; loop++) {

        /* read (n1,l1,j1), (n2,l2,j2), (n3,l3,j3), (n4,l4,j4), MT, 2*J */

     if(read_int_number(file_ptr, 14, int_data) == FALSE) {
       printf("\n\nError in function read_two_part_matr_elem();");
       printf("\nSomething is wrong with the (n,l,j) - values");
       printf(" for matrix element number = %d\n",loop);
       MPI_Abort(MPI_COMM_WORLD,Rank);
     }
     // read two-particle matrix element

     if(read_float_number(file_ptr, 1, &dbl_val) == FALSE) {
       printf("\n\nError in function read_two_part_matr_elem();");
       printf("\nWrong value for two-particle matrix element no %d",loop);
       printf("\nfor the case without Center_of_Mass effect");
       printf("\nread from file %s\n", filename);
       MPI_Abort(MPI_COMM_WORLD,Rank);
     }

            /*
	    ** converts (n,l.j) input data to orbit numbers from
            ** jbas[] and store the result in structure orb_ident
            */ 

     orb_ident = id_matr_elem_identity(sp_bas, int_data);

            /*
	    **  Only matrix elements relevant for the present
            **  set of single-particle orbits are saved
            */

     if(  (orb_ident.a == -1) || (orb_ident.b== -1)
	||(orb_ident.c == -1) || (orb_ident.d == -1)) continue;
                  
            /*
	    ** if necessary reorder orbit numbers and change normalization
            ** factors which may modify the matrix element. Calculate and
            ** store in structure orb_ident left and right identity number.
            */

     dbl_val *= id_orbit_interchange(sp_bas, int_data[13], &orb_ident, length);

            /* if necessary, interchange right and
	    ** left orbit identity such that
            **         orb_ident.left <= orb_ident.right
            */
       
     if( orb_ident.left > orb_ident.right) {
       temp          = orb_ident.left;
       orb_ident.left  = orb_ident.right;
       orb_ident.right = temp; 
     }
            /* store config number together with matrix element */

     ptr->config   = (orb_ident.left << (2*length)) + orb_ident.right;
     ptr->J        = int_data[13];
     ptr->value[0] = dbl_val;
     ptr++;
     new_num_elem++;

   } /* end loop through all input two-particle matrix elements */

   fclose(file_ptr);

   if(new_num_elem < num_elem) {           /* reduce size of j_pot.two_part[] */
      two_part = REALLOC(two_part, new_num_elem, TWO_PART, func,"realloc j_pot.two_part[]");
   }
   j_pot->num      = new_num_elem;                        /* save number of matrix element */
   j_pot->two_part = two_part;                           /* and the corresponding pointer */
 
            /* 
            ** All matrix elements are sorted into groups of same total J
            ** and for each group after increasing conf_ident.
            */

   num_J = sp_bas->jbas[0].j + 1;
   sort_matrix_elements(num_J, j_pot);

            /* check that all necessary two-particle matrix elements are stored */
 
   check_id_matr_elem(sp_bas,num_J, *j_pot);

} /* End: function  read_id_two_part_matr_elem() */

     /*
     ** The function 
     **         id_matr_elem_identity()
     ** converts particle matrix elements identification from
     ** (n,l,j) form into an orbital number specified according
     ** to JBAS jbas[]. The function returns four orbital
     ** numbers in structure ORB_ID orb;  
     */

static ORB_ID id_matr_elem_identity(SP_BAS *sp_bas, int *int_data)
{
   int        k;
   JBAS       *ptr;
   ORB_ID     orb;
 
   orb.a= -1;                                          /* local initialization */
   orb.b= -1;
   orb.c= -1;
   orb.d= -1;

   for(k = 0, ptr = sp_bas->jbas; k < sp_bas->numj_orb; k++, ptr++)     {
      if(   (ptr->osc == int_data[0]) && (ptr->l == int_data[1]) 
         && (ptr->j == int_data[2]))  orb.a = k;
      if(   (ptr->osc == int_data[3]) && (ptr->l == int_data[4]) 
         && (ptr->j == int_data[5]))  orb.b = k;
       if(   (ptr->osc == int_data[6]) && (ptr->l == int_data[7]) 
         && (ptr->j == int_data[8]))  orb.c = k;
      if(   (ptr->osc == int_data[9]) && (ptr->l == int_data[10]) 
         && (ptr->j == int_data[11])) orb.d = k;
   }

   return orb;

} /* End: function id_matr_elem_identity() */

     /*
     ** The function
     **          id_orbit_interchange()
     ** reorder in case of identical particles orb_a, orb_b
     ** orb_c and orb_d such that    
     **        orb_a <= orb_b , orb_c <= orb_d                   
     ** if(orb_a == orb_b) a normalization factor of sqrt2 is produced
     ** to obtain unnormalized two-particle matrix elements. A similar
     ** factor is added if(orb_c == orb_d). Calculate and store in
     ** structure orb_id left and right identity number.
     ** The function returns the corresponding phase_norm factor.               
     */

static double id_orbit_interchange(SP_BAS *sp_bas, int tot_J, 
                                             ORB_ID *orb, int length)
{
   int        temp;
   double     phase, norm, sqr2 = 1.41421356237; /* square root of 2 */;

   phase = +1.0;      
   norm  = 1.0;                                   /* initialization */

   if(orb->a > orb->b)  {
      temp    = orb->a;
      orb->a  = orb->b;
      orb->b  = temp;
      phase *= (-1.0)*PHASE((sp_bas->jbas[orb->a].j + sp_bas->jbas[orb->b].j - tot_J) >> 1);
   }

   if(orb->c > orb->d)  {
      temp    = orb->c;
      orb->c  = orb->d;
      orb->d  = temp;
      phase *= (-1.0)*PHASE((sp_bas->jbas[orb->c].j + sp_bas->jbas[orb->d].j - tot_J) >> 1);
   }
             /* not normalized matrix elements */

   if(orb->a == orb->b)  norm *= sqr2;
   if(orb->c == orb->d)  norm *= sqr2;

           /* calculate  config number for matrix element */

   orb->left  = (ULL)((orb->a << length) + orb->b);
   orb->right = (ULL)((orb->c << length) + orb->d);       

   return (phase * norm);

} /* End: function id_orbit_interchange() */

     /* 
     ** The function
     **        sort_matrix_elements()             
     ** sorts all matrix elements into groups of same total J and
     ** for each group sort after increasing configuration number
     */

static void sort_matrix_elements(int num_J, J_POT *j_pot)
{
   char    *func = {"sort_matrix_elements: "};    
   int     jtot, count, temp; 

        /* sort matrix elements into groups after increasing total J-value */

   qsort(j_pot->two_part, (size_t) j_pot->num, (size_t)sizeof(TWO_PART),
                      (int(*)(const void *, const void *)) total_J_comp);

             /* memory for j_pot->start[] */

   j_pot->start = MALLOC(num_J + 1, int, func, "j_pot.start[]");

   for(jtot = 0, count = 0; jtot < num_J; jtot++) {
      j_pot->start[jtot] = count;
      while((count < j_pot->num) && (j_pot->two_part[count].J == (jtot << 1))) count++;

         /* sort the matrix elements with J = tot_J after increasing config */

      if((temp = count - j_pot->start[jtot]) > 0)   {
         qsort(j_pot->two_part + j_pot->start[jtot], (size_t)(temp), (size_t)sizeof(TWO_PART),
                       (int(*)(const void *, const void *)) config_comp);
      }
   }
   j_pot->start[jtot] = count;
   
} /* End: function sort_matrix_elements() */

     /*
     ** The function                         
     **        int total_J_comp()                  
     ** is a utility function for the library function qsort() in order to sort the
     ** two-body matrix elements into groups of increasing value of the total J-value.
     */

static int total_J_comp(const TWO_PART *one, const TWO_PART *two)
{
  if(one->J > two->J)       return +1;
  else  if(one->J < two->J) return -1;
  else                        return  0;

} /* End: function total_J_comp() */

     /*
     ** The function                         
     **        int config_comp()                  
     ** is a utility function for the library function qsort() in order
     ** to sort the two-body matrix elements for fixed total J into 
     ** groups of increasing conf_ident number.
     */

static int config_comp(const TWO_PART *one, const TWO_PART *two)
{
  if(one->config > two->config)       return +1;
  else  if(one->config < two->config) return -1;
  else                      return  0;

} /* End: function config_comp() */

     /*
     ** The function 
     **        check_id_matr_elem()
     ** checks current identical two-particle matrix elements to
     ** see if there are any missing compared to the available
     ** single-particle orbits in the model.
     */

static void check_id_matr_elem(SP_BAS *sp_bas, int num_J, J_POT j_pot)
{
   char     *func = {"check_id_matr_elem: "};
   int      ja_orb, jb_orb, j_min, j_max, jtot, temp, control,
            parity,  *states_J_plus, *states_J_minus, *ptr;

   states_J_plus  = CALLOC(num_J, int, func, "id_num_J_plus[]");    /* local memory */
   states_J_minus = CALLOC(num_J, int, func, "id_num_J_minus[]");

   for(ja_orb = 0; ja_orb < sp_bas->numj_orb; ja_orb++) {
      for(jb_orb = ja_orb; jb_orb <sp_bas->numj_orb; jb_orb++) {
         j_min = (abs(sp_bas->jbas[ja_orb].j - sp_bas->jbas[jb_orb].j)) >> 1;
	 j_max = (sp_bas->jbas[ja_orb].j + sp_bas->jbas[jb_orb].j) >> 1;
         parity = ((sp_bas->jbas[ja_orb].l +sp_bas->jbas[jb_orb].l) % 2) ? -1 : +1;
         ptr    = (parity == +1) ? states_J_plus : states_J_minus;
         temp = (ja_orb == jb_orb) ? TRUE : FALSE;
         for(jtot = j_min; jtot <= j_max; jtot++)  {
	    if(temp && (jtot % 2)) continue;                    /* no |(jj)J = odd> */
            ptr[jtot]++;
         }
      } /* end orb_b */
   } /* end orb_a */

   for(jtot = 0; jtot < num_J; jtot++)  {
      control = j_pot.start[jtot + 1] - j_pot.start[jtot];
      temp    = (states_J_plus[jtot] * (states_J_plus[jtot] + 1)) >> 1;
      temp   += (states_J_minus[jtot] * (states_J_minus[jtot] + 1)) >> 1;

      if(control != temp) {
         printf("\n\nError in function check_id_matr_elem(): ");
         printf("\nNot enough two-particle matrix elements input file");
         printf(" with total_J = %d", jtot);
         printf("\nRead number = %d - should have been = %d\n", control, temp);
       MPI_Abort(MPI_COMM_WORLD,Rank);
      }
   }     

   free(states_J_minus);                         /* release local memory */
   free(states_J_plus);
 
} /* End: function check_id_matr_elem() */

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
     ** <k.l| INT |k.l>  for identical particles 
     ** The diagonal matrix elements  for k < l are  stored  in vector: 
     **         table[((2*num_m - k - 3) * k)/2 + l].diag_val
     */

static void id_m_pot_diag(SP_BAS *sp_bas, MATR_OP *op_int, J_POT j_pot)
{
   int
             k, l, limit;
   double
             *id_diag;
      
   limit   = sp_bas->numm_orb - 1;                                     /* initialization */
   id_diag = op_int->id_diag;
   for(k = 0; k < limit; k++) {
      for(l = k+1; l <= limit; l++, id_diag++) {
	*id_diag = id_m_scheme_diag_int( k,l, sp_bas->mbas, j_pot);
      }   /* end l index */
   }  /* end k index */

} /* End: function id_m_pot_diag() */

     /*
     ** The function                                  
     **               id_m_scheme_diag_int()                   
     ** calculates and return for identical particles diagonal effective interaction
     ** in m-scheme
     **  <m_orb_a, m_orb_b | Int | m_orb_a, m_orb_b > =                        
     **    SUM(jtot) (  
     **       clebsch(j_a, m_a, j_b, m_b, jtot) * clebsch(j_a, m_a, j_b, m_b, jtot)         
     **     * <j_orb_a, j_orb_b, jtot| Int | j_orb_a, j_orb_b, jtot >
     */

static double  id_m_scheme_diag_int(int m_orb_a, int m_orb_b, MBAS *mbas, J_POT j_pot)
{
   int       j_a, m_a, j_b, m_b, j_min, j_max, jtot;
   ULL        j_orb_a, j_orb_b, temp, length, config, high, low, mid;
   double    int_m, cleb;

   int_m = 0.0;                                /* initialization */
   j_orb_a = (ULL) mbas[m_orb_a].orb;                        /* orbit numbers */
   j_orb_b = (ULL) mbas[m_orb_b].orb;
   j_a = mbas[m_orb_a].j;                 /* j and m values */
   m_a = mbas[m_orb_a].m;
   j_b = mbas[m_orb_b].j;
   m_b = mbas[m_orb_b].m;
   j_min = (abs(j_a - j_b)) >>1;
   j_max = (j_a + j_b) >> 1;
   if(j_orb_a > j_orb_b )  {     /* interchange left and right side of matr. elem. */
      temp    = j_orb_a;
      j_orb_a = j_orb_b;
      j_orb_b = temp;
   }
   length = sizeof(ULL) << 1;
   config = (j_orb_a << (3 * length)) + (j_orb_b << (2 * length)) 
            + (j_orb_a << length) + j_orb_b;
   temp = (j_orb_a == j_orb_b) ? TRUE : FALSE;
   for(jtot = j_min; jtot <= j_max; jtot++)  {
      if(temp && (jtot % 2)) continue;  /* no contribution for |(jj)J = odd> */
      low     = j_pot.start[jtot];
      high    = j_pot.start[jtot + 1];
      for( ; ; ) {   // binary search for configuration
         mid = (low + high) >> 1;
         if(config < j_pot.two_part[mid].config)      high = mid -1;
         else if(config >j_pot.two_part[mid].config)  low = mid + 1;
         else                                         break;
      }    /* end of search for configuration identifier  */
      cleb =  clebsch_gordan(j_a, j_b, jtot << 1,m_a, m_b);
      int_m += cleb *cleb * j_pot.two_part[mid].value[0];
   }  /* end loop over all jtot */
   return int_m;

} /* End: function m_scheme_diag_int() */

     /*
     ** The function                                          
     **             id_m_pot_nondiag(,...)                    
     ** calculates all nondiagonal m-scheme two-particle matrix elements
     ** <i,j| INT |k.l>  for identical particles of
     ** The nondiagonal matrix elements  are  stored in a table pointed to by  
     ** table[((2*num_m - k - 3) * k)/2 + (l - 1)]
     ** The function returns number of calculated matrix elements.
     */

static int id_m_pot_nondiag(SP_BAS *sp_bas, MATR_OP *op_int,  J_POT j_pot, ID_INT *nondiag_bas)
{
  int
            i, j, k, l, m_kl, parity_kl, limit, count, number;
  
  ID_INT
            *nondiag, **nondiag_table = NULL;

  nondiag_table = op_int->id_nondiag_table;
  op_int->maxNondiagElem = 0;

  limit   =  sp_bas->numm_orb - 1;                           /* initialization */
  nondiag = nondiag_bas;
  number  = 0;

  for(k = 0; k < limit; k++) {                            /* first right */ 
    for(l = k+1; l <= limit; l++, nondiag_table++)  {                       /* second right */
      m_kl = sp_bas->mbas[k].m + sp_bas->mbas[l].m;
      parity_kl = sp_bas->mbas[k].par * sp_bas->mbas[l].par;

      *nondiag_table = nondiag;                /* store pointer */

      for(count = 0,j = l; j <= limit; j++)  {                    /* second left*/
	for(i = (j == l) ? k + 1: 0; i < j; i++){      /* first left  */
	  if(   (m_kl == sp_bas->mbas[i].m + sp_bas->mbas[j].m) /* test mj and parity */
		&& (parity_kl == sp_bas->mbas[i].par * sp_bas->mbas[j].par))  {

   nondiag->val = id_m_scheme_nondiag_int(i,j,k,l, sp_bas->mbas, j_pot);

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
	op_int->maxNondiagElem = MAX(op_int->maxNondiagElem, count);
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
     **               id_m_scheme_nondiag_int()                   
     ** calculates and return effective interaction in m-scheme (a <= b and c <= d) 
     **   <m_orb_a, m_orb_b | Int | m_orb_c, m_orb_d > =                        
     **      SUM(jtot) (  clebsch(j_a, m_a, j_b, m_b, jtot)         
     **                  * clebsch(j_c, m_c, j_d, m_d, jtot)         
     **                  * <j_orb_a, j_orb_b, jtot| Int | j_orb_c, j_orb_d, jtot >
     */

static double  id_m_scheme_nondiag_int(int m_orb_a, int m_orb_b, int m_orb_c,
                                             int m_orb_d, MBAS *mbas, J_POT j_pot)
{
   int       j_a, m_a, j_b, m_b, j_c, m_c, j_d, m_d,
             phase, j_min, j_max, jtot, type, length;
   ULL        j_orb_a, j_orb_b, j_orb_c, j_orb_d, conf_bra,
             temp, conf_ket, config, high, low, mid;
   double    int_m;

   int_m = 0.0;                                /* initialization */
   
   j_orb_a = (ULL)mbas[m_orb_a].orb;                        /* orbit numbers */
   j_orb_b = (ULL)mbas[m_orb_b].orb;
   j_orb_c = (ULL)mbas[m_orb_c].orb;
   j_orb_d = (ULL)mbas[m_orb_d].orb;  

   j_a = mbas[m_orb_a].j;                 /* j and m values */
   m_a = mbas[m_orb_a].m;
   j_b = mbas[m_orb_b].j;
   m_b = mbas[m_orb_b].m;
   j_c = mbas[m_orb_c].j;
   m_c = mbas[m_orb_c].m;
   j_d = mbas[m_orb_d].j;
   m_d = mbas[m_orb_d].m;

   j_min = (MAX(abs(j_a - j_b), abs(j_c - j_d))) >> 1; /* max and min values of total J*2 */ 
   j_max = (MIN(j_a + j_b, j_c + j_d)) >> 1;

   if(j_min > j_max)  return int_m ;

   type  = 0;
   phase = +1;

   if(j_orb_a > j_orb_b )  {               /* interchange left side of matr. elem. */
      temp    = j_orb_a;
      j_orb_a = j_orb_b;
      j_orb_b = temp;
      type    = 1;
      phase   = PHASE((j_a + j_b - 2*j_min)/2 + 1);
   }
   if(j_orb_c > j_orb_d )  {               /* interchange right side of matr. elem. */
      temp    = j_orb_c;
      j_orb_c = j_orb_d;
      j_orb_d = temp;
      type    = type ? 0 : 1;
      phase  *= PHASE((j_c + j_d - 2*j_min)/2 + 1);
   }
   length   = sizeof(ULL) << 1;
   conf_bra = (j_orb_a << length) + j_orb_b;
   conf_ket = (j_orb_c << length) + j_orb_d;

   if (conf_bra > conf_ket) {     /* possible interchange left and right of matr. elem. */
      temp     = conf_bra ;
      conf_bra = conf_ket ;
      conf_ket = temp ;
   }
   config = (conf_bra << (2 * length)) + conf_ket;

   temp = ((j_orb_a == j_orb_b) || (j_orb_c == j_orb_d)) ? TRUE : FALSE;
   for(jtot = j_min; jtot <= j_max; jtot++)  {
      if(!(temp && (jtot % 2)))    {    /* no contribution for |(jj)J = odd> */
         low     = j_pot.start[jtot];
         high    = j_pot.start[jtot + 1];
         for( ; ; )  { //  binary search for configuration
            mid = (low + high) >> 1;
            if(config < j_pot.two_part[mid].config)      high = mid -1;
            else if(config >j_pot.two_part[mid].config)  low = mid + 1;
            else                                             break;
         }    /* end of search for configuration identifier  */

         int_m += phase * clebsch_gordan(j_a, j_b, jtot << 1,m_a, m_b)
                         * clebsch_gordan(j_c, j_d, jtot << 1,m_c, m_d)
                         * j_pot.two_part[mid].value[0];
      } /* end if - loop */

      if(type) phase = -phase;     /* for type = 1 change phase for next jtot */

   }  /* end loop over all jtot */

            /* "Whitehead" phase is included */
 
   return (int_m * mbas[m_orb_a].phase * mbas[m_orb_b].phase
                  * mbas[m_orb_c].phase * mbas[m_orb_d].phase);

} /* End: function id_m_scheme_nondiag_int() */
