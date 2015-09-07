/*******************  The module lanc-pn-pot.c  ******************/

#include "shell.h"

     /*
     ** The module entrance function 
     **  void proton_neutron_effective_interaction(FILE_NAME *file_name, PN_SHELL *model)
     **    1. read angular momentum coupled effective
     **       two-particle matrix elements for from file 
     **    2. calculate and store all m-scheme effective
     **       two-particle matrix elements
     **    3. calculate and store all m-scheme matrix 
     **       elements of J**2
     */

             /**** data definitions ****/

       /* Storage of the spherical j-j coupled effective
        * matrix elements for identical particles
        *         <n_a,n_b | V | n_c, n_d>
        */
 
   typedef struct  {
         int   
                       J;    /* twice total J                       */
          UL
                  config;    /* = n_a << 3*length + n_b << 2*length */
			     /* + n_c<<length + n_d,                */
	                     /* with length = sizeof(UL) >> 2       */
          double
                   value;    /* value of the matrix element         */
   } TWO_PART;

   typedef struct  {
           int
                     num,    /* total number of matrix elements */
                  *start;    /* start position of given J group */
           TWO_PART
               *two_part;    /* pointer to the matrix elements  */
   } J_POT;

   typedef   struct {
           int
                       a,    /* left orbit  - particle one */
                       b,    /* left orbit  - particle two */
                       c,    /* right orbit - particle one */
                       d;    /* right orbit - particle two */
           UL
                    left,    /* left two-particle config   */
                   right;    /* right part-particle config */
   } ORB_ID;

                  /*
                  ** local function declarations in 
		  **      file lanc-input.c
                  */

static void pn_id_effective_interaction(char *file_name, int type, PN_SHELL *model, GR_BASIS *gr_basis);
      /*
      ** reads from file angular momentum coupled two-particle
      ** identical particle matrix elements, transfer the elements 
      ** into m-scheme representation and store the result in tables:
      ** ID_TABLE id_table[type][], type = 0(1) - protons(neutrons)
      */

static void pn_id_ang_mom_interaction(int type, PN_SHELL *model, GR_BASIS *gr_basis);
      /*
      ** calculates m-scheme matrix elements of J**2 and store the
      ** result in tables: ID_TABLE id_table[type][],
      ** type = 0(1) - protons(neutrons)
      */

static void pn_pn_effective_interaction(char *file_name, PN_SHELL *model, GR_BASIS *gr_basis);
      /*
      ** reads from file angular momentum coupled two-particle
      ** proton-neutron matrix elements, transfer the elements into
      ** m-scheme representation and store the result in table:
      **                 PN_TABLE pn_table[].
      */

static void pn_pn_ang_mom_interaction(PN_SHELL *model, GR_BASIS *gr_basis);
      /*
      ** calculates m-scheme proton/neutron matrix elements of J**2
      ** and store the result in tables: PN_TABLE pn_table[].
      */

static void pn_read_id_two_part_matr_elem(char *file_name, int type,int numj_orb,
                                                            JBAS *jbas, J_POT *j_pot);
      /*
      ** reads angular momentum coupled two-particle effective matrix elements
      ** for identical particles, type = 0(1) - proton (neutron). Each matrix element
      ** is given an identity number. Then they are  sorted into groups of  total_J 
      ** in a monotonic increasing sequence and finally stored in structure J_POT j_pot[].
      */
static int read_int_number(FILE *in_data, int number, int *vector);
      /*
      ** reads (number) integers from FILE *in_data stream and
      ** store them in vector[]. If too few integers are read
      ** the function returns  FALSE, otherwise TRUE.
      */

static int read_float_number(FILE *in_data, int number, double *vector);
      /*   
      ** reads (number) float from FILE *in_data stream,
      ** convert them to double and store them in vector[].
      ** If too few doubles are read the function returns 
      ** FALSE, otherwise TRUE.
      */
static int new_read_float_number(FILE *in_data, int number, double *vector);
      /*   
      ** reads (number) float from FILE *in_data stream,
      ** convert them to double and store them in vector[].
      ** If too few doubles are read the function returns 
      ** FALSE, otherwise TRUE.
      */

static int read_data_string(FILE *in_data, char *data_string);
      /*
      ** reads one piece of input data as a text_string containing no blanks. 
      ** Whitespaces in front of the data_string is removed. Data_string must
      ** not contain the character '\'. The function returns FALSE if no 
      ** proper data_string is found, otherwise TRUE.
      */

static ORB_ID pn_id_matr_elem_identity(int numj_orb, JBAS *jbas, int *int_data);
      /*
      ** converts particle matrix elements identification from
      ** (n,l,j) form into an orbital number specified according
      ** to JBAS jbas[]. The function returns four orbital
      ** numbers in structure ORB_ID orb;  
      */

static double pn_id_orbit_interchange(int tot_J, int numj_orb, JBAS *jbas, 
                                                     ORB_ID *orb, int length);
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

static void pn_sort_matrix_elements(int num_J, J_POT *j_pot);
      /*
      ** sorts all matrix elements into groups of same total J and
      ** for each group sort after increasing configuration number
      */

static int pn_total_J_comp(const TWO_PART *one, const TWO_PART *two);
      /*
      ** is a utility function for the library function qsort() in order to sort the
      ** two-body matrix elements into groups of increasing value of the total J-value.
      */

static int pn_config_comp(const TWO_PART *one, const TWO_PART *two);
      /*
      ** is a utility function for the library function qsort() in order to sort the two-body
      ** matrix elements for fixed total J into groups of increasing conf_ident number.
      */

static void pn_check_id_matr_elem(int type, int numj_orb, JBAS *jbas, 
                                                  int num_J, J_POT j_pot);
      /*
      ** checks current identical two-particle matrix elements to
      ** see if there are any missing compared to the available
      ** single-particle orbits in the model.
      */

static void pn_read_pn_two_part_matr_elem(char *file_name, int numj_orbZ, JBAS *jbasZ,
                                                int numj_orbN, JBAS *jbasN, J_POT *j_pot);
      /*
      ** reads angular momentum coupled two-particle effective matrix elements
      ** for proton/neutron case.  Each matrix element is given an identity number.
      **  Then they are  sorted into groups of  total_J 
      ** in a monotonic increasing sequence and finally stored in structure J_POT j_pot[].
      ** and store the result in structure J_POT j_pot[].
      */

static ORB_ID pn_pn_matr_elem_identity(int numj_orbZ, JBAS *jbasZ, int numj_orbN,
                                                         JBAS *jbasN, int *int_data);
      /*
      ** converts particle matrix elements identification from
      ** (n,l,j) form into an orbital number specified according
      ** to JBAS jbas[] for protons and neutrons. The function
      ** returns four orbital numbers in structure ORB_ID orb;  
      */

static void pn_check_pn_matr_elem(int numj_orbZ, JBAS *jbasZ, int numj_orbN, 
                                           JBAS *jbasN, int num_J, J_POT j_pot);
      /*
      ** checks current identical two-particle matrix elements to
      ** see if there are any missing compared to the available
      ** single-particle orbits.
      */

static int pn_max_nondiag_id_elem(int num_m, MBAS *mbas);
      /*
      ** calculates and returns the maximum  possible number non-diagonal
      ** two-body matrix elements for identical particles.
      */

static void pn_id_m_pot_diag(int interaction, int num_m, MBAS *mbas, J_POT j_pot, double *id_diag);
      /*
      ** calculates all diagonal m-scheme two-particle matrix elements
      ** <k.l| V |k.l>  for identical particles of an interaction
      ** determined by the parameter 
      **   interaction = 0:  effective interaction 
      **               = 1:  angular momentum interaction
      ** The diagonal matrix elements  for k < l are  stored  in vector: 
      **         table[((2*num_m - k - 3) * k)/2 + l].diag_val
      */

static double  pn_id_m_scheme_diag_veff(int m_orb_a, int m_orb_b, MBAS *mbas,
                                                                    J_POT j_pot);
      /*
      ** calculates and return for identical particles diagonal effective interaction
      ** in m-scheme
      **  <m_orb_a, m_orb_b | Veff | m_orb_a, m_orb_b > =                        
      **    SUM(jtot) (  
      **       clebsch(j_a, m_a, j_b, m_b, jtot) * clebsch(j_a, m_a, j_b, m_b, jtot)         
      **     * <j_orb_a, j_orb_b, jtot| Veff | j_orb_a, j_orb_b, jtot >
      */

static double  pn_id_m_scheme_diag_ang_mom(int k, int l, MBAS *mbas);
      /*
      ** calculates the diagonal two-particle matrix  element in m-scheme of the
      ** angular momentum operator <k,l| J^2 |k,l> for for identical particles.
      */

static void pn_add_id_single_particle_terms(int interaction, int part_num, 
                            int num_m, MBAS *mbas, double *id_diag, int tot_M);
      /*
      ** adds contributions from single-particle terms to the effective
      ** two-particle matrix elements  <k.l | V | k.l> for cases:
      **   interaction = 0:  single-particle energies only  
      **               = 1:  single-particle-terms-from-J**2
      */

static int pn_id_m_pot_nondiag(int interaction, int num_m, MBAS *mbas, J_POT j_pot,
                                                       ID_INT **id_nondiag, ID_INT *id_int);
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
  
static double  pn_id_m_scheme_nondiag_veff(int m_orb_a, int m_orb_b, int m_orb_c,
                                               int m_orb_d, MBAS *mbas, J_POT j_pot);
      /*
      ** calculates and return effective interaction in m-scheme (a <= b and c <= d) 
      **   <m_orb_a, m_orb_b | Veff | m_orb_c, m_orb_d > =                        
      **      SUM(jtot) (  clebsch(j_a, m_a, j_b, m_b, jtot)         
      **                  * clebsch(j_c, m_c, j_d, m_d, jtot)         
      **                  * <j_orb_a, j_orb_b, jtot| Veff | j_orb_c, j_orb_d, jtot >
      */

static double pn_id_m_scheme_nondiag_ang_mom(int i, int j, int k, int l, MBAS *mbas);
      /*
      ** calculates and returns non-diagonal two-particle matrix element
      ** in m-scheme of the angular momentum  operator
      **             <i,j| J**2 |k,l>
      ** for the identical particle case. It is assumed that i < j and k < l.         
      */

static int pn_max_nondiag_pn_elem(int parZ, int num_mZ, MBAS *mbasZ, 
                                                int num_mN, MBAS *mbasN);
      /*
      ** calculates and returns the maximum  possible number of 
      ** proton-neutron  non-diagonal two-body matrix elements.
      */

static void pn_pn_m_pot_diag(int interaction, int num_mZ, MBAS *mbasZ, int num_mN,
                                           MBAS *mbasN, J_POT j_pot, double *pn_diag);
      /*
      ** calculates all diagonal m-scheme two-particle matrix elements
      ** < p(k), n(l) | V | p(k), l(n) > for proton-neutron case with 
      ** an interaction determined by the parameter 
      **   interaction = 0:  effective interaction 
      **               = 1:  angular momentum interaction
      ** The diagonal matrix elements are  stored  in vector: 
      **         table[k * num_mN + l].diag_val
      */

static double  pn_pn_m_scheme_diag_veff(int m_orb_a, int m_orb_b, MBAS *mbasZ,
                                                        MBAS *mbasN, J_POT j_pot);
      /*
      ** calculates and return diagonal proton-neutron effective interaction 
      ** in m-scheme
      **   <m_orb_a(p), m_orb_b(n) | Veff | m_orb_a(p), m_orb_b(n) > =                        
      **    SUM(jtot) ( 
      **     clebsch(j_a, m_a, j_b, m_b, jtot) * clebsch(j_a, m_a, j_b, m_b, jtot)         
      **    * <j_orb_a(p), j_orb_b(n), jtot| Veff | j_orb_a(p), j_orb_b(n), jtot >
      */

static void pn_add_pn_single_particle_terms(int interaction, int type, int part_num,
                                   int num_mZ, int num_mN, MBAS *mbas, double *pn_diag);
      /*
      ** adds single-particle contribution for protons(neutrons) when type = 0(1)
      ** to two-particle matrix elements  <k(p).l(n) | V | k(p).l(n)> for cases:
      **   interaction = 0:  single-particle energies only  
      **               = 1:  single-particle-terms-from-J**2              
      */

static int pn_pn_m_pot_nondiag(int interaction,GR_BASIS *gr_basis, J_POT pot, 
			                PN_INT ***pn_nondiag_int, PN_INT *pn_int);
      /*
      ** calculates all nondiagonal m-scheme proton-neutron two-particle
      ** matrix elements < p(i), n(j) | V | p(k), l(n)> of an interaction
      ** determined by the parameter 
      **   interaction = 0:  effective interaction 
      **               = 1:  angular momentum interaction
      ** The nondiagonal matrix elements  are  stored in a table pointed to by  
      ** table[((k * num_mN+ l].start_nondiag
      ** The function returns number of calculated matrix elements.
      */

static double pn_pn_m_scheme_nondiag_veff(int m_orb_a, int m_orb_b, int m_orb_c, 
                                 int m_orb_d, MBAS *mbasZ, MBAS *mbasN, J_POT j_pot);
      /*
      ** calculates and return proton-neutron effective interaction in 
      ** m-scheme ( i <= j and k <= l) 
      **   < m_orb_a(p) m_orb_b(n) | Veff | m_orb_c(p) m_orb_d(n) > =                        
      **          SUM(J) (  clebsch(j_a m_a j_b m_b J)         
      **                 * clebsch(j_c m_c j_d m_d J)         
      **                 * phase( ji + jk -jj - jl)           
      **                 * < j_orb_a(p) j_orb_b(n) J| Veff | j_orb_c(p) j_orb_d(n) J > 
      */   
 
static double  pn_pn_m_scheme_nondiag_ang_mom(int m_orb_a, int m_orb_b, 
                        int m_orb_c, int m_orb_d, MBAS *mbasZ, MBAS *mbasN);
      /*
      ** calculates and return  non-diagonal two-particle matrix
      ** element in m-scheme of angular momentum  operator                                    
      **             <i,j| J**2 |k,l>                
      ** for proton-neutron case.
      */
                /**** End: function declarations ****/


               /**** The function definitions  ****/ 

     /*
     ** The function 
     **          proton_neutron_effective_interaction()
     **    1. read angular momentum coupled effective
     **       two-particle matrix elements for from file 
     **    2. calculate and store all m-scheme effective
     **       two-particle matrix elements
     **    3. calculate and store all m-scheme matrix 
     **       elements of J**2
     */

void proton_neutron_effective_interaction(INPUT *input, PN_SHELL *model, GR_BASIS *gr_basis)
{

   if(model->part_type >= 4)  {                                    /* p-p interaction */
      pn_id_effective_interaction(input->veff_pp, PROTON, model, gr_basis);
      pn_id_ang_mom_interaction(PROTON, model, gr_basis);
   }
   if((model->part_type == 3) || (model->part_type == 5))   {        /* n-n interaction */
      pn_id_effective_interaction(input->veff_nn, NEUTRON, model, gr_basis);
      pn_id_ang_mom_interaction(NEUTRON, model, gr_basis);
   }
   if((model->part_type >= 2)) {                                    /* p-n interaction */
      pn_pn_effective_interaction(input->veff_pn, model, gr_basis);
      pn_pn_ang_mom_interaction(model, gr_basis);
   }

} /* End: function proton_neutron_effective_interaction() */

     /*
     ** The function 
     **         pn_id_effective_interaction()
     ** reads from file angular momentum coupled two-particle
     ** identical particle matrix elements, transfer the elements 
     ** into m-scheme representation and store the result in tables:
     ** ID_TABLE id_table[type][], type = 0(1) - protons(neutrons)
     */

static void pn_id_effective_interaction(char *file_name, int type, PN_SHELL *model, GR_BASIS *gr_basis)
{
   char     *func = {"pn_id_effective_interaction(): "};
   int      num_elem, num;
   J_POT    j_pot;
   ID_INT   *id_nondiag_elem_bas, *id_nondiag_elem_new, *id_nondiag_elem_old;

           /* memory for  kl-table */

   num_elem = (gr_basis->m_orb[type] * (gr_basis->m_orb[type] - 1)) >> 1;
   gr_basis->op_int.id_diag[type]
            = MALLOC(num_elem, double, func, "id_kl_int_diag[]");  
   gr_basis->op_int.id_nondiag[type]
            = MALLOC(num_elem, ID_INT *, func,"id_kl_int_nondiag_table[]");

   num_elem       = pn_max_nondiag_id_elem(gr_basis->m_orb[type], gr_basis->mbas[type]);
   id_nondiag_elem_bas = MALLOC(num_elem, ID_INT, func,"id_nondiag_int_elements[]");
   
           /* read and store two-particle matrix elements in J-scheme */

   pn_read_id_two_part_matr_elem(file_name, type, model->num_j[type],
                                                 model->jbas[type], &j_pot);

          /* calculate and store diag m-scheme two-particle matrix elements */

   pn_id_m_pot_diag(VEFF, gr_basis->m_orb[type], gr_basis->mbas[type], j_pot, 
                                             gr_basis->op_int.id_diag[type]);

   pn_add_id_single_particle_terms(VEFF, model->part[type], gr_basis->m_orb[type],
                   gr_basis->mbas[type], gr_basis->op_int.id_diag[type], model->MJ);

/**************  test output  ***************

   printf("\n\ndiagonal m_scheme %s effective two-particle matrix elements:\n",
	                 (type == 0) ? "proton-proton" : "neutron-neutron"); 
   table_ptr = gr_basis->op_int.id_table[type];
   for(k = 0; k < gr_basis->m_orb[type]; k++) {
      for(l = k + 1; l < gr_basis->m_orb[type]; l++, table_ptr++) {
         printf("\nk = %2d  l = %2d   diag_value = %8.4f", k, l,table_ptr->diag_val); 
      }
   }

********************  end test output    *************/

       /* calculate and store nondiag m-scheme effective two-particle matrix elements */

   num = pn_id_m_pot_nondiag(VEFF, gr_basis->m_orb[type], gr_basis->mbas[type], j_pot, 
			            gr_basis->op_int.id_nondiag[type], id_nondiag_elem_bas); 


   free(j_pot.start);        /* release temporary memory in pp_pot */
   free(j_pot.two_part);  

   if(num > num_elem) {
      printf("\n\nError in function pn_id_effective_interaction(): ");
      printf("\nfor type = %d matrix elements: memory reservation", type);
      printf(" = %d - calculated = %d",num_elem, num); 
      printf("\n This may not occur!!!!!\n");
      exit(1);
   }

   if(num < num_elem) {           /* reduce size of nondiag[] to new number */
      id_nondiag_elem_old = id_nondiag_elem_bas;
      id_nondiag_elem_new= REALLOC(id_nondiag_elem_old, num, ID_INT, func,
                                                     "new id_nondiag_elements[]");
      if(id_nondiag_elem_new != id_nondiag_elem_old) { 
	printf("\n\nError in function pn_id_effective_interaction():");
	printf("\nFunction REALLOC() does not work properly - move memory location");
	printf("\nof id_nondiag matrix elements to a new place\n");
	exit(1);
      }
   }

/**************  test output  ***************

   printf("\n\nm_scheme %s two-particle matrix elements:\n",
	                 (type == 0) ? "proton-proton" : "neutron-neutron"); 
   table_ptr = gr_basis->op_int.id_table[type];
   for(k = 0; k < gr_basis->m_orb[type]; k++) {
      for(l = k + 1; l < gr_basis->m_orb[type]; l++, table_ptr++) {
         printf("\nk = %2d  l = %2d   diag_value = %8.4f", k, l,table_ptr->diag_val); 
         printf("   adress nondiag = %x ",table_ptr->start_nondiag);
         int_ptr = table_ptr->start_nondiag;
         if(int_ptr == NULL_PTR) continue; 
         do  {
	       printf("\n         ");
               printf("one = %4X  two = %4X  val = %8.4f", int_ptr->one, 
                                                        int_ptr->two, int_ptr->val);
	 } while((++int_ptr)->one);
      }
   } 

********************  end test output    *************/

} /* End: function pn_id_effective_interaction() */

     /*
     ** The function 
     **         pn_id_ang_mom_interaction()
     ** calculates m-scheme matrix elements of J**2 and store the
     ** result in tables: ID_TABLE id_table[type][],
     ** type = 0(1) - protons(neutrons)
     */

static void pn_id_ang_mom_interaction(int type, PN_SHELL *model, GR_BASIS *gr_basis)
{
   char     *func = {"pn_id_ang_mom_interaction(): "};
   int      num_elem, num;
   J_POT    pot;                                      /* dummy */
   ID_INT   *id_nondiag_elem_bas, *id_nondiag_elem_new, *id_nondiag_elem_old;

           /* memory for  kl-table */

   num_elem = (gr_basis->m_orb[type] * (gr_basis->m_orb[type] - 1)) >> 1;
   gr_basis->op_ang.id_diag[type] = MALLOC(num_elem, double, func, "id_kl_ang_table[]");  
   gr_basis->op_ang.id_nondiag[type]
            = MALLOC(num_elem, ID_INT *, func,"kl_ang_nondiag_table[]");

   num_elem       = pn_max_nondiag_id_elem(gr_basis->m_orb[type], gr_basis->mbas[type]);
   id_nondiag_elem_bas = MALLOC(num_elem,ID_INT,func,
                                          "pn_id_effective_nondiag_ang_elements[]");
   
             /* calculate and store diag m-scheme of J**2 */

   pn_id_m_pot_diag(ANG, gr_basis->m_orb[type], gr_basis->mbas[type], pot,
                                      gr_basis->op_ang.id_diag[type]);

   /**************  test output  ***************

   printf("\n\ndiagonal m_scheme %s J**2 matrix elements:\n",
	                 (type == 0) ? "proton-proton" : "neutron-neutron"); 
   table_ptr = gr_basis->id_ang_table[type];
   for(k = 0; k < gr_basis->m_orb[type]; k++) {
      for(l = k + 1; l < gr_basis->m_orb[type]; l++, table_ptr++) {
         printf("\nk = %2d  l = %2d   diag_value = %8.4f", k, l,table_ptr->diag_val); 
      }
   }

********************  end test output    *************/


   pn_add_id_single_particle_terms(ANG, model->part[type], gr_basis->m_orb[type],
                              gr_basis->mbas[type], gr_basis->op_ang.id_diag[type], model->MJ);

          /* calculate and store nondiag m-scheme J**2 matrix elements */

   num = pn_id_m_pot_nondiag(ANG, gr_basis->m_orb[type], gr_basis->mbas[type], pot, 
			      gr_basis->op_ang.id_nondiag[type], id_nondiag_elem_bas); 

   if(num > num_elem) {
      printf("\n\nError in function pn_id_effective_interaction(): ");
      printf("\nfor type = %d J**2 matrix elements: memory reservation", type);
      printf(" = %d - calculated = %d",num_elem, num); 
      printf("\n This may not occur!!!!!\n");
      exit(1);
   }

   if(num < num_elem) {           /* reduce size of nondiag[] to new number */
      id_nondiag_elem_old = id_nondiag_elem_bas;
      id_nondiag_elem_new= REALLOC(id_nondiag_elem_old, num, ID_INT, func,
                               "new id_nondiag_J**2_elements[]");
      if(id_nondiag_elem_new != id_nondiag_elem_old) { 
	printf("\n\nError in function  pn_id_ang_mom_interaction():");
	printf("\nFunction REALLOC() does not work properly - move memory location");
	printf("\nof id_nondiag matrix elements to a new place\n");
	exit(1);
      }
   }

   /**************  test output  ***************

   printf("\n\nm_scheme %s J**2 matrix elements:\n",
	                 (type == 0) ? "proton-proton" : "neutron-neutron"); 
   table_ptr = gr_basis->op_ang.id_table[type];
   for(k = 0; k < gr_basis->m_orb[type]; k++) {
      for(l = k + 1; l < gr_basis->m_orb[type]; l++, table_ptr++) {
         printf("\nk = %2d  l = %2d   diag_value = %8.4f", k, l,table_ptr->diag_val); 
         printf("   adress nondiag = %x",table_ptr->start_nondiag);
         int_ptr = table_ptr->start_nondiag;
         if(int_ptr == NULL_PTR) continue; 
         do  {
	       printf("\n         ");
               printf("one = %4X  two = %4X  val = %8.4f", int_ptr->one, 
                                                        int_ptr->two, int_ptr->val);
	 } while((++int_ptr)->one);
      }
   } 

********************  end test output    *************/

} /* End: function pn_id_ang_mom_interaction() */

     /*
     ** The function 
     **         pn_pn_effective_interaction()
     ** reads from file angular momentum coupled two-particle
     ** proton-neutron matrix elements, transfer the elements into
     ** m-scheme representation and store the result in table:
     **                 PN_TABLE pn_table[].
     */

static void pn_pn_effective_interaction(char *file_name, PN_SHELL *model, GR_BASIS *gr_basis)
{
  char     *func = {"pn_pn_effective_interaction(): "};
  int      k, num_elem, num;
  J_POT    j_pot;
  PN_INT   *pn_nondiag_elem_bas, *pn_nondiag_elem_new, *pn_nondiag_elem_old;

           /* proton-neutron diagonal matrix element kl-table */

  num_elem                 = gr_basis->m_orb[0] * gr_basis->m_orb[1];
  gr_basis->op_int.pn_diag = MALLOC(num_elem, double, func, "pn_kl_int_diag[]");  
 
          /* list of pointers to nondiagomal matrix elements kl-table */

  num = (gr_basis->mbas[0]->m + 1) * (gr_basis->parZ ? 1 : 2);
  gr_basis->op_int.pn_nondiag
            = MALLOC(num, PN_INT **, func,"pn_int_nondiag_table[]");

         /* set of  proton-neutron nondiagonal matrix element kl-table */

  for(k = 0; k< num; k++) {
    gr_basis->op_int.pn_nondiag[k]
      = MALLOC(num_elem, PN_INT *, func,"pn_int_nondiag_table[]");
  }

          /* table of proton-neutron nondiagonal PN_INT matrix elements */

  num_elem = pn_max_nondiag_pn_elem(gr_basis->parZ, gr_basis->m_orb[0], gr_basis->mbas[0],
                                                  gr_basis->m_orb[1], gr_basis->mbas[1]);
  pn_nondiag_elem_bas = MALLOC(num_elem, PN_INT, func, "pn_nondiag_int_elem_bas[]");
   
           /* read and store two-particle matrix elements in J-scheme */

  pn_read_pn_two_part_matr_elem(file_name, model->num_j[0], model->jbas[0],
                                         model->num_j[1], model->jbas[1], &j_pot);

          /* calculate and store diag m-scheme effective two-particle matrix elements */

  pn_pn_m_pot_diag(VEFF, gr_basis->m_orb[0], gr_basis->mbas[0], gr_basis->m_orb[1],
                                  gr_basis->mbas[1], j_pot, gr_basis->op_int.pn_diag);
/**************  test output  ***************

   printf("\n\ndiagonal m_scheme proton_neutron effective two-particle matrix elements:\n");
   table_ptr = gr_basis->pn_int_table;
   for(k = 0; k < gr_basis->m_orb[0]; k++) {
      for(l = 0; l < gr_basis->m_orb[1]; l++, table_ptr++) {
         printf("\nk = %2d  l = %2d   diag_value = %8.4f", k, l,table_ptr->diag_val); 
      }
   }

********************  end test output    *************/

        /*
	** single-particle contributions to effective interation are
        ** calculated together with proton Veff and neutron Veff
        ** except when proton and/or neutron numbers are ONE.
        ** Then contributions are added here!
        */

  if(model->part[0] == 1) {
    pn_add_pn_single_particle_terms(VEFF, PROTON, model->part[1], gr_basis->m_orb[0],
              gr_basis->m_orb[1], gr_basis->mbas[0], gr_basis->op_int.pn_diag);
  }
  if(model->part[1] == 1) {
    pn_add_pn_single_particle_terms(VEFF, NEUTRON, model->part[0], gr_basis->m_orb[0],
                                 gr_basis->m_orb[1], gr_basis->mbas[1], gr_basis->op_int.pn_diag);
  }
          /* calculate and store nondiag m-scheme effective two-particle matrix elements */

  num = pn_pn_m_pot_nondiag(VEFF, gr_basis, j_pot, gr_basis->op_int.pn_nondiag,
                                                           pn_nondiag_elem_bas); 
    
  free((void *)j_pot.two_part);

  free((void *)j_pot.start);        /* release temporary memory in pp_pot */

  if(num > num_elem) {
    printf("\n\nError in function pn_pn_effective_interaction(): ");
    printf("\n p-n matrix elements: memory reservation = %d - calculated = %d",
                                                                   num_elem, num);
    printf("\n This may not occur!!!!!\n");
    exit(1);
  }
  if(num < num_elem) {          /* reduce size of nondiag[] to actual number num. */
    pn_nondiag_elem_old = pn_nondiag_elem_bas;
    pn_nondiag_elem_new = REALLOC(pn_nondiag_elem_old, num, PN_INT, func,
                                             "new ang_mom pn_nondiag_elements[]");
    if(pn_nondiag_elem_new != pn_nondiag_elem_old) { 
      printf("\n\nError in function pn_pn_effective_interaction():");
      printf("\nFunction REALLOC() does not work properly - move memory location");
      printf("\nof pn_nondiag matrix elements to a new place\n");
      exit(1);
    }
  }

/**************  test output  ***************

   printf("\n\nm_scheme proton/neutron effective two-particle matrix elements:\n");

   table_ptr = gr_basis->op_int.pn_table;
   for(k = 0; k < gr_basis->m_orb[0]; k++) {
      for(l = 0; l < gr_basis->m_orb[1]; l++, table_ptr++) {
         printf("\nk = %2d  l = %2d   diag_value = %8.4f", k, l,table_ptr->diag_val); 
         printf("   adress nondiag = %x",table_ptr->start_nondiag);
         int_ptr = table_ptr->start_nondiag;
         if(int_ptr == NULL_PTR) continue; 
         do  {
	       printf("\n         ");
               printf("one_1 = %4X  one_2 = %4X two_1 = %4X  two_2 = %4X   val = %8.4f",
                    int_ptr->one[0], int_ptr->one[1], int_ptr->two[0], int_ptr->two[1],
                                                                       int_ptr->val);
               int_ptr++;
	 } while(int_ptr->one[1]);
      }
   } 

********************  end test output    *************/

} /* End: function pn_pn_effective_interaction() */

     /*
     ** The function 
     **         pn_pn_ang_mom_interaction()
     ** calculates m-scheme proton/neutron matrix elements of J**2
     ** and store the result in tables: PN_TABLE pn_table[].
     */

static void pn_pn_ang_mom_interaction(PN_SHELL *model, GR_BASIS *gr_basis)
{
   char     *func = {"pn_pn_ang_mom_interaction(): "};
   int      k, num_elem, num;
   J_POT    pot;                                             /* dummy */
   PN_INT   *pn_nondiag_elem_bas, *pn_nondiag_elem_new, *pn_nondiag_elem_old;

           /* memory for proton-neutron kl-table */

   num_elem                  = gr_basis->m_orb[0] * gr_basis->m_orb[1];
   gr_basis->op_ang.pn_diag = MALLOC(num_elem, double, func, "pn_kl_ang_table[]");  

   num = (gr_basis->mbas[0]->m + 1) * (gr_basis->parZ ? 1 : 2);
  gr_basis->op_ang.pn_nondiag
            = MALLOC(num, PN_INT **, func,"pn_ang_nondiag_table[]");
  for(k = 0; k< num; k++) {
    gr_basis->op_ang.pn_nondiag[k]
      = MALLOC(num_elem, PN_INT *, func,"pn_int_nondiag_table[]");
  }
  num_elem = pn_max_nondiag_pn_elem(gr_basis->parZ, gr_basis->m_orb[0], gr_basis->mbas[0],
                                                  gr_basis->m_orb[1], gr_basis->mbas[1]);
  pn_nondiag_elem_bas = MALLOC(num_elem, PN_INT, func, "pn_nondiag_ang_elem_bas[]");
   
          /* calculate and store diag m-scheme J**2 matrix elements */

   pn_pn_m_pot_diag(ANG, gr_basis->m_orb[0], gr_basis->mbas[0], gr_basis->m_orb[1],
                                     gr_basis->mbas[1], pot, gr_basis->op_ang.pn_diag);

   /**************  test output  ***************

   printf("\n\ndiagonal m_scheme proton_neutron J**2 two-particle matrix elements:\n");
   table_ptr = gr_basis->pn_ang_table;
   for(k = 0; k < gr_basis->m_orb[0]; k++) {
      for(l = 0; l < ->m_orb[1]; l++, table_ptr++) {
         printf("\nk = %2d  l = %2d   diag_value = %8.4f", k, l,table_ptr->diag_val); 
      }
   }

********************  end test output    *************/

        /*
	** single-particle contributions to J**2  interation are
        ** calculated together with proton Veff and neutron Veff
        ** except when proton and/or neutron numbers are ONE.
        ** Then contributions are added here!
        */

   if(model->part[0] == 1) {
      pn_add_pn_single_particle_terms(ANG, PROTON, model->part[1], gr_basis->m_orb[0],
                       gr_basis->m_orb[1], gr_basis->mbas[0], gr_basis->op_ang.pn_diag);
   }
   if(model->part[1] == 1) {
      pn_add_pn_single_particle_terms(ANG, NEUTRON, model->part[0], gr_basis->m_orb[0],
                                  gr_basis->m_orb[1], gr_basis->mbas[1], gr_basis->op_ang.pn_diag);
   }
          /* calculate and store nondiag m-scheme J**2 two-particle matrix elements */

   num = pn_pn_m_pot_nondiag(ANG, gr_basis, pot, gr_basis->op_ang.pn_nondiag,
                                                            pn_nondiag_elem_bas);

   if(num > num_elem) {
      printf("\n\nError in function pn_pn_ang_mom_interaction(): ");
      printf("\n p-n matrix elements: memory reservation = %d - calculated = %d", 
                                                                   num_elem, num);
      printf("\n This may not occur!!!!!\n");
      exit(1);
   }
   if(num < num_elem) {          /* reduce size of nondiag[] to actual number num. */
      pn_nondiag_elem_old = pn_nondiag_elem_bas;
      pn_nondiag_elem_new = REALLOC(pn_nondiag_elem_old, num, PN_INT, func,
                                   "new ang_mom pn_nondiag_elements[]");
      if(pn_nondiag_elem_new != pn_nondiag_elem_old) { 
	if(pn_nondiag_elem_new != pn_nondiag_elem_old) { 
	  printf("\n\nError in function  pn_pn_ang_mom_interaction():");
	  printf("\nFunction REALLOC() does not work properly - move memory location");
	  printf("\nof pn_nondiag matrix elements to a new place\n");
	  exit(1);
	}
      }
   }
   /**************  test output  ***************

   printf("\n\nm_scheme proton/neutron angular momentum two-particle matrix elements:\n");

   table_ptr = gr_basis->op_ang.pn_table;
   for(k = 0; k < gr_basis->m_orb[0]; k++) {
      for(l = 0; l < gr_basis->m_orb[1]; l++, table_ptr++) {
         printf("\nk = %2d  l = %2d   diag_value = %8.4f", k, l,table_ptr->diag_val); 
         printf("   adress nondiag = %x",table_ptr->start_nondiag);
         int_ptr = table_ptr->start_nondiag;
         if(int_ptr == NULL_PTR) continue; 
         do  {
	       printf("\n         ");
               printf("one_1 = %4X  one_2 = %4X two_1 = %4X  two_2 = %4X   val = %8.4f",
                    int_ptr->one[0], int_ptr->one[1], int_ptr->two[0], 
                                        int_ptr->two[1], int_ptr->val);
               int_ptr++;
	 } while(int_ptr->one[1]);
      }
   } 

********************  end test output    *************/

} /* End: function pn_pn_ang_mom_interaction() */

     /*
     ** The function
     **      pn_read_id_two_part_matr_elem()
     ** reads angular momentum coupled two-particle effective matrix elements
     ** for identical particles, type = 0(1) - proton (neutron). Each matrix element
     ** is given an identity number. Then they are  sorted into groups of  total_J 
     ** in a monotonic increasing sequence and finally stored in structure J_POT j_pot[].
     */

static void pn_read_id_two_part_matr_elem(char *file_name, int type,int numj_orb,
                                                JBAS *jbas, J_POT *j_pot)
{
   char      ch, *func = {"pn_read_id_two_part_matr_elem(): "};
   int       loop, int_data[14], length, max_orb, num_elem, new_num_elem, num_J;
   UL        temp;
   double    dbl_val;
   ORB_ID    orb_ident;
   TWO_PART  *two_part, *ptr;
   FILE      *file_ptr;

/**********  test output  ********
   int       orb_a, orb_b, orb_c, orb_d;
**************   end output   ******/

   length  = sizeof(UL) << 1;       /* number of bits for indentity number */
   max_orb = (UL_ONE << length);     /* maximum number of spherical orbits */

   if(numj_orb > max_orb)   {
      printf("\n\nError in function pn_read_id_two_part_matr_elem():");
      printf("\nToo many spherical %s particle single-particle orbits", 
                                 (type == 0) ? "proton":"neutron");
      printf("\nnumj_orb = %d  -- must not exceed %d!",numj_orb, max_orb);
      printf("\n since an integer is only %d bytes\n", sizeof(int));
      exit(1);
   }

   if( (file_ptr = fopen(file_name,"r")) == NULL) {              /* open input file */
      printf("\n\nError in function pn_read_id_two_part_matr_elem();");
      printf("\nWrong file = %s for input of %s eff.interaction.\n",
                           file_name,(type == 0) ? "proton":"neutron");
      exit(1);
   }

   if((fscanf(file_ptr,"%d%c", &num_elem, &ch) != 2) | (ch !=  '\n')) { /* num matr. elem */
      printf("\n\n Error in function read_two_part_matr_elem(): ");
      printf("\nFirst number in file %s", file_name);
      printf(" - number of matrix element - is not correct\n");
      exit(1);
   } 
            /* memory to store two-particle matrix elements */

   two_part = MALLOC(num_elem, TWO_PART, func,"j_pot.two_part[]");
   
   for(loop = 0, new_num_elem = 0, ptr = two_part; loop < num_elem; loop++) {

            /* read (n1,l1,j1), (n2,l2,j2), (n3,l3,j3), (n4,l4,j4), MT, 2*J */

      if(read_int_number(file_ptr, 14, int_data) == FALSE) {
         printf("\n\nError in function read_two_part_matr_elem();");
         printf("\nSomething is wrong with the (n,l,j) - values");
         printf(" for matrix element number = %d\n",loop);
         exit(1);
      }
            /* read effective matrix element */

      if(new_read_float_number(file_ptr, 1, &dbl_val) == FALSE) {
         printf("\n\nError in function read_two_part_matr_elem();");
         printf("\nWrong value for two-particle matrix element no %d",loop);
         printf("\n read from file %s\n", file_name);
         exit(1);
      }
            /* 
	    ** converts (n,l.j) input data to orbit numbers from
            ** jbas[] and store the result in structure orb_ident
            */ 

      orb_ident = pn_id_matr_elem_identity(numj_orb, jbas, int_data);

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

      dbl_val *= pn_id_orbit_interchange(int_data[13], numj_orb, jbas, &orb_ident, length);

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

      ptr->config     = (orb_ident.left << (2*length)) + orb_ident.right;
      ptr->J          = int_data[13];
      (ptr++)->value  = dbl_val;
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

   num_J = jbas[0].j + 1;
   pn_sort_matrix_elements(num_J, j_pot);

/**************  test output  ***********

   printf("\n\nSpherical %s twobody matrix elements:\n", (type == 0) ? "proton" : "neutron");
   for(loop = 0; loop < j_pot->num; loop++) {
      temp = 0XFF;
      orb_a = (j_pot->two_part[loop].config >> (3 * length)) & temp;
      orb_b = (j_pot->two_part[loop].config >> (2 * length)) & temp;
      orb_c = (j_pot->two_part[loop].config >> (length)) & temp;
      orb_d = (j_pot->two_part[loop].config ) & temp;
      printf("\nloop = %2d", loop);
      printf("  orb_a = %d  orb_b = %d  orb_c = %d  orb_d = %d", orb_a, orb_b, orb_c, orb_d);
      printf("   2*J = %d  value = %8.4f", j_pot->two_part[loop].J, 
                                       j_pot->two_part[loop].value);
   }

*************  end test output  *******/

            /* check that all necessary two-particle matrix elements are stored */
 
   pn_check_id_matr_elem(type, numj_orb, jbas, num_J, *j_pot);

} /* End: function  pn_read_id_two_part_matr_elem() */

    /*
    ** The function                            
    **             read_int_number()                   
    ** reads (number) integers from FILE *in_data stream and
    ** store them in vector[]. If too few integers are read
    ** the function returns  FALSE, otherwise TRUE.
    */

static int read_int_number(FILE *in_data, int number, int *vector)
{
   char   data_string[ONE_LINE], *ptr;
   int    loop, phase, val;

   for(loop = 0; loop < number; loop++) vector[loop] = 0;       /* initialization */

   for(loop = 0; loop < number; loop++)  {

      if(!read_data_string(in_data, data_string))  return FALSE;

        /* convert data_string to an integer number */ 
 
      ptr = data_string;                 /* initialization */
      val = 0;
      phase = +1;
      do  {
         switch(*ptr) {
	      case '+':  phase = +1;
	                 break;
              case '-':  phase = -1;
	                 break;
              default : if(!isdigit((int)(*ptr)))  return FALSE;
                        val *= 10;
                        val += (int) (*ptr - '0');
         }
      } while(*(++ptr));
      vector[loop] = phase * val;
   }  /* end loop through all integer data */
   if(loop < number)  return FALSE;
   else               return TRUE;

}  /* End: function read_int_number()  */


    /*
    ** The function                            
    **             read_float_number()                   
    ** reads (number) float from FILE *in_data stream,
    ** convert them to double and store them in vector[].
    ** If too few doubles are read the function returns 
    ** FALSE, otherwise TRUE.
    */

static int read_float_number(FILE *in_data, int number, double *vector)
{
   char   data_string[ONE_LINE], *ptr;
   int    loop;

   for(loop = 0; loop < number; loop++) vector[loop] = 0.0;       /* initialization */

   for(loop = 0; loop < number; loop++)  {

      if(!read_data_string(in_data, data_string))  return FALSE;

         /* convert data_string to an integer number */ 
  
      ptr = data_string;              /* initialization */
      while(*ptr) {
         if(!(   (*ptr == '+') || (*ptr == '-') 
              || (*ptr == '.') || isdigit(((int)(*ptr))))) return FALSE;
         ptr++;
      }
      vector[loop] = atof(data_string);     /* convert and save float number */

   }  /* end loop through all floating poin numbers */

   if(loop < number)   return FALSE;
   else                 return TRUE;

} /* End: function read_float_number() */

    /*
    ** The function                            
    **             new_read_float_number()                   
    ** reads (number) float from FILE *in_data stream,
    ** convert them to double and store them in vector[].
    ** If too few doubles are read the function returns 
    ** FALSE, otherwise TRUE.
    */

static int new_read_float_number(FILE *in_data, int number, double *vector)
{
   char   data_string[ONE_LINE], *ptr;
   int    loop;

   for(loop = 0; loop < number; loop++) vector[loop] = 0.0;       /* initialization */

   for(loop = 0; loop < number; loop++)  {

      if(!read_data_string(in_data, data_string))  return FALSE;

      vector[loop] = atof(data_string);     /* convert and save float number */

   }  /* end loop through all floating poin numbers */

   if(loop < number)   return FALSE;
   else                 return TRUE;

} /* End: function new_read_float_number() */

    /*
     * The function                            
     *           read_data_string()                   
     * reads one piece of input data as a text_string containing no blanks. 
     * Whitespaces in front of the data_string is removed. Data_string must
     * not contain the character '\'. The function returns FALSE if no 
     * proper data_string is found, otherwise TRUE.
     */

static int read_data_string(FILE *in_data, char *data_string)
{
   int   val;

   *data_string =  0x0;                      /* initialization */

   while((val = getc(in_data)) && isspace(val)); /* remove whitespace before data string */

   if((val == EOF) || (((char) val) == '<')) return FALSE;

   if( (char) val == '\\')  {   /* data_string must not contain the character '\' */
      ungetc(val, in_data);     /* return it to the stream */
      return FALSE;
   }
   *(data_string++) = (char) val;     /* first character in data_string */

           /* read remaining data string */

   while((val = getc(in_data)) && (!isspace(val) && (val != EOF))
                               && ((char) val != '\\') && ((char) val != '<')) 
                                               *(data_string++) = (char) val;

   if( (char) val == '\\')  {   /* data_string must not contain the character '\' */
      ungetc(val, in_data);     /* return it to the stream */
      return FALSE;
   }
   *data_string = 0x0;         /* terminate the string */

   return TRUE;

}  /* End: function read_data_string() */

     /*
     ** The function 
     **         pn_id_matr_elem_identity()
     ** converts particle matrix elements identification from
     ** (n,l,j) form into an orbital number specified according
     ** to JBAS jbas[]. The function returns four orbital
     ** numbers in structure ORB_ID orb;  
     */

static ORB_ID pn_id_matr_elem_identity(int numj_orb, JBAS *jbas, int *int_data)
{
   int        k;
   JBAS       *ptr;
   ORB_ID     orb;
 
   orb.a= -1;                                          /* local initialization */
   orb.b= -1;
   orb.c= -1;
   orb.d= -1;

   for(k = 0, ptr = jbas; k < numj_orb; k++, ptr++)     {
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

} /* End: function pn_id_matr_elem_identity() */

     /*
     ** The function
     **          pn_id_orbit_interchange()
     ** reorder in case of identical particles orb_a, orb_b
     ** orb_c and orb_d such that    
     **        orb_a <= orb_b , orb_c <= orb_d                   
     ** if(orb_a == orb_b) a normalization factor of sqrt2 is produced
     ** to obtain unnormalized two-particle matrix elements. A similar
     ** factor is added if(orb_c == orb_d). Calculate and store in
     ** structure orb_id left and right identity number.
     ** The function returns the corresponding phase_norm factor.               
     */

static double pn_id_orbit_interchange(int tot_J, int numj_orb, JBAS *jbas, 
                                           ORB_ID *orb, int length)
{
   int        phase, temp;
   double     norm, sqr2 = 1.41421356237; /* square root of 2 */;

   phase = +1;      
   norm  = 1.0;                                   /* initialization */

   if(orb->a > orb->b)  {
      temp    = orb->a;
      orb->a  = orb->b;
      orb->b  = temp;
      phase *= -PHASE((jbas[orb->a].j + jbas[orb->b].j - tot_J) >> 1);
   }

   if(orb->c > orb->d)  {
      temp    = orb->c;
      orb->c  = orb->d;
      orb->d  = temp;
      phase *= -PHASE((jbas[orb->c].j + jbas[orb->d].j - tot_J) >> 1);
   }
             /* not normalized matrix elements */

   if(orb->a == orb->b)  norm *= sqr2;
   if(orb->c == orb->d)  norm *= sqr2;

           /* calculate  config number for matrix element */

   orb->left  = (UL)((orb->a << length) + orb->b);
   orb->right = (UL)((orb->c << length) + orb->d);       

   return (phase * norm);

} /* End: function pn_id_orbit_interchange() */

     /* 
     ** The function
     **        pn_sort_matrix_elements()             
     ** sorts all matrix elements into groups of same total J and
     ** for each group sort after increasing configuration number
     */

static void pn_sort_matrix_elements(int num_J, J_POT *j_pot)
{
   char    *func = {"pn_sort_matrix_elements: "};    
   int     jtot, count, temp; 

        /* sort matrix elements into groups after increasing total J-value */

   qsort(j_pot->two_part, (UL) j_pot->num, sizeof(TWO_PART),
                      (int(*)(const void *, const void *)) pn_total_J_comp);

             /* memory for j_pot->start[] */

   j_pot->start = MALLOC(num_J + 1, int, func, "j_pot.start[]");

   for(jtot = 0, count = 0; jtot < num_J; jtot++) {
      j_pot->start[jtot] = count;
      while((count < j_pot->num) && (j_pot->two_part[count].J == (jtot << 1))) count++;

         /* sort the matrix elements with J = tot_J after increasing config */

      if((temp = count - j_pot->start[jtot]) > 0)   {
         qsort(j_pot->two_part + j_pot->start[jtot], (UL)(temp),
            sizeof(TWO_PART),(int(*)(const void *, const void *)) pn_config_comp);
      }
   }
   j_pot->start[jtot] = count;
   
} /* End: function pn_sort_matrix_elements() */

     /*
     ** The function                         
     **        int pn_total_J_comp()                  
     ** is a utility function for the library function qsort() in order to sort the
     ** two-body matrix elements into groups of increasing value of the total J-value.
     */

static int pn_total_J_comp(const TWO_PART *one, const TWO_PART *two)
{
  if(one->J > two->J)       return +1;
  else  if(one->J < two->J) return -1;
  else                        return  0;

} /* End: function pn_total_J_comp() */

     /*
     ** The function                         
     **        int pn_config_comp()                  
     ** is a utility function for the library function qsort() in order
     ** to sort the two-body matrix elements for fixed total J into 
     ** groups of increasing conf_ident number.
     */

static int pn_config_comp(const TWO_PART *one, const TWO_PART *two)
{
  if(one->config > two->config)       return +1;
  else  if(one->config < two->config) return -1;
  else                      return  0;

} /* End: function pn_config_comp() */

     /*
     ** The function 
     **        pn_check_id_matr_elem()
     ** checks current identical two-particle matrix elements to
     ** see if there are any missing compared to the available
     ** single-particle orbits in the model.
     */

static void pn_check_id_matr_elem(int type, int numj_orb, JBAS *jbas, int num_J, J_POT j_pot)
{
   char     *func = {"pn_check_id_matr_elem: "};
   int      ja_orb, jb_orb, j_min, j_max, jtot, temp, control,
            parity,  *states_J_plus, *states_J_minus, *ptr;

   states_J_plus  = CALLOC(num_J, int, func, "id_num_J_plus[]");    /* local memory */
   states_J_minus = CALLOC(num_J, int, func, "id_num_J_minus[]");

   for(ja_orb = 0; ja_orb < numj_orb; ja_orb++) {
      for(jb_orb = ja_orb; jb_orb < numj_orb; jb_orb++) {
         j_min = (abs(jbas[ja_orb].j - jbas[jb_orb].j)) >> 1;
	 j_max = (jbas[ja_orb].j + jbas[jb_orb].j) >> 1;
         parity = ((jbas[ja_orb].l + jbas[jb_orb].l) % 2) ? -1 : +1;
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
         printf("\n\nError in function pn_check_id_matr_elem(): ");
         printf("Not enough two-particle matrix elements input for %s",
	                          ((type) ? "neutrons":"protons"));    
         printf(" with total_J = %d", jtot);
         printf("\nRead number = %d - should have been = %d\n", control, temp);
         exit(1);
      }
   }     

   free(states_J_minus);                         /* release local memory */
   free(states_J_plus);
 
} /* End: function pn_check_id_matr_elem() */

     /*
     ** The function
     **      pn_read_pn_two_part_matr_elem()
     ** reads angular momentum coupled two-particle effective matrix elements
     ** for proton/neutron case.  Each matrix element is given an identity number.
     **  Then they are  sorted into groups of  total_J 
     ** in a monotonic increasing sequence and finally stored in structure J_POT j_pot[].
     ** and store the result in structure J_POT j_pot[].
     */

static void pn_read_pn_two_part_matr_elem(char *file_name, int numj_orbZ, JBAS *jbasZ,
                                          int numj_orbN, JBAS *jbasN, J_POT *j_pot)
{
   char       ch, *func = {"pn_read_pn_two_part_matr_elem(): "};
   int       loop, int_data[14], num_elem, new_num_elem,length, max_orb, num_J;
   UL        temp;
   double    dbl_val;
   ORB_ID    orb_ident;
   TWO_PART  *two_part, *ptr;
   FILE      *file_ptr;

/**********  test output  ********
   int       orb_a, orb_b, orb_c, orb_d;
**************   end output   ******/

   length  = sizeof(UL) << 1;       /* number of bits for orbital number */
   max_orb = (UL_ONE << length);     /* maximum number of spherical orbits */

   if(numj_orbZ > max_orb)   {
      printf("\n\nError in function pn_read_pn_two_part_matr_elem():");
      printf("\nToo many spherical proton  particle single-particle orbits"); 
      printf("\nnumj_orb = %d  -- must not exceed %d!",numj_orbZ, max_orb);
      printf("\n since an integer is only %d bytes\n", sizeof(int));
      exit(1);
   }
   if(numj_orbN > max_orb)   {
      printf("\n\nError in function pn_read_pn_two_part_matr_elem():");
      printf("\nToo many spherical neutron single-particle orbits");
      printf("\nnumj_orb = %d  -- must not exceed %d!",numj_orbN, max_orb);
      printf("\n since an integer is only %d bytes\n", sizeof(int));
      exit(1);
   }
  
   if( (file_ptr = fopen(file_name,"r")) == NULL) {              /* open input file */
      printf("\n\nError in function pn_read_pn_two_part_matr_elem();");
      printf("\nWrong file = %s for input of eff.interaction.\n",file_name);
      exit(1);
   }
   
   
   if((fscanf(file_ptr,"%d%c", &num_elem, &ch) != 2) | (ch !=  '\n')) { /* num of elem */
   printf("\n\n Error in function pn_read_pn_two_part_matr_elem(): ");
    printf("\nFirst number in file %s", file_name);
     printf(" - number of matrix element - is not correct\n");
       exit(1);
   }

           /* memory to store two-particle matrix elements */

   two_part = MALLOC(num_elem, TWO_PART, func,"j_pot.two_part[]");
   
   for(loop = 0, new_num_elem = 0, ptr = two_part; loop < num_elem; loop++) {

           /* read (n1,l1,j1), (n2,l2,j2), (n3,l3,j3), (n4,l4,j4), MT, 2*J */

      if(read_int_number(file_ptr, 14, int_data) == FALSE) {
         printf("\n\nError in function read_two_part_matr_elem();");
         printf("\nSomething is wrong with the (n,l,j) - values");
         printf(" for matrix element number = %d\n",loop);
         exit(1);
      }
           /* read effective matrix element */

      if(new_read_float_number(file_ptr, 1, &dbl_val) == FALSE) {
         printf("\n\nError in function read_two_part_matr_elem();");
         printf("\nWrong value for two-particle matrix element no %d",loop);
         printf("\n read from file %s\n", file_name);
         exit(1);
      }
            /*
	    ** converts (n,l.j) input data to orbit numbers from
            ** jbas[] and store the result in structure orb_ident
            */ 

      orb_ident = pn_pn_matr_elem_identity(numj_orbZ, jbasZ, numj_orbN, jbasN, int_data);
           /*
	    **  Only matrix elements relevant for the present
            **  set of single-particle orbits are saved
            */

      if(  (orb_ident.a == -1) || (orb_ident.b == -1) 
         ||(orb_ident.c == -1) || (orb_ident.d == -1)) continue;
 
           /* calculate  config number for matrix element */

      orb_ident.left  = (UL)((orb_ident.a << length) + orb_ident.b);
      orb_ident.right = (UL)((orb_ident.c << length) + orb_ident.d);       

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

      ptr->config     = (orb_ident.left << (2*length)) + orb_ident.right;
      ptr->J          = int_data[13];
      (ptr++)->value  = dbl_val;          // reduce pn matrix elem
      new_num_elem++;

   } /* end loop through all input two-particle matrix elements */

   fclose(file_ptr);

   if(new_num_elem < num_elem) {           /* reduce size of j_pot.two_part[] */
      two_part = REALLOC(two_part, new_num_elem, TWO_PART, func,"realloc j_pot.two_part[]");
   }
   j_pot->num      = new_num_elem;                        /* save number of matrix element */
   j_pot->two_part = two_part;                            /* and the corresponding pointer */

         /* 
         ** All matrix elements are sorted into groups of same total J
         ** and for each group after increasing conf_ident.
         */

   num_J = (jbasZ[0].j + jbasN[0].j + 2) >> 1;
   pn_sort_matrix_elements(num_J, j_pot);


/**************  test output  ***********

   printf("\n\nSpherical proton/neutron twobody matrix elements:\n");
   for(loop = 0; loop < j_pot->num; loop++) {
      temp = 0XFF;
      orb_a = (j_pot->two_part[loop].config >> (3 * length)) & temp;
      orb_b = (j_pot->two_part[loop].config >> (2 * length)) & temp;
      orb_c = (j_pot->two_part[loop].config >> (length)) & temp;
      orb_d = (j_pot->two_part[loop].config ) & temp;
      printf("\nloop = %2d", loop);
      printf("  orb_a = %d  orb_b = %d  orb_c = %d  orb_d = %d",
                                            orb_a, orb_b, orb_c, orb_d);
      printf("   2*J = %d  value = %8.4f", j_pot->two_part[loop].J,
                                                 j_pot->two_part[loop].value);
  }
*************  end test output  *******/

         /* check that all necessary two-particle matrix elements are stored */
 
   pn_check_pn_matr_elem(numj_orbZ, jbasZ, numj_orbN, jbasN, num_J, *j_pot);

} /* End: function  pn_read_pn_two_part_matr_elem() */

     /*
     ** The function 
     **         pn_pn_matr_elem_identity()
     ** converts particle matrix elements identification from
     ** (n,l,j) form into an orbital number specified according
     ** to JBAS jbas[] for protons and neutrons. The function
     ** returns four orbital numbers in structure ORB_ID orb;  
     */

static ORB_ID pn_pn_matr_elem_identity(int numj_orbZ, JBAS *jbasZ, int numj_orbN,
                                                 JBAS *jbasN, int *int_data)
{
   int        k;
   JBAS       *ptr;
   ORB_ID     orb;
 
   orb.a= -1;                                          /* local initialization */
   orb.b= -1;
   orb.c= -1;
   orb.d= -1;

   for(k = 0, ptr = jbasZ; k < numj_orbZ; k++, ptr++)     {
      if(   (ptr->osc == int_data[0]) && (ptr->l == int_data[1]) 
         && (ptr->j == int_data[2]))  orb.a = k;
      if(   (ptr->osc == int_data[6]) && (ptr->l == int_data[7]) 
         && (ptr->j == int_data[8]))  orb.c = k;
   } 
   for(k = 0, ptr = jbasN; k < numj_orbN; k++, ptr++)     {
      if(   (ptr->osc == int_data[3]) && (ptr->l == int_data[4]) 
         && (ptr->j == int_data[5]))  orb.b = k;
      if(   (ptr->osc == int_data[9]) && (ptr->l == int_data[10]) 
         && (ptr->j == int_data[11])) orb.d = k;
   }

   return orb;

} /* End: function pn_id_matr_elem_identity() */

     /*
     ** The function 
     **        pn_check_pn_matr_elem()
     ** checks current identical two-particle matrix elements to
     ** see if there are any missing compared to the available
     ** single-particle orbits.
     */

static void pn_check_pn_matr_elem(int numj_orbZ, JBAS *jbasZ, int numj_orbN,
                                         JBAS *jbasN, int num_J, J_POT j_pot)
{
   char     *func = {"pn_check_pn_matr_elem: "};
   int      ja_orb, jb_orb, j_min, j_max, jtot, temp, control,
            parity,  *states_J_plus, *states_J_minus, *ptr;

   states_J_plus  = CALLOC(num_J, int, func, "pn_num_J_plus[]");    /* local memory */
   states_J_minus = CALLOC(num_J, int, func, "pn_num_J_minus[]");

   for(ja_orb = 0; ja_orb < numj_orbZ; ja_orb++) {
      for(jb_orb = 0; jb_orb < numj_orbN; jb_orb++) {
         j_min = (abs(jbasZ[ja_orb].j - jbasN[jb_orb].j)) >> 1;
	 j_max = (jbasZ[ja_orb].j + jbasN[jb_orb].j) >> 1;
         parity = ((jbasZ[ja_orb].l + jbasN[jb_orb].l) % 2) ? -1 : +1;
         ptr    = (parity == +1) ? states_J_plus : states_J_minus;
         for(jtot = j_min; jtot <= j_max; jtot++)  ptr[jtot]++;

      } /* end orb_b */
   } /* end orb_a */

   for(jtot = 0; jtot < num_J; jtot++)  {
      control = j_pot.start[jtot + 1] - j_pot.start[jtot];
      temp    = (states_J_plus[jtot] * (states_J_plus[jtot] + 1)) >> 1;
      temp   += (states_J_minus[jtot] * (states_J_minus[jtot] + 1)) >> 1;

      if(control != temp) {
         printf("\n\nError in function pn_check_pn_matr_elem(): ");
         printf("Not enough two-particle proton/neutron matrix input elements");
         printf(" for total_J = %d", jtot);
         printf("\nRead number = %d - should have been = %d\n", control, temp);
         exit(1);
      }
   }     

} /* End: function pn_check_pn_matr_elem() */

     /*
     ** The function
     **        pn_max_nondiag_id_elem()
     ** calculates and returns the maximum  possible number non-diagonal
     ** two-body matrix elements for identical particles.
     */

static int pn_max_nondiag_id_elem(int num_m, MBAS *mbas)
{
   int i, j, k, l, limit, m_kl, parity_kl, count;

   count = 0;                     /* initialization */
   limit = num_m - 1;

   for(k = 0; k < limit; k++)  {         /* first index right  */
      for(l = k+1; l <= limit; l++) {     /* second index right */
         m_kl = mbas[k].m + mbas[l].m;        /* m value and parity */
         parity_kl = mbas[k].par * mbas[l].par;
         for(j = l; j <= limit; j++) {                   /* second index left */
            for(i = (j == l) ? k + 1: 0; i < j; i++){     /* first index left */   
               if(   (m_kl == (mbas[i].m + mbas[j].m)) /* test m value and parity */
                  && (parity_kl == (mbas[i].par * mbas[j].par))) count++;
            } /* end of index i */
         } /* end of index j */
         count++; /* add space for an extra ZERO */
      } /* end of index l */
   } /* end of index k */
  return count;
} /* End: function pn_max_id_nondiag_elem() */

     /*
     ** The function                                          
     **              pn_id_m_pot_diag(...)                    
     ** calculates all diagonal m-scheme two-particle matrix elements
     ** <k.l| V |k.l>  for identical particles of an interaction
     ** determined by the parameter 
     **   interaction = 0:  effective interaction 
     **               = 1:  angular momentum interaction
     ** The diagonal matrix elements  for k < l are  stored  in vector: 
     **         table[((2*num_m - k - 3) * k)/2 + l].diag_val
     */

static void pn_id_m_pot_diag(int interaction, int num_m, MBAS *mbas,
                                           J_POT j_pot, double *id_diag)
{
   int          k, l, limit;

   limit = num_m - 1;                                     /* initialization */

   for(k = 0; k < limit; k++) {
      for(l = k+1; l <= limit; l++, id_diag++) {
	switch(interaction) { 

              case 0: *id_diag                          /* eff-int */
                           = pn_id_m_scheme_diag_veff( k,l, mbas, j_pot);
                      break;
              case 1: *id_diag                    /* ang-mom. interaction */
                           = pn_id_m_scheme_diag_ang_mom(k,l, mbas);
                      break;
              default:  printf("\nError in function pn_id_m_pot_diag(): ");
                        printf("\nWrong interaction = %d\n", interaction);
                        exit(1);
                        break;
        } /* end switch */

      }   /* end l index */
   }  /* end k index */

} /* End: function pn_id_m_pot_diag() */

     /*
     ** The function                                  
     **               pn_id_m_scheme_diag_veff()                   
     ** calculates and return for identical particles diagonal effective interaction
     ** in m-scheme
     **  <m_orb_a, m_orb_b | Veff | m_orb_a, m_orb_b > =                        
     **    SUM(jtot) (  
     **       clebsch(j_a, m_a, j_b, m_b, jtot) * clebsch(j_a, m_a, j_b, m_b, jtot)         
     **     * <j_orb_a, j_orb_b, jtot| Veff | j_orb_a, j_orb_b, jtot >
     */

static double  pn_id_m_scheme_diag_veff(int m_orb_a, int m_orb_b, MBAS *mbas, J_POT j_pot)
{
   int       j_a, m_a, j_b, m_b, j_min, j_max, jtot;
   UL        j_orb_a, j_orb_b, temp, length, config, high, low, mid;
   double    veff_m, cleb;

   veff_m = 0.0;                                /* initialization */
   j_orb_a = (UL) mbas[m_orb_a].orb;                        /* orbit numbers */
   j_orb_b = (UL) mbas[m_orb_b].orb;
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
   length = sizeof(UL) << 1;
   config = (j_orb_a << (3 * length)) + (j_orb_b << (2 * length)) 
            + (j_orb_a << length) + j_orb_b;
   temp = (j_orb_a == j_orb_b) ? TRUE : FALSE;
   for(jtot = j_min; jtot <= j_max; jtot++)  {
      if(temp && (jtot % 2)) continue;  /* no contribution for |(jj)J = odd> */
      low     = j_pot.start[jtot];
      high    = j_pot.start[jtot + 1];
      while(1)   {                         /*  binary search for configuration */
         mid = (low + high) >> 1;
         if(config < j_pot.two_part[mid].config)      high = mid -1;
         else if(config >j_pot.two_part[mid].config)  low = mid + 1;
         else                                         break;
      }    /* end of search for configuration identifier  */
      cleb =  clebsch_gordan(j_a, j_b, jtot << 1,m_a, m_b);
      veff_m += cleb *cleb * j_pot.two_part[mid].value;
   }  /* end loop over all jtot */
   return veff_m;

} /* End: function m_scheme_diag_veff() */

     /*
     ** The function                                
     **          pn_id_m_scheme_diag_ang_mom()            
     ** calculates the diagonal two-particle matrix  element in m-scheme of the
     ** angular momentum operator <k,l| J^2 |k,l> for for identical particles.
     */

static double  pn_id_m_scheme_diag_ang_mom(int k, int l, MBAS *mbas)
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

} /* End: function pn_id_m_scheme_diag_ang_mom() */

     /*
     ** The function
     **       pn_add_id_single_particle_terms()
     ** adds contributions from single-particle terms to the effective
     ** two-particle matrix elements  <k.l | V | k.l> for cases:
     **   interaction = 0:  single-particle energies only  
     **               = 1:  single-particle-terms-from-J**2
     */

static void pn_add_id_single_particle_terms(int interaction, int part_num,
                         int num_m, MBAS *mbas, double *id_diag, int tot_M)
{
   int         k, l, limit;
   double      factor, num;

   limit  = num_m - 1;  /* max. k - values  */
   factor = (tot_M * (tot_M + 2))/(2.0 * part_num * (part_num - 1));
   num    = 1.0 /((double) (part_num - 1));
   for(k = 0; k < limit; k++) { 
      for(l = k+1; l <= limit; l++, id_diag++) {   
         switch(interaction)  { 
                case 0: *id_diag += (mbas[k].e + mbas[l].e) * num;
                        break;
	        case 1: *id_diag +=   0.25 * (double) (mbas[k].j * (mbas[k].j + 2)
                                       + mbas[l].j * (mbas[l].j + 2)) * num;
		        break;
                default:  printf("\nError in function  pn_add_id_single_particle_terms(): ");
                          printf("\nWrong interaction type = %d\n", interaction);
                          exit(1);
                         break;
         } /* end switch */
      } /* end l-loop */
   } /* end k-loop */

} /* End: function  pn_add_id_single_particle_terms() */

     /*
     ** The function                                          
     **              pn_id_m_pot_nondiag(,...)                    
     ** calculates all nondiagonal m-scheme two-particle matrix elements
     ** <i,j| V |k.l>  for identical particles of an interaction
     ** determined by the parameter 
     **   interaction = 0:  effective interaction 
     **               = 1:  angular momentum interaction
     ** The nondiagonal matrix elements  are  stored in a table pointed to by  
     ** table[((2*num_m - k - 3) * k)/2 + l].start_nondiag
     ** The function returns number of calculated matrix elements.
     */
     
static int pn_id_m_pot_nondiag(int interaction, int num_m, MBAS *mbas,
                           J_POT j_pot, ID_INT **id_nondiag, ID_INT *id_int)
{
  int      i, j, k, l, m_kl, parity_kl, limit, count, number;
  ID_INT   *nondiag;

/************** test output   *********

   printf("\n\nnondiagonal identical particle m-scheme matrix elements: interaction = %d\n",
                                                                               interaction); 
************ end test output    ****/

  limit   = num_m - 1;                           /* initialization */
  nondiag = id_int;

  number = 0;
  for(k = 0; k < limit; k++) {                            /* first right */ 
    for(l = k+1; l <= limit; l++, id_nondiag++)  {                       /* second right */
      m_kl = mbas[k].m + mbas[l].m;
      parity_kl = mbas[k].par * mbas[l].par;

      *id_nondiag = nondiag;                /* store pointer */

	 /******************
         table->start_nondiag = nondiag;
	 *******************/

/************** test output   *********

   printf("\n k = %2d  l = %2d  start address = %X\n", k, l, table->start_nondiag);

************ end test output    ****/


      for(count = 0,j = l; j <= limit; j++)  {                    /* second left*/
	for(i = (j == l) ? k + 1: 0; i < j; i++){      /* first left  */
	  if(   (m_kl == mbas[i].m + mbas[j].m) /* test mj and parity */
		&& (parity_kl == mbas[i].par * mbas[j].par))  {
	
	    switch(interaction) {
		      case 0: nondiag->val          /* effective interaction */
				      = pn_id_m_scheme_nondiag_veff(i,j,k,l, mbas, j_pot);
			      break;
		      case 1: nondiag->val          /* angular momentum interaction */
			          = pn_id_m_scheme_nondiag_ang_mom(i,j,k,l, mbas);
		               break;
                      default: printf("\nError in function pn_id_m_pot_nondiag(). ");
                               printf(" Wrong interaction = %d\n",interaction);
                               exit(1);
	    } /* end switch */

	    if(fabs(nondiag->val) > MATRIX_LIMIT) {
	      nondiag->one    = UL_ONE<<j ^ UL_ONE<<i;
	      nondiag->two    = (UL_ONE<<j) - (UL_ONE<<(i+1));
	      count++;
	      number++;  

/**************** test output  **********
printf("\n  number = %3d: <i = %d j = %d|V_nn|k = %d l = %d>",  number, i, j, k, l);
printf("\ncount = %3d  nondiag.one = %4X   nondiag.two = %4X   value = %8.4f",
			    count, nondiag->one, nondiag->two, nondiag->val); 
*********** end test output   *****/

	      nondiag++;
	    }
	  } /* end mj and parity test */
		     	 	    
	}  /* end of index i, second particle left */
      } /* end of index j, first particle left */
      if(!count)  {                  /* no matrix elements for present (k,l) */
	*id_nondiag = NULL_PTR;
      }
      else   {                      /* count matrix elements for present (k,l) */
	number++; 
	(nondiag++)->one = UL_ZERO;    /* terminate each (ij)-set with a ZERO */
      }
    } /* end of index l, second right */
  } /* end of  index k, first right */

  return number;     /* return the total number of elements in nondiag[] */

} /* End: function pn_id_m_pot_nondiag() */

     /*
     ** The function                                  
     **               pn_id_m_scheme_nondiag_veff()                   
     ** calculates and return effective interaction in m-scheme (a <= b and c <= d) 
     **   <m_orb_a, m_orb_b | Veff | m_orb_c, m_orb_d > =                        
     **      SUM(jtot) (  clebsch(j_a, m_a, j_b, m_b, jtot)         
     **                  * clebsch(j_c, m_c, j_d, m_d, jtot)         
     **                  * <j_orb_a, j_orb_b, jtot| Veff | j_orb_c, j_orb_d, jtot >
     */

static double  pn_id_m_scheme_nondiag_veff(int m_orb_a, int m_orb_b, int m_orb_c,
                                             int m_orb_d, MBAS *mbas, J_POT j_pot)
{
   int       j_a, m_a, j_b, m_b, j_c, m_c, j_d, m_d,
             phase, j_min, j_max, jtot, type, length;
   UL        j_orb_a, j_orb_b, j_orb_c, j_orb_d, conf_bra,
             temp, conf_ket, config, high, low, mid;
   double    veff_m;

   veff_m = 0.0;                                /* initialization */
   
   j_orb_a = (UL)mbas[m_orb_a].orb;                        /* orbit numbers */
   j_orb_b = (UL)mbas[m_orb_b].orb;
   j_orb_c = (UL)mbas[m_orb_c].orb;
   j_orb_d = (UL)mbas[m_orb_d].orb;  

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

   if(j_min > j_max)  return veff_m ;

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
   length   = sizeof(UL) << 1;
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
         while(1)   {                         /*  binary search for configuration */
            mid = (low + high) >> 1;
            if(config < j_pot.two_part[mid].config)      high = mid -1;
            else if(config >j_pot.two_part[mid].config)  low = mid + 1;
            else                                             break;
         }    /* end of search for configuration identifier  */

         veff_m += phase * clebsch_gordan(j_a, j_b, jtot << 1,m_a, m_b)
                         * clebsch_gordan(j_c, j_d, jtot << 1,m_c, m_d)
                         * j_pot.two_part[mid].value;
      } /* end if - loop */

      if(type) phase = -phase;     /* for type = 1 change phase for next jtot */

   }  /* end loop over all jtot */

            /* "Whitehead" phase is included */
 
   return (veff_m * mbas[m_orb_a].phase * mbas[m_orb_b].phase
                  * mbas[m_orb_c].phase * mbas[m_orb_d].phase);

} /* End: function pn_id_m_scheme_nondiag_veff() */

     /*
     ** The function                                
     **       pn_id_m_scheme_nondiag_ang_mom()            
     ** calculates and returns non-diagonal two-particle matrix element
     ** in m-scheme of the angular momentum  operator
     **             <i,j| J**2 |k,l>
     ** for the identical particle case. It is assumed that i < j and k < l.         
     */

static double pn_id_m_scheme_nondiag_ang_mom(int i, int j, int k, int l, MBAS *mbas)
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

               /* "Whitehead" phase is included */

   phase =   mbas[i].phase * mbas[j].phase
           * mbas[k].phase * mbas[l].phase;

   if(   ((n_i == n_k) && (l_i == l_k) && (j_i == j_k)) 
      && ((n_j == n_l) && (l_j == l_l) && (j_j == j_l)))  {

      /* (n_i,j_i ) = ( n_k,j_k) & (n_j,j_j) = (n_l,j_l) */

      if((m_i == (m_k - 2)) && (m_j == (m_l + 2)))  { /* m_i = m_k - 1  m_j = m_l + 1 */
         value = (j_k + m_k)*(j_k - m_k + 2)*(j_l - m_l)*(j_l + m_l + 2);
      }
      else if((m_i == (m_k + 2)) && (m_j == (m_l - 2))) { /* m_i = m_k + 1  m_j = m_l - 1 */
         value = (j_k - m_k)*(j_k + m_k + 2)*(j_l + m_l)*(j_l - m_l + 2);
      }
   }

   if(   ((n_i == n_l) && (l_i == l_l) && (j_i == j_l)) 
      && ((n_j == n_k) && (l_j == l_k) && (j_j == j_k)))  {

      /** (n_i,j_i ) = ( n_l,j_l) & (n_j,j_j) = (n_k,j_k) **/

      if((m_i == (m_l + 2)) && (m_j == (m_k - 2)))  { /* m_i = m_l + 1  m_j = m_k - 1 */
         value = (j_k + m_k)*(j_k - m_k + 2)*(j_l - m_l)*(j_l + m_l + 2);
         phase *= -1;
      }
      if((m_i == (m_l - 2)) && (m_j == (m_k + 2)))  { /* m_i = m_l - 1  m_j = m_k + 1 */
         value = (j_k - m_k)*(j_k + m_k + 2)*(j_l + m_l)*(j_l - m_l + 2);
         phase *= -1;
      }
   }
   if(abs(value))      return (0.25 * phase * sqrt((double) value));
   else                return 0.0;

} /* End: function pn_id_m_scheme_nondiag_ang_mom() */

     /*
     ** The function
     **        pn_max_nondiag_pn_elem()
     ** calculates and returns the maximum  possible number of 
     ** proton-neutron  non-diagonal two-body matrix elements.
     */

static int pn_max_nondiag_pn_elem(int parZ, int num_mZ, MBAS *mbasZ, int num_mN, MBAS *mbasN)
{
  register int   total_num, i, j, k, l, m_k, m_l,par_k, par_l,
                 max_del_mp, del_m, del_par, par_p, limit, count, number;
    
/* Calculate the number of non-diagonal  matrix elements with i(p) >= k(p) and j(n) >= l(n),*/

  max_del_mp = mbasZ[0].m;

  limit = parZ ? 1 : 2;    /* number of parity change in matr.elem. */

  total_num = 0;                   /* initialization */

  for(k = 0; k < num_mZ; k++) {            /* right proton particle */
    m_k   = mbasZ[k].m;
    par_k = mbasZ[k].par; 
    for(l = 0; l < num_mN; l++) {          /* right neutron particle */
      m_l   = mbasN[l].m;
      par_l = mbasN[l].par;
         
      number = 0;              /* counter for each (kl) - group */
      
      for(del_m = 0; del_m <= max_del_mp; del_m++)  {
         
         for(del_par = 0, par_p = +1; del_par < limit; del_par++, par_p = -1)   {

            count = 0; /* element counter for each block with given del_par */ 
         
            total_num++;    /* header element */
            
            for(i = ((del_m == 0) && (del_par == 1)) ? 0 : k; 
                                                 i < num_mZ; i++)   {
               for( ;  !(  (del_m == ((m_k   - mbasZ[i].m) / 2))
                         &&(par_p == (par_k * mbasZ[i].par)))      
                     && ( i < num_mZ); i++);
               if(i == num_mZ) break;
               for(j = (i == k) ? l + 1: 0; j < num_mN; j++)   {
                  for( ;   !(  (del_m == ((mbasN[j].m - m_l) / 2)) 
                            &&(par_p == (mbasN[j].par * par_l)))
                        && ( j < num_mN); j++);
                  if(j == num_mN) break; 
                  
                  count++;
                  number++;                 
                  total_num++;   /* matrix element */
               
               } /* end loop through j-neutrons */
            } /* end loop through i-protons */  
            if(!count)  {
                total_num--;  /*no matrix element in the present group */ 
	    }
         } /* end del_par loop */

      } /* end del _m for the two parity blocks of matrix elements.*/
         
      if(number) total_num++; /* extra bottom element */   

    } /* end of index l, neutron particle right */
  } /* end of  index k, proton particle right */

  return total_num;

} /* End: function pn_max_nondiag_pn_elem() */

     /*
     ** The function                                          
     **              pn_pn_m_pot_diag(...)                    
     ** calculates all diagonal m-scheme two-particle matrix elements
     ** < p(k), n(l) | V | p(k), l(n) > for proton-neutron case with 
     ** an interaction determined by the parameter 
     **   interaction = 0:  effective interaction 
     **               = 1:  angular momentum interaction
     ** The diagonal matrix elements are  stored  in vector: 
     **          pn_diag[k * num_mN + l].diag_val
     */

static void pn_pn_m_pot_diag(int interaction, int num_mZ, MBAS *mbasZ,
                      int num_mN, MBAS *mbasN, J_POT j_pot, double *pn_diag)
{
   int          k, l;

   for(k = 0; k < num_mZ; k++) {                                  /* proton right */
      for(l = 0; l < num_mN; l++, pn_diag++) {                      /* neutron right */
         switch(interaction)   {
	      case 0 :  *pn_diag = pn_pn_m_scheme_diag_veff(k,l, mbasZ, mbasN, j_pot);
		        break;
	      case 1:   *pn_diag = 0.5 * ((double) (mbasZ[k].m * mbasN[l].m));
		        break;
              default:  printf("\nError in function pn_pn_m_pot_diag() :");
                        printf(" Wrong interaction = %d\n",interaction);
                        exit(1);
                        break;
         } /* end switch() */
      } /* end index l, neutron right */
   } /* end index k, proton right */
   
} /* End: function pn_pn_m_pot_diag() */

     /*
     ** The function                                  
     **               pn_pn_m_scheme_diag_veff()                   
     ** calculates and return diagonal proton-neutron effective interaction 
     ** in m-scheme
     **   <m_orb_a(p), m_orb_b(n) | Veff | m_orb_a(p), m_orb_b(n) > =                        
     **    SUM(jtot) ( 
     **     clebsch(j_a, m_a, j_b, m_b, jtot) * clebsch(j_a, m_a, j_b, m_b, jtot)         
     **    * <j_orb_a(p), j_orb_b(n), jtot| Veff | j_orb_a(p), j_orb_b(n), jtot >
     */

static double  pn_pn_m_scheme_diag_veff(int m_orb_a, int m_orb_b, MBAS *mbasZ,
                                                      MBAS *mbasN, J_POT j_pot)
{
   int       j_a, m_a, j_b, m_b, j_min, j_max, jtot;
   UL        j_orb_a, j_orb_b, length, config, high, low, mid;
   double    veff_m, cleb;


   veff_m = 0.0;                                                   /* initialization */
   j_orb_a = (UL) mbasZ[m_orb_a].orb;                              /* orbit numbers */
   j_orb_b = (UL) mbasN[m_orb_b].orb;
   j_a = mbasZ[m_orb_a].j;                                         /* j and m values */
   m_a = mbasZ[m_orb_a].m;
   j_b = mbasN[m_orb_b].j;
   m_b = mbasN[m_orb_b].m;
   j_min = (abs(j_a - j_b)) >> 1;
   j_max = (j_a + j_b) >> 1;
   length = sizeof(UL) << 1;
   config =   (j_orb_a << (3 * length)) + (j_orb_b << (2 * length))
            + (j_orb_a << length) + j_orb_b;

   for(jtot = j_min; jtot <= j_max; jtot++)  {
      low  = j_pot.start[jtot];
      high = j_pot.start[jtot + 1];
      while(1)   {                         /*  binary search for configuration */
            mid = (low + high)>>1;
            if(config < j_pot.two_part[mid].config)      high = mid -1;
            else if(config >j_pot.two_part[mid].config)  low = mid + 1;
            else                                         break;
      }    /* end of search for configuration identifier  */

      cleb =  clebsch_gordan(j_a, j_b, jtot << 1,m_a, m_b);
      veff_m += cleb *cleb * j_pot.two_part[mid].value;

   }            /* End of angular mom loop         */

  return veff_m;

} /* End: function pn_pn_m_scheme_diag_veff() */

     /*
     ** The function
     **                 pn_add_pn_single_particle_terms()
     ** adds single-particle contribution for protons(neutrons) when type = 0(1)
     ** to two-particle matrix elements  <k(p).l(n) | V | k(p).l(n)> for cases:
     **   interaction = 0:  single-particle energies only  
     **               = 1:  single-particle-terms-from-J**2              
     */

static void pn_add_pn_single_particle_terms(int interaction, int type, int part_num, 
                                 int num_mZ, int num_mN, MBAS *mbas, double *pn_diag)
{
   int       k, l;
   double    sp_energy, ang_mom, num;

   num = 1.0/((double)(part_num));
   for(k = 0; k < num_mZ; k++) {
      if(type == 0)  {
         sp_energy = mbas[k].e * num;
         ang_mom = (0.25 * mbas[k].j * (mbas[k].j + 2.0))* num;
      }
      for(l = 0; l < num_mN; l++, pn_diag++) {
         if(type == 1)  {
            sp_energy = mbas[l].e * num;
            ang_mom = (0.25 * mbas[l].j * (mbas[l].j + 2.0)) * num;
         }
         switch(interaction)  { 
               case 0: *pn_diag += sp_energy; 
                       break;
	       case 1:
		       *pn_diag += ang_mom;
                       break;
	       default: printf("\nError in function pn_add_pn_single_particle_terms():"); 
		 printf("\nWrong interaction type = %d\n",interaction);
                 exit(1);
                 break;
         } /* end switch() */  
      } /* end l-loop */
   } /* end k-loop */

} /* End: function pn_add_pn_single_particle_terms() */

     /*
     ** The function                                          
     **              pn_pn_m_pot_nondiag(,...)                    
     ** calculates all nondiagonal m-scheme proton-neutron two-particle
     ** matrix elements < p(i), n(j) | V | p(k), l(n)> of an interaction
     ** determined by the parameter 
     **   interaction = 0:  effective interaction 
     **               = 1:  angular momentum interaction
     ** The nondiagonal matrix elements  are  stored in a list
     ** pn_nondiag_elem[] pointed to by  
     ** table[((k * num_mN+ l].start_nondiag
     ** The function returns number of calculated matrix elements.
     */

static int pn_pn_m_pot_nondiag(int interaction, GR_BASIS *gr_basis, J_POT pot, 
                                        PN_INT ***pn_nondiag, PN_INT *pn_nondiag_elem)
{
  int      
           num_mZ, num_mN, i, j, k, l, m_k, m_l,par_k, par_l, 
           del_m, del_par, num_parity, par_p, number, tot_num;
  MBAS
           *mbasZ, *mbasN;
  PN_INT
           **kl_table, *nondiag;

  num_mZ     = gr_basis->m_orb[0];
  num_mN     = gr_basis->m_orb[1];
  num_parity = gr_basis->parZ ? 1 : 2;      /* number of parity change in matr.elem. */
  mbasZ      = gr_basis->mbas[0];  
  mbasN      = gr_basis->mbas[1];

  nondiag = pn_nondiag_elem;
  tot_num = 0;

  for(del_m = 0; del_m <=  mbasZ[0].m; del_m++)  {
    for(del_par = 0, par_p = +1; del_par < num_parity; del_par++, par_p = -1, pn_nondiag++) {

       /* 
       ** for each del_m and del_p: calculate non-diagonal
       **  matrix elements with i(p) >= k(p) and j(n) >= l(n),
       */
    
      kl_table = *pn_nondiag;
      for(k = 0; k < num_mZ; k++) {            /* right proton particle */
	m_k   = mbasZ[k].m;
	par_k = mbasZ[k].par;
	for(l = 0; l < num_mN; l++, kl_table++) {          /* right neutron particle */
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

	      switch(interaction)  {
                         case 0: 
			    nondiag->val = pn_pn_m_scheme_nondiag_veff(i,j,k,l, mbasZ,
                                                                         mbasN, pot);
                            break;
                         case 1: 
			    nondiag->val = pn_pn_m_scheme_nondiag_ang_mom(i,j,k,l,
                                                                         mbasZ, mbasN);
                            break;
                         default: 
                            printf("\nError in function pn_pn_m_pot_nondiag(): ");
                            printf(" Wrong interaction type = %d\n",interaction);
                            exit(1);
	      } /* end switch */
	      if(fabs(nondiag->val) > MATRIX_LIMIT) {
		nondiag->one[0] =  UL_ONE<<i;
		nondiag->one[1] =  UL_ONE<<j;
		nondiag->two[0]
                            = (i == k) ? 0 : (UL_ONE<<MAX(i,k)) - (UL_ONE<<(MIN(i,k) + 1));
		nondiag->two[1] 
                              = (j == l) ? 0 : (UL_ONE<<MAX(j,l)) - (UL_ONE<<(MIN(j,l)+1));
		number++;
		tot_num++;
		nondiag++;
	      }
	    }  /* end loop through j-neutrons */
	  } /* end loop through i-protons */  
	  if(number)   { 
	    nondiag->one[0]     = 0;  /* values to terminate a block of matrix elemnts.*/
	    (nondiag++)->one[1] = 0;
	    tot_num++;
	  }     
	  else   {                /* no matrix elements in the present block.*/
	    *kl_table = NULL_PTR;
	  }
	} /* end of index l, neutron particle right */
      } /* end of  index k, proton particle right */

    } /* end del_par loop */
  } /* end del _m for the two parity blocks of matrix elements.*/

  return tot_num;       /* return total number of matrix elements */

} /* End: function pn_pn_m_pot_nondiag() */

     /*
     ** The function                                  
     **               pn_pn_m_scheme_nondiag_veff()                   
     ** calculates and return proton-neutron effective interaction in 
     ** m-scheme ( i <= j and k <= l) 
     **   < m_orb_a(p) m_orb_b(n) | Veff | m_orb_c(p) m_orb_d(n) > =                        
     **          SUM(J) (  clebsch(j_a m_a j_b m_b J)         
     **                 * clebsch(j_c m_c j_d m_d J)         
     **                 * phase( ji + jk -jj - jl)           
     **                 * < j_orb_a(p) j_orb_b(n) J| Veff | j_orb_c(p) j_orb_d(n) J > 
     */

static double  pn_pn_m_scheme_nondiag_veff(int m_orb_a, int m_orb_b, int m_orb_c,
                                int m_orb_d, MBAS *mbasZ, MBAS *mbasN, J_POT j_pot)
{
   int       j_a, m_a, j_b, m_b, j_c, m_c, j_d, m_d,
             j_min, j_max, jtot, length;
   UL        j_orb_a, j_orb_b, j_orb_c, j_orb_d, conf_bra,
             temp, conf_ket, config, high, low, mid;

  double    veff_m;

   veff_m = 0.0;                                          /* initialization */

   j_orb_a = (UL)mbasZ[m_orb_a].orb;                        /* orbit numbers */
   j_orb_b = (UL)mbasN[m_orb_b].orb;
   j_orb_c = (UL)mbasZ[m_orb_c].orb;
   j_orb_d = (UL)mbasN[m_orb_d].orb;  

   j_a = mbasZ[m_orb_a].j;                 /* j and m values */
   m_a = mbasZ[m_orb_a].m;
   j_b = mbasN[m_orb_b].j;
   m_b = mbasN[m_orb_b].m;
   j_c = mbasZ[m_orb_c].j;
   m_c = mbasZ[m_orb_c].m;
   j_d = mbasN[m_orb_d].j;
   m_d = mbasN[m_orb_d].m;

   j_min = MAX(abs(j_a - j_b), abs(j_c - j_d));  /* max and min values of  total J*2 */  
   j_max = MIN(j_a + j_b,  j_c + j_d );

   if( j_min > j_max )  return veff_m ;

   length   = sizeof(UL) << 1;
   conf_bra = (j_orb_a << length) + j_orb_b;
   conf_ket = (j_orb_c << length) + j_orb_d;

   if (conf_bra > conf_ket) {     /* possible interchange left and right of matr. elem. */
      temp     = conf_bra ;
      conf_bra = conf_ket ;
      conf_ket = temp ;
   }
   config = (conf_bra << (2 * length)) + conf_ket;

   for(jtot = j_min; jtot <= j_max; jtot +=2)  {
      low     = j_pot.start[jtot>>1];
      high    = j_pot.start[(jtot>>1) + 1];
      while(1)   { 
         mid = (low + high)>>1;
         if(config < j_pot.two_part[mid].config)      high = mid -1;
         else if(config >j_pot.two_part[mid].config)  low = mid + 1;
         else                                         break;
      }
      veff_m +=   clebsch_gordan(j_a,j_b,jtot,m_a,m_b)
                * clebsch_gordan(j_c,j_d,jtot,m_c,m_d)
                * j_pot.two_part[mid].value;
   }   
                /* "Whitehead time-reversal" phase is included */

   veff_m *=  mbasZ[m_orb_a].phase * mbasN[m_orb_b].phase
     * mbasZ[m_orb_c].phase * mbasN[m_orb_d].phase;

  return veff_m;

} /* End: function pn_pn_m_scheme_nondiag_veff() */

     /*
     ** The function                                
     **       pn_pn_m_scheme_nondiag_ang_mom_()            
     ** calculates and return  non-diagonal two-particle matrix
     ** element in m-scheme of angular momentum  operator                                    
     **             <i,j| J**2 |k,l>                
     ** for proton-neutron case.
     */

static double  pn_pn_m_scheme_nondiag_ang_mom(int m_orb_a, int m_orb_b,
                        int m_orb_c, int m_orb_d, MBAS *mbasZ, MBAS *mbasN)
{
   int       j_a, m_a, j_b, m_b, j_c, m_c, j_d, m_d, value;

   if(  (mbasZ[m_orb_a].orb != mbasZ[m_orb_c].orb) 
     || (mbasN[m_orb_b].orb != mbasN[m_orb_d].orb)) {
     return D_ZERO;     /* J**2 matrix element is ZERO */
   }

      j_a = mbasZ[m_orb_a].j;      /* m_orb_a and m_orb_c are particle one */
   m_a = mbasZ[m_orb_a].m;
   j_c = mbasZ[m_orb_c].j;
   m_c = mbasZ[m_orb_c].m;

   j_b = mbasN[m_orb_b].j;    /* m_orb_b and m_orb_d are particle one */
   m_b = mbasN[m_orb_b].m;
   j_d = mbasN[m_orb_d].j;
   m_d = mbasN[m_orb_d].m;

   value = 0;

   if((m_a == (m_c - 2)) && (m_b == (m_d + 2)))  {       /* m_a = m_c - 1  m_b = m_d + 1 */
      value = (j_c + m_c) * (j_c - m_c + 2) * (j_d - m_d) * (j_d + m_d + 2);
   }
   else if((m_a == (m_c + 2)) && (m_b == (m_d - 2))) {   /* m_a = m_c + 1  m_b = m_d - 1 */
      value = (j_c - m_c) * (j_c + m_c + 2) * (j_d + m_d) * (j_d - m_d + 2);
   }
   if(value)   {
      return (  0.25 * sqrt((double) value)
            * mbasZ[m_orb_a].phase * mbasZ[m_orb_c].phase /* Whitehead time-reversal phase */
              * mbasN[m_orb_b].phase * mbasN[m_orb_d].phase);
   }  
   else       return D_ZERO;

} /* End: function pn_pn_m_scheme_nondiag_ang_mom() */
