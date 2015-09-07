
/*******************  The module pn_two_part_intJ.c  ******************/

#include  "PAR-pnShellModel.h"

     /*
     ** The entrance function 
     **         pn_two_part_intJ()
     ** reads from file angular momentum coupled two-particle
     ** identical particle matrix elements, transfer the elements
     ** into m-scheme representation and store diagonal matrix
     ** elements in 
     **         MATR_OP op_int->id_diag[] 
     ** and pointers to the nondiag matrix elements in 
     **         MATR_OP op_int->id_nondiag_table[]
     */

// data definitions for two-particle matrix elements

     /* Storage of the spherical j-j coupled effective
     ** matrix elements for identical particles
     **         <n_a,n_b | V | n_c, n_d>
     */
 
/*******************   local data structure  ***********/

   typedef struct  {
          int           J;  // twice total J
          UL       config;  // = n_a << 3*length + n_b << 2*length
			    // + n_c<<length + n_d,
	                    // with length = sizeof(ULL) 
          double value[2];  // value of veff element[0] and 
                            // center_of_mass element [1]
   } TWO_PART;

   typedef struct  {
           int
                     num,    /* total number of matrix elements */
                  *start;    /* start position of given J group */
           TWO_PART
               *two_part;    /* pointer to the matrix elements  */
   } J_POT;

                /**** end data definitions ****/

    /******  local function declarations  ********/

static int pn_max_nondiag_elem(int parZ, SP_BAS *spBasZ, SP_BAS *spBasN);
     /*
     ** calculates and returns the maximum  possible number of 
     ** proton-neutron  non-diagonal two-body matrix elements.
     */

 static void pn_read_two_part_matr_elem(char *file_name, SP_BAS *spBasZ, 
			       SP_BAS *spBasN, int calc_CM,J_POT *j_pot);
     /*
     ** reads angular momentum coupled two-particle effective matrix elements
     ** for proton/neutron case.  Each matrix element is given an identity number.
     **  Then they are  sorted into groups of  total_J 
     ** in a monotonic increasing sequence and finally stored in structure J_POT j_pot[].
     ** and store the result in structure J_POT j_pot[].
     */

static ORB_ID pn_matr_elem_identity(SP_BAS *spBasZ, SP_BAS *spBasN,int *int_data);
     /*
     ** converts particle matrix elements identification from
     ** (n,l,j) form into an orbital number specified according
     ** to JBAS jbas[] for protons and neutrons. The function
     ** returns four orbital numbers in structure ORB_ID orb;  
     */

static void pn_sort_matrix_elements(int num_J, J_POT *j_pot);
     /*
     ** sorts all matrix elements into groups of same total J and
     ** for each group sort after increasing configuration number
     */

static int pn_total_J_comp(const TWO_PART *one, const TWO_PART *two);
     /*     
     ** is a utility function for the library function qsort() 
     ** in order to sort the two-body matrix elements into groups
     ** of increasing value of the total J-value.
     */

static int pn_config_comp(const TWO_PART *one, const TWO_PART *two);
     /*
     ** is a utility function for the library function qsort() in order
     ** to sort the two-body matrix elements for fixed total J into 
     ** groups of increasing conf_ident number.
     */

 static void pn_check_matr_elem(SP_BAS *spBasZ, SP_BAS *spBasN, int num_J,
				   J_POT j_pot); 
     /*
     ** checks current identical two-particle matrix elements to
     ** see if there are any missing compared to the available
     ** single-particle orbits.
     */

static void pn_m_pot_diag(SP_BAS *spBasZ, SP_BAS *spBasN,
			  J_POT j_pot, int calc_CM, MATR_PN *op_int);
     /*
     ** calculates all diagonal m-scheme two-particle matrix elements
     ** < p(k), n(l) | VEFF | p(k), l(n) > for proton-neutron case with 
     ** The diagonal matrix elements are  stored  in vector: 
     **          diag[k * num_mN + l].diag_val
     */

static void  pn_m_scheme_diag_int(int m_orb_a, int m_orb_b, MBAS *mbasZ,
		   MBAS *mbasN, J_POT j_pot, int calc_CM, double *pnDiag);
     /*
     ** calculates and return diagonal proton-neutron effective interaction 
     ** in m-scheme
     **   <m_orb_a(p), m_orb_b(n) | Veff | m_orb_a(p), m_orb_b(n) > =                        
     **    SUM(jtot) ( 
     **     clebsch(j_a, m_a, j_b, m_b, jtot) * clebsch(j_a, m_a, j_b, m_b, jtot)         
     **    * <j_orb_a(p), j_orb_b(n), jtot| Veff | j_orb_a(p), j_orb_b(n), jtot >
     */

static int pn_m_pot_nondiag(int parZ, SP_BAS *spBasZ, SP_BAS *spBasN,
                               int calc_CM,J_POT pot, MATR_PN *op_int);

     /*
     ** calculates all nondiagonal m-scheme proton-neutron two-particle
     ** matrix elements < p(i), n(j) | VEFF | p(k), l(n)>
     ** The nondiagonal matrix elements  are  stored in a list
     ** nondiag_elem[] pointed to by  
     ** table[((k * num_mN+ l].start_nondiag
     ** The function returns number of calculated matrix elements.
     */

static void pn_m_scheme_nondiag_veff(int m_orb_a, int m_orb_b, int m_orb_c,
                                int m_orb_d, MBAS *mbasZ, MBAS *mbasN, 
				     int calc_CM,J_POT j_pot, double *matrElem);
     /*
     ** calculates and return proton-neutron effective interaction in 
     ** m-scheme ( i <= j and k <= l) 
     **   < m_orb_a(p) m_orb_b(n) | Veff | m_orb_c(p) m_orb_d(n) > =                        
     **          SUM(J) (  clebsch(j_a m_a j_b m_b J)         
     **                 * clebsch(j_c m_c j_d m_d J)         
     **                 * phase( ji + jk -jj - jl)           
     **                 * < j_orb_a(p) j_orb_b(n) J| Veff | j_orb_c(p) j_orb_d(n) J > 
     */

static void addPN_1to2SinglePartTermVeff(int type, int partNum,int num_mZ,
                       int num_mN,MBAS *mbas, int calc_CM,MATR_PN *op_int);
     /*
     ** adds for proton-neutron particles  contributions from 
     ** single-particle terms to the effective two-particle 
     ** matrix elements  <k.l |int| k.l>
     */

    /******  end local function declarations  ********/

     /*
     ** The entrance function 
     **         pn_two_part_intJ()
     ** reads from file angular momentum coupled two-particle
     ** proton-neutron matrix elements, transfer the elements
     ** into m-scheme representation and store diagonal matrix
     ** elements in 
     **         MATR_PN op_int->id_diag[] 
     ** and pointers to the nondiag matrix elements in 
     **         MATR_PN op_int->id_nondiag_table[]
     */

void pn_two_part_intJ(char *filename, int parZ, SP_BAS *spBasZ,
                  SP_BAS *spBasN, int calc_CM, MATR_PN *op_int)
{
  char     *func = {" pn_two_part_int(): "};
  int      num_kl_elem, num_kl_tables, k, numCalc;
  J_POT    j_pot;
  PN_INT   *nondiag_bas;

  num_kl_elem = spBasZ->numm_orb * spBasN->numm_orb;

  op_int->pn_diag = MALLOC(num_kl_elem,double, func,
			          "pn_kl_diag[]");  
  if(calc_CM == 1) {
    op_int->pn_CM_diag = MALLOC(num_kl_elem,double, func,
		                    "pn_kl_CM_diag[]");
  }

  // num del_m and del_p kl-tables

  num_kl_tables = (spBasZ->mbas[0].m + 1) * (parZ ? 1 : 2);

  op_int->kl_list = MALLOC(num_kl_tables, PN_INT **, func,
                                      "op_int->kl_list[]");
  for(k = 0; k < num_kl_tables; k++) {
    op_int->kl_list[k]
      = MALLOC(num_kl_elem, PN_INT *,func,"pn_int_kl_table[]");
  }

  op_int->numNondiagElem  = pn_max_nondiag_elem(parZ, spBasZ, spBasN);

  op_int->nondiag_bas = MALLOC(op_int->numNondiagElem, PN_INT, func,
                                        "pn_nondiag_int_elem_bas[]");
  // read and store proton-neutron matrix elements in J-scheme

  pn_read_two_part_matr_elem(filename,spBasZ,spBasN,calc_CM,&j_pot);

    // calculate and store diag m-scheme two-particle matr elem

  pn_m_pot_diag(spBasZ, spBasN, j_pot,calc_CM, op_int);

        /*
	** If number of proton and/or neutron is less
	** than two the single-particle contributions
	** to effective interation are added here 
        */
  
  if(spBasZ->part == 1) { 
    addPN_1to2SinglePartTermVeff(PROTON,spBasN->part,
	spBasZ->numm_orb,spBasN->numm_orb,
             spBasZ->mbas, calc_CM,op_int);
  }
  if(spBasN->part == 1) {
    addPN_1to2SinglePartTermVeff(NEUTRON,spBasZ->part,
	spBasZ->numm_orb,spBasN->numm_orb,
        spBasN->mbas, calc_CM,op_int);
  }
  
  // calculate and store nondiag m-scheme two-particle matr elem

  numCalc = pn_m_pot_nondiag(parZ, spBasZ, spBasN, calc_CM,j_pot, op_int);

  if(numCalc > op_int->numNondiagElem) {
    printf("\n\nRank%d: Error in pn_two_part_intJ(): ",Rank);
    printf("\n pn_matrElem: mem reserve = %d - calc = %d",
	                  op_int->numNondiagElem, numCalc);
    printf("\n This may not occur!!!!!\n");
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }
  if(numCalc == 0) {
    free(op_int->nondiag_bas);
    op_int->nondiag_bas      = (PN_INT *)NULL_PTR;
    op_int->numNondiagElem= 0;
    for(k = num_kl_tables - 1; k >= 0; k--) {
      free(op_int->kl_list[k]);
    }
    free(op_int->kl_list);
  }
  else if(numCalc < op_int->numNondiagElem) {// reduce nondiag[]
    op_int->numNondiagElem = numCalc;
    nondiag_bas = op_int->nondiag_bas;
    op_int->nondiag_bas = REALLOC(nondiag_bas,op_int->numNondiagElem,
                                    PN_INT, func,"new nondiag_bas[]");
    if(op_int->nondiag_bas != nondiag_bas) { 
      printf("\n\nRank%d: Error in pn_two_part_intJ(): ",Rank);
      printf("\nREALLOC() does not work properly - moved mem loc");
      printf("\nof nondiag matrix elements to a new place\n");
      MPI_Abort(MPI_COMM_WORLD,Rank);
    }
  }
  free(j_pot.start); // release temp mem in op_pn
  free(j_pot.two_part);  

  return;

} // End: function pn_two_part_intJ()

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
     **      pn_read_two_part_matr_elem()
     ** reads angular momentum coupled two-particle effective matrix elements
     ** for proton/neutron case.  Each matrix element is given an identity number.
     **  Then they are  sorted into groups of  total_J 
     ** in a monotonic increasing sequence and finally stored in structure J_POT j_pot[].
     ** and store the result in structure J_POT j_pot[].
     */

static void pn_read_two_part_matr_elem(char *file_name, SP_BAS *spBasZ, 
			       SP_BAS *spBasN, int calc_CM,J_POT *j_pot)
{
   char       ch, *func = {"read_two_part_matr_elem(): "};
   int       loop, int_data[14], num_elem, new_num_elem,length, max_orb, num_J;
   UL        temp;
   ORB_ID    orb_ident;
   TWO_PART  *two_part, *ptr;
   FILE      *file_ptr;

   length  = 2*sizeof(UL);   // number of bits for orbital number
   max_orb = (UL_ONE << length); // maximum number of spherical orbits

   if(spBasZ->numj_orb > max_orb)   {
     printf("\n\nRank%d: Error in function pn_read_pn_two_part_matr_elem():",Rank);
      printf("\nToo many spherical proton  particle single-particle orbits"); 
      printf("\nnumj_orb = %d  -- must not exceed %d!",spBasZ->numj_orb, max_orb);
      printf("\n since an integer is only %d bytes\n", (int)sizeof(int));
      MPI_Abort(MPI_COMM_WORLD,Rank);
   }
   if(spBasN->numj_orb > max_orb)   {
     printf("\n\nRank%d: Error in function pn_read_pn_two_part_matr_elem():",Rank);
      printf("\nToo many spherical neutron single-particle orbits");
      printf("\nnumj_orb = %d  -- must not exceed %d!",spBasN->numj_orb, max_orb);
      printf("\n since an integer is only %d bytes\n", (int)sizeof(int));
      MPI_Abort(MPI_COMM_WORLD,Rank);
   }
  
   if( (file_ptr = fopen(file_name,"r")) == NULL) { // open input file
     printf("\n\nRank%d: Error in function pn_read_pn_two_part_matr_elem();",Rank);
      printf("\nWrong file = %s for input of eff.interaction.\n",file_name);
      MPI_Abort(MPI_COMM_WORLD,Rank);
   }

   if((fscanf(file_ptr,"%d%c", &num_elem, &ch) != 2) | (ch !=  '\n')) { // num of elem
     printf("\n\n Rank%d: Error in function pn_read_pn_two_part_matr_elem(): ",Rank);
      printf("\nFirst number in file %s", file_name);
      printf(" - number of matrix element - is not correct\n");
      MPI_Abort(MPI_COMM_WORLD,Rank);
   } 

   // memory to store two-particle matrix elements

   two_part = MALLOC(num_elem, TWO_PART, func,"j_pot.two_part[]");
   
   for(loop = 0, new_num_elem = 0, ptr = two_part; loop < num_elem; loop++) {

     // read (n1,l1,j1), (n2,l2,j2), (n3,l3,j3), (n4,l4,j4), MT, 2*J */

     if(read_int_number(file_ptr, 14, int_data) == FALSE) {
       printf("\n\nRank%d: Error in function pn_read_two_part_matr_elem();",Rank);
       printf("\nSomething is wrong with the (n,l,j) - values");
       printf(" for matrix element number = %d\n",loop);
       MPI_Abort(MPI_COMM_WORLD,Rank);
     }

     if(calc_CM == 0) { // read two-part VEFF only
       if(read_float_number(file_ptr, 1,ptr->value) == FALSE) {
	 printf("\n\nRank%d: Error in function pn_read_two_part_matr_elem();",Rank);
	 printf("\nWrong value for two-particle matrix element no %d",loop);
	 printf("\n read from file %s\n", file_name);
	 MPI_Abort(MPI_COMM_WORLD,Rank);
       }
     } // end read two-particle matrix elem CM = 0
     else  { // read both two-part VEFF and CM 
       if(read_float_number(file_ptr, 2, ptr->value) == FALSE) {
	 printf("\n\nRank%d: Error in function pn_read_two_part_matr_elem();",Rank);
	 printf("\nWrong value for the case of CENTER_OF_MASS %d",loop);
	 printf("\n read from file %s\n", file_name);
	 MPI_Abort(MPI_COMM_WORLD,Rank);
       }

       /********************  test run  ***********/
       ptr->value[1] *= -1.0;
       /*******************   rnf ttest run  *******/

     } // end read two-particle matrix elem CM = 1

        /*
	** converts (n,l.j) input data to orbit numbers from
	** jbas[] and store the result in structure orb_ident
	*/ 

     orb_ident = pn_matr_elem_identity(spBasZ, spBasN, int_data);

            /*
	    **  Only matrix elements relevant for the present
            **  set of single-particle orbits are saved
            */

     if(  (orb_ident.a == -1) || (orb_ident.b == -1) 
        ||(orb_ident.c == -1) || (orb_ident.d == -1)) continue;
 
     // calculate  config number for matrix element

     orb_ident.left  = (UL)((orb_ident.a << length) + orb_ident.b);
     orb_ident.right = (UL)((orb_ident.c << length) + orb_ident.d);       

              /* if necessary, interchange right and
	      ** left orbit identity such that
              **         orb_ident.left <= orb_ident.right
              */

     if( orb_ident.left > orb_ident.right) {
       temp            = orb_ident.left;
       orb_ident.left  = orb_ident.right;
       orb_ident.right = temp; 
     }

     // store config number together with matrix element

     ptr->config   = (orb_ident.left << (2*length)) + orb_ident.right;
     ptr->J        = int_data[13];
     ptr++;
     new_num_elem++;

   } // end loop through all input two-particle matrix elements

   fclose(file_ptr);

   if(new_num_elem < num_elem) { // reduce size of j_pot.two_part[]
     two_part = REALLOC(two_part, new_num_elem, TWO_PART, func,"realloc j_pot.two_part[]");
   }
   j_pot->num      = new_num_elem;  // save number of matrix element
   j_pot->two_part = two_part;      // and the corresponding pointer

         /* 
         ** All matrix elements are sorted into groups of same total J
         ** and for each group after increasing conf_ident.
         */

   num_J = (spBasZ->jbas[0].j + spBasN->jbas[0].j + 2) >> 1;
 
  pn_sort_matrix_elements(num_J, j_pot);

   // check that all necessary two-particle matrix elements are stored
 

   pn_check_matr_elem(spBasZ, spBasN, num_J, *j_pot);

} // End: function  pn_read_two_part_matr_elem()
                       
     /*
     ** The function 
     **         pn_matr_elem_identity()
     ** converts particle matrix elements identification from
     ** (n,l,j) form into an orbital number specified according
     ** to JBAS jbas[] for protons and neutrons. The function
     ** returns four orbital numbers in structure ORB_ID orb;  
     */

static ORB_ID pn_matr_elem_identity(SP_BAS *spBasZ, SP_BAS *spBasN,int *int_data)
{

   int        k;
   JBAS       *ptr;
   ORB_ID     orb;
 
   orb.a= -1; // local initialization
   orb.b= -1;
   orb.c= -1;
   orb.d= -1;

   for(k = 0, ptr = spBasZ->jbas; k < spBasZ->numj_orb; k++, ptr++)     {
     if(   (ptr->osc == int_data[0]) && (ptr->l == int_data[1]) 
	&& (ptr->j == int_data[2]))  orb.a = k;
     if(   (ptr->osc == int_data[6]) && (ptr->l == int_data[7]) 
        && (ptr->j == int_data[8]))  orb.c = k;
   } 
   for(k = 0, ptr = spBasN->jbas; k < spBasN->numj_orb; k++, ptr++)     {
     if(  (ptr->osc == int_data[3]) && (ptr->l == int_data[4]) 
        &&(ptr->j == int_data[5]))  orb.b = k;
     if(  (ptr->osc == int_data[9]) && (ptr->l == int_data[10]) 
        &&(ptr->j == int_data[11])) orb.d = k;
   }

   return orb;

} // End: function pn_matr_elem_identity()

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

   // sort matrix elements into groups after increasing total J-value

   qsort(j_pot->two_part, (UL) j_pot->num, sizeof(TWO_PART),
                      (int(*)(const void *, const void *)) pn_total_J_comp);

   // memory for j_pot->start[]

   j_pot->start = MALLOC(num_J + 1, int, func, "j_pot.start[]");

   for(jtot = 0, count = 0; jtot < num_J; jtot++) {
      j_pot->start[jtot] = count;
      while((count < j_pot->num) && (j_pot->two_part[count].J == (jtot << 1))) count++;

      // sort the matrix elements with J = tot_J after increasing config

      if((temp = count - j_pot->start[jtot]) > 0)   {
         qsort(j_pot->two_part + j_pot->start[jtot], (UL)(temp),
            sizeof(TWO_PART),(int(*)(const void *, const void *)) pn_config_comp);
      }
   }
   j_pot->start[jtot] = count;
   
} // End: function pn_sort_matrix_elements()

     /*
     ** The function                         
     **         pn_total_J_comp()                  
     ** is a utility function for the library function qsort() 
     ** in order to sort the two-body matrix elements into groups
     ** of increasing value of the total J-value.
     */

static int pn_total_J_comp(const TWO_PART *one, const TWO_PART *two)
{
  if(one->J > two->J)       return +1;
  else  if(one->J < two->J) return -1;
  else                        return  0;

} // End: function pn_total_J_comp()

     /*
     ** The function                         
     **        pn_config_comp()                  
     ** is a utility function for the library function qsort() in order
     ** to sort the two-body matrix elements for fixed total J into 
     ** groups of increasing conf_ident number.
     */

static int pn_config_comp(const TWO_PART *one, const TWO_PART *two)
{
  if(one->config > two->config)       return +1;
  else  if(one->config < two->config) return -1;
  else                      return  0;

} // End: function pn_config_comp()

     /*
     ** The function 
     **        pn_check_pn_matr_elem()
     ** checks current identical two-particle matrix elements to
     ** see if there are any missing compared to the available
     ** single-particle orbits.
     */

static void pn_check_matr_elem(SP_BAS *spBasZ, SP_BAS *spBasN, int num_J,
                                                               J_POT j_pot) 
{
   char     *func = {"pn_check_matr_elem: "};
   int      ja_orb, jb_orb, j_min, j_max, jtot, temp, control,
            parity,  *states_J_plus, *states_J_minus, *ptr;

   states_J_plus  = CALLOC(num_J, int, func, "pn_num_J_plus[]");    // local memory
   states_J_minus = CALLOC(num_J, int, func, "pn_num_J_minus[]");

   for(ja_orb = 0; ja_orb < spBasZ->numj_orb; ja_orb++) {
     for(jb_orb = 0; jb_orb < spBasN->numj_orb; jb_orb++) {
       j_min = (abs(spBasZ->jbas[ja_orb].j - spBasN->jbas[jb_orb].j)) >> 1;
       j_max = (spBasZ->jbas[ja_orb].j + spBasN->jbas[jb_orb].j) >> 1;
       parity = ((spBasZ->jbas[ja_orb].l + spBasN->jbas[jb_orb].l) % 2) ? -1 : +1;
       ptr    = (parity == +1) ? states_J_plus : states_J_minus;
       for(jtot = j_min; jtot <= j_max; jtot++)  ptr[jtot]++;

     } // end orb_b
   } // end orb_a 

   for(jtot = 0; jtot < num_J; jtot++)  {
     control = j_pot.start[jtot + 1] - j_pot.start[jtot];
     temp    = (states_J_plus[jtot] * (states_J_plus[jtot] + 1)) >> 1;
     temp   += (states_J_minus[jtot] * (states_J_minus[jtot] + 1)) >> 1;

     if(control != temp) {
       printf("\n\nRank%d: Error in function pn_check_pn_matr_elem(): ",Rank);
       printf("Not enough two-particle proton/neutron matrix input elements");
       printf(" for total_J = %d", jtot);
       printf("\nRead number = %d - should have been = %d\n", control, temp);
       MPI_Abort(MPI_COMM_WORLD,Rank);
     }
   }     

} // End: function pn_check_matr_elem()

     /*
     ** The function                                          
     **              pn_m_pot_diag(...)                    
     ** calculates all diagonal m-scheme two-particle matrix elements
     ** < p(k), n(l) | VEFF | p(k), l(n) > for proton-neutron case with 
     ** The diagonal matrix elements are  stored  in vector: 
     **          pn_diag[k * num_mN + l].diag_val
     */

static void pn_m_pot_diag(SP_BAS *spBasZ, SP_BAS *spBasN,
                           J_POT j_pot, int calc_CM, MATR_PN *op_int)
{
  int  k, l;
  double     pnDiag[2], *diagPtr,*CM_diagPtr;

  diagPtr    = op_int->pn_diag;
  CM_diagPtr = op_int->pn_CM_diag; 

  for(k = 0; k < spBasZ->numm_orb; k++) { // proton right
    for(l = 0; l < spBasN->numm_orb; l++) { // neutron right
      pn_m_scheme_diag_int(k,l, spBasZ->mbas, spBasN->mbas,
                                     j_pot,calc_CM, pnDiag);
      *diagPtr = pnDiag[0];
      diagPtr++;
      if(calc_CM == 1) {
	*CM_diagPtr = pnDiag[1];
	CM_diagPtr++;
      }
    } // end index l, neutron right
  } // end index k, proton right
   
} // End: function pn_m_pot_diag()

     /*
     ** The function                                  
     **               pn_m_scheme_diag_int()                   
     ** calculates and return diagonal proton-neutron effective interaction 
     ** in m-scheme
     **   <m_orb_a(p), m_orb_b(n) | Veff | m_orb_a(p), m_orb_b(n) > =                        
     **    SUM(jtot) ( 
     **     clebsch(j_a, m_a, j_b, m_b, jtot) * clebsch(j_a, m_a, j_b, m_b, jtot)         
     **    * <j_orb_a(p), j_orb_b(n), jtot| Veff | j_orb_a(p), j_orb_b(n), jtot >
     */

static void  pn_m_scheme_diag_int(int m_orb_a, int m_orb_b, MBAS *mbasZ,
                   MBAS *mbasN, J_POT j_pot, int calc_CM, double *pnDiag)
{
  int       j_a, m_a, j_b, m_b, j_min, j_max, jtot;
  UL        j_orb_a, j_orb_b, length, config, high, low, mid;
  double    cleb;

  pnDiag[0] = 0.0; // initialization
  pnDiag[1] = 0.0;

  j_orb_a = (UL) mbasZ[m_orb_a].orb;  // orbit numbers
  j_orb_b = (UL) mbasN[m_orb_b].orb;
  j_a = mbasZ[m_orb_a].j; // j and m values
  m_a = mbasZ[m_orb_a].m;
  j_b = mbasN[m_orb_b].j;
  m_b = mbasN[m_orb_b].m;
  j_min = (abs(j_a - j_b)) >> 1;
  j_max = (j_a + j_b) >> 1;
  length = 2*sizeof(UL);
  config =   (j_orb_a << (3 * length)) + (j_orb_b << (2 * length))
            + (j_orb_a << length) + j_orb_b;

  for(jtot = j_min; jtot <= j_max; jtot++)  {
    low  = j_pot.start[jtot];
    high = j_pot.start[jtot + 1];
    while(1)   {  //  binary search for configuration
      mid = (low + high)>>1;
      if(config < j_pot.two_part[mid].config)      high = mid -1;
      else if(config >j_pot.two_part[mid].config)  low = mid + 1;
      else                                         break;
    } // end of search for configuration identifier
    
    cleb =  clebsch_gordan(j_a, j_b, jtot << 1,m_a, m_b);
    pnDiag[0] += cleb *cleb * j_pot.two_part[mid].value[0];
    if(calc_CM == 1) {
      pnDiag[1] += cleb *cleb * j_pot.two_part[mid].value[1];
    }
  }            // end of angular mom loop

} // End: function pn_m_scheme_diag_int()

     /*
     ** The function                                          
     **              pn_m_pot_nondiag(,...)                    
     ** calculates all nondiagonal m-scheme proton-neutron two-particle
     ** matrix elements < p(i), n(j) | VEFF | p(k), l(n)>
     ** The nondiagonal matrix elements  are  stored in a list
     ** pn_nondiag_elem[] pointed to by  
     ** table[((k * num_mN+ l].start_nondiag
     ** The function returns number of calculated matrix elements.
     */

static int pn_m_pot_nondiag(int parZ, SP_BAS *spBasZ, SP_BAS *spBasN,
			    int calc_CM,J_POT pot, MATR_PN *op_int)
{
  int       num_mZ, num_mN, i, j, k, l, m_k, m_l,par_k, par_l, 
            del_m, del_par, num_m_par, num_parity, par_p,
            number, tot_num;
  MBAS      *mbasZ, *mbasN;
  PN_INT    **kl_table, *nondiag;

  num_mZ     = spBasZ->numm_orb;
  num_mN     = spBasN->numm_orb;
  num_parity = parZ ? 1 : 2;  // num parity change in matr.elem
  mbasZ      = spBasZ->mbas;  
  mbasN      = spBasN->mbas;

  nondiag   = op_int->nondiag_bas;

  tot_num = 0;

  num_m_par = 0;
  for(del_m = 0; del_m <=  spBasZ->mbas[0].m; del_m++){
 
   for(del_par = 0, par_p = +1; del_par < num_parity; 
                       del_par++, par_p = -1) {

       /* 
       ** for each del_m and del_p: calculate non-diagonal
       **  matrix elements with i(p) >= k(p) and j(n) >= l(n),
       */
    
     kl_table =  op_int->kl_list[num_m_par];
     num_m_par++;
     
      for(k = 0; k < num_mZ; k++) { // right proton part
	m_k   = mbasZ[k].m;
	par_k = mbasZ[k].par;
	for(l = 0; l < num_mN; l++, kl_table++) {  // right neutron part
	  m_l   = mbasN[l].m;
	  par_l = mbasN[l].par;

	  *kl_table = nondiag;
	  number = 0;   // counter for each (kl) - group


	  for(i = ((del_m == 0) && (del_par == 1)) ? 0 : k;

                                           i < num_mZ; i++) {
	    for( ;  !(  (del_m == ((m_k   - mbasZ[i].m) / 2))
                       &&(par_p == (par_k * mbasZ[i].par)))      
                       && ( i < num_mZ); i++);
	    if(i == num_mZ) break;
	    for(j = (i == k) ? l + 1: 0; j < num_mN; j++)   {
	      for( ;   !(  (del_m == ((mbasN[j].m - m_l) / 2)) 
                          &&(par_p == (mbasN[j].par * par_l)))
                          && ( j < num_mN); j++);
	      if(j == num_mN) break; 

              pn_m_scheme_nondiag_veff(i,j,k,l, mbasZ, mbasN, calc_CM,
                                                     pot,nondiag->val);

	      if(  (fabs(nondiag->val[0])
                 + fabs(nondiag->val[1])) > MATRIX_LIMIT) {
		nondiag->one[0] =  ULL_ONE<<i;
		nondiag->one[1] =  ULL_ONE<<j;
		nondiag->two[0] 
                   = (i == k) ? 0 
                              : (ULL_ONE<<MAX(i,k))-(ULL_ONE<<(MIN(i,k) + 1));
		nondiag->two[1] 
                   = (j == l) ? 0 
                              : (ULL_ONE<<MAX(j,l)) - (ULL_ONE<<(MIN(j,l)+1));
		number++;
		tot_num++;
		nondiag++;
	      }
	    }  // end loop through j-neutrons
	  } // end loop through i-protons
	  if(number)   { 
	    nondiag->one[0]     = 0;// values to terminate a block of matr elem
	    (nondiag++)->one[1] = 0;
	    tot_num++;
	  }     
	  else   {  // no matr elem in present block
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
     **               pn_pn_m_scheme_nondiag_veff()                   
     ** calculates and return proton-neutron effective interaction in 
     ** m-scheme ( i <= j and k <= l) 
     **   < m_orb_a(p) m_orb_b(n) | Veff | m_orb_c(p) m_orb_d(n) > =                        
     **          SUM(J) (  clebsch(j_a m_a j_b m_b J)         
     **                 * clebsch(j_c m_c j_d m_d J)         
     **                 * phase( ji + jk -jj - jl)           
     **                 * < j_orb_a(p) j_orb_b(n) J| Veff | j_orb_c(p) j_orb_d(n) J > 
     */

static void pn_m_scheme_nondiag_veff(int m_orb_a, int m_orb_b, int m_orb_c,
                                int m_orb_d, MBAS *mbasZ, MBAS *mbasN, 
                          int calc_CM,J_POT j_pot, double *matrElem)
{
   int     j_a, m_a, j_b, m_b, j_c, m_c, j_d, m_d,
           j_min, j_max, jtot, length;
   UL      j_orb_a, j_orb_b, j_orb_c, j_orb_d, conf_bra,
           temp, conf_ket, config, high, low, mid;
   double  factor;

   matrElem[0] = 0.0; // initialization
   matrElem[1] = 0.0;

  j_orb_a = (UL)mbasZ[m_orb_a].orb; // orbit numbers
   j_orb_b = (UL)mbasN[m_orb_b].orb;
   j_orb_c = (UL)mbasZ[m_orb_c].orb;
   j_orb_d = (UL)mbasN[m_orb_d].orb;  

   j_a = mbasZ[m_orb_a].j; // j and m values
   m_a = mbasZ[m_orb_a].m;
   j_b = mbasN[m_orb_b].j;
   m_b = mbasN[m_orb_b].m;
   j_c = mbasZ[m_orb_c].j;
   m_c = mbasZ[m_orb_c].m;
   j_d = mbasN[m_orb_d].j;
   m_d = mbasN[m_orb_d].m;

   // max and min values of total J*2

   j_min = MAX(abs(j_a - j_b), abs(j_c - j_d));
   j_max = MIN(j_a + j_b,  j_c + j_d );

   if(j_min > j_max )  return;

   length   = sizeof(UL) << 1;
   conf_bra = (j_orb_a << length) + j_orb_b;
   conf_ket = (j_orb_c << length) + j_orb_d;

   if (conf_bra > conf_ket) { 

     // possible interchange left and right of matr elem

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

      factor =  clebsch_gordan(j_a,j_b,jtot,m_a,m_b)
	      * clebsch_gordan(j_c,j_d,jtot,m_c,m_d);
      matrElem[0] += factor * j_pot.two_part[mid].value[0];
      if(calc_CM == 1) {
	matrElem[1] += factor * j_pot.two_part[mid].value[1];
      }
   } // end loop over all jtot

   // "Whitehead time-reversal" phase is included 

   matrElem[0] *=  mbasZ[m_orb_a].phase * mbasN[m_orb_b].phase
                 * mbasZ[m_orb_c].phase * mbasN[m_orb_d].phase;


   if(calc_CM == 1) {
     matrElem[1] *=  mbasZ[m_orb_a].phase * mbasN[m_orb_b].phase
                   * mbasZ[m_orb_c].phase * mbasN[m_orb_d].phase;
   }
} // End: function pn_m_scheme_nondiag_veff()

     /*
     ** The function
     **       addPN_1to2SinglePartTermVeff()
     ** adds for proton-neutron particles  contributions from 
     ** single-particle terms to the effective two-particle 
     ** matrix elements  <k.l |int| k.l>
     */
static void addPN_1to2SinglePartTermVeff(int type, int partNum,int num_mZ, 
                      int num_mN, MBAS *mbas, int calc_CM,MATR_PN *op_int)
{
   int     k, l;
   double  spEnergy,  spEnergyCM, factor1to2,
           *pnDiag,*pn_CM_diag;

   pnDiag = op_int->pn_diag;
   if(calc_CM == 1) {
     pn_CM_diag = op_int->pn_CM_diag;
   }
   else {
     pn_CM_diag = (double *)NULL_PTR;
   }

   factor1to2 = 1.0/((double)(partNum));

   for(k = 0; k < num_mZ; k++) {
     if(type == PROTON)  {
       spEnergy = mbas[k].e * factor1to2;
       if(calc_CM == 1) {
	  spEnergyCM = (2*mbas[k].osc + mbas[k].l + 1.5)
	              * factor1to2;
       }
       else {
	 spEnergyCM = D_ZERO;
       }

     }
     for(l = 0; l < num_mN; l++) {
       if(type == NEUTRON) {
	 spEnergy = mbas[l].e * factor1to2;


	 if(calc_CM == 1) {
	   spEnergyCM = (2*mbas[l].osc + mbas[l].l + 1.5)
	               * factor1to2;
	 }
	 else {
	   spEnergyCM = D_ZERO;
	 }
       }

       *(pnDiag++) += spEnergy;
       if(calc_CM == 1) {
	 *(pn_CM_diag++) += spEnergyCM;
       }
     } // end l-loop
   } //end k-loop

} // End: function  addPN_1to2SinglePartTermVeff()
