/************************  File lanc-pn-store.c *************/

#include "shell.h"

     /*
     ** The module entrance function 
     **       pn_store_SD_matr_elem()
     ** calculates and stores on file all diagonal and as many
     ** as possible of the nondiagonal SD matrix elements
     **   <SD(Z)g',p';SD(N)g',n'|OP|SD(Z)g,p;SD(N)g,n>
     ** The parameter op_type determines the two-particle 
     ** interaction OP.
     */


            /**** local data definitions ****/

typedef   struct {
  UL
                         sd;
      int
                     sd_num,
                      phase,
                       mask;
   } RESULT;

     typedef   struct {
       int
                      init_gr,
	             final_gr,
                  stored_diag,
	       stored_nondiag;
     } HEADER;

	    /**** local function declarations ****/

static void pn_occupied_orbits(int part_type, GR_BASIS *gr_basis, GROUP *gr);
    /*
    ** calculates and stores orbit numbers of all occupied orbits  in all 
    ** |SD(Z)> (part_type = 0) and |SD(N)> (part_type = 1) in present group gr[] 
    */

 static int pn_diag_SD_elem(int op_type, GR_BASIS *gr_basis, GROUP *gr);
    /*
    ** calculates the diagonal SD matrix elements 
    **   <SD(Z)g,p;SD(N)g,n|OP|SD(Z)g,p;SD(N)g,n>
    ** with OP = OP(pp) + OP(nn) and OP(pn).
    ** The result is stored in memory pointed to by gr_basis->mem_free.
    ** The function returns the number of matrix elements stored.
    */

static void pn_diag_pp_matr_elem(int op_type, GR_BASIS *gr_basis, GROUP *gr);
    /*
    ** calculates diagonal matrix elements 
    **           <SD(N),SD(Z)|OP(pp)| SD(Z),SD(N)> 
    ** where the two-particle matrix elements of OP(pp) are stored in table[0]
    ** using the formula
    **        op(k,l) = ((num - k) * k)/2 -1 +l
    ** where num = 2 * gr.num_orb[type] - 3. op_type = 0 refers to VEFF and 
    ** op_type = 1 to ANG. The result is stores in diag_SD[]. 
    */

static void  pn_diag_nn_matr_elem(int op_type, GR_BASIS *gr_basis, GROUP *gr);
    /*
    ** calculates diagonal matrix elements 
    **           <SD(N),SD(Z)|OP(nn)| SD(Z),SD(N)> 
    ** where the two-particle matrix elements of OP(nn)
    ** are stored in table[1] using the formula
    **        op(k,l) = ((num - k) * k)/2 -1 +l
    ** where num = 2 * gr.num_orb[type] - 3.
    ** op_type = 0 refers to VEFF and op_type = 1 to ANG.
    ** The result is stores in diag_SD[]. 
    */

static void  pn_diag_pn_matr_elem(int op_type, GR_BASIS *gr_basis, GROUP *gr);
    /*
    ** calculates diagonal SD matrix elements 
    **    <SD(N),SD(Z)|OP(pn)| SD(Z),SD(N)> 
    ** The operator matrix elements are stored in table[] by the formula
    **        op(k,l) = num_orbN * k + l
    ** and the result is stores in diag_SD[]. 
    */

static int pn_nondiag_SD_elem(int op_type, GR_BASIS *gr_basis, GROUP *group);
    /*
    ** calculates  all nondiagonal SD matrix elements
    **           <SD(Z)g,p';SD(N)g,n'|OP|SD(Z)g,p;SD(N)g,n>
    ** for a single group (g) if enough memory pointed to by gr_basis->mem_free
    ** is available. The two-particle interactions is defined through
    ** the parameter op_type and includes  OP(pp), OP(nn) and OP(pn).
    ** The function returns -1 if not enough memory is available to store
    ** the matrix elements, Otherwise it returns the number of nondiagonal
    **  matrix elements stored.
    */

static int pn_nondiag_id_SD_elem(int op_type, int part_type, GR_BASIS *gr_basis,
                                                            GROUP *gr, int *pos);
    /*
    ** calculates all nondiagonal identical particle matrix elements 
    **      <SD(N)g.n',SD(Z)g.p'|OP| SD(Z)g,p,SD(N)g,n> 
    ** within one proton/neutron group only. Pointers to all two-particle
    ** matrix elements for particle type = PROTON(NEUTRON) are stored in 
    ** structure gr_basis->op_int.id_table[] by the formula
    **        op(k,l) = ((num - k) * k)/2 -1 +l
    ** where num = 2 * gr.num_orb[type] - 3. The nondiag matrix elements
    ** are  saved in memory pointed to by gr_basis->mem_free.
    ** The function returns number of matrix elements stored.
    */

static int sort_nondiag(const STORE *one, const STORE *two);
    /*
    ** is a utility function for the library function qsort() in order to
    ** sort nondiagSD matrix elements of type STORE store[]
    ** after increasing store[].final.
    */

static int compress_nondiag_SD_elem(int num, STORE *store);
    /*
    ** takes num calculated nondiag SD elem which are sorted
    ** after increasing store[].final and add together all
    ** elements  with same store[].final.
    ** Note: The function must have num > 1;
    ** The final number of elements is returned.
    */

static int pn_nondiag_pn_SD_elem_A(int op_type, GR_BASIS *gr_basis, GROUP *gr,
                                  int *pp_pos, STORE *matr_pp, STORE *matr_nn);
    /*
    ** calculates all nondiagonal pn matrix elements 
    **  <SD(Z)g,p';SD(N)g,n'|OP(pp) + OP(nn) + OP(pn)|SD(Z)g,p,SD(N)g,n>
    ** within a single group g. The two-particle matrix elements are
    ** specified through the parameter op_type.
    ** The calculation is performed in following steps. For each 
    ** init |SD(Z),SD(N)>:
    **   1. Convert previous calculated nondiagonal pp-matrix elements
    **      into pn format and store them in memory pointed to by 
    **      gr_basis->mem_free. Update gr_basis->mem_free.
    **   2. Convert previous calculated nondiagonal nn-matrix elements
    **      into pn format and store them in memory pointed to by 
    **      gr_basis->mem_free. Update gr_basis->mem_free.
    **   3. Calculate nondiagonal pn-matrix elements and store them
    **      in memory pointed to by gr_basis->mem_free.
    **   4. If the remaining free memory is too small for the matrix
    **      elements the function terminates and return -1.
    **   5. Sort after final |SD(Z),SD(N)> and compress all nondiag
    **      matrix elements.
    **   6. The function returns number of matrix elements stored.
    */   

static void pn_write_SD_matr_elem(char *file_name, int new_file,HEADER header,char * mem_ptr);
    /*
    ** opens the file "file_name" and write/append diagonal
    ** and nondiagonal matrix elements 
    **       <SD(Z)g',p: SD(N)g',n|OP|SD(Z)g,p: SD(N)g,n>
    ** for (g,p), (g,n) --> (g',p'),(g',n')
    ** if(header.store_nondiag == -1) no nondiagonal matrix
    ** elements are calculated due to lack of memory.
    ** If(new_file == YES) a new file is created.
    */

static int pn_nondiag_pn_SD_elem_B(int op_type, GR_BASIS *gr_basis, GROUP *init_gr,
                                                                    GROUP *final_gr);
    /*
    ** calculates all nondiagonal pn-matrix elements 
    **  <SD(Z)g',p';SD(N)g',n'|OP(pp) + OP(nn) + OP(pn)|SD(Z)g,p,SD(N)g,n>
    ** between two groups g -> g'. The two-particle matrix elements are
    ** specified through the parameter op_type.
    ** The calculation is performed in following steps. For each 
    ** init |SD(Z),SD(N)>:
    **   1. Calculate nondiagonal pn-matrix elements and store them
    **      in memory pointed to by gr_basis->mem_free.
    **   2. If the free memory is too small for the matrix
    **      elements the function terminates and return -1.
    **   3. Sort after final |SD(Z),SD(N)> and compress all nondiag
    **      matrix elements.
    **   4. The function returns number of matrix elements stored.
    */   

                /**** End: function declarations ****/

               /**** The function definitions  ****/ 
 
     /*
     ** The entrance function
     **       pn_store_SD_matr_elem()
     ** calculates and stores on file all diagonal and as many
     ** as possible of the nondiagonal SD matrix elements
     **   <SD(Z)g',p';SD(N)g',n'|OP|SD(Z)g,p;SD(N)g,n>
     ** The parameter op_type determines the two-particle 
     ** matrix elements in  operator OP.
     */

void pn_store_SD_matr_elem(int op_type, GR_BASIS *gr_basis, GROUP *init_gr)
{
   char
            *func = {"pn_store_SD_matr_elem(): "},
            *file_name; 

   int
            new_file, init_gr_num, final_gr_num,
            stored_diag, stored_nondiag, dim, max_delM_Z;
   HEADER
            header;
   GROUP
            *final_gr;

            /* maximum change in M_Z between initial and final groups */

   max_delM_Z = 2 * gr_basis->mbas[0]->m;

             /* file to store SD matrix elements */

   file_name = op_type ? gr_basis->store_ang : gr_basis->store_veff;

        /* memory to store occupied orbits in all |SD(Z)> and SD(N)> in a group */

   gr_basis->occ[0] = MALLOC(gr_basis->maxSD[0] * gr_basis->part[0], int,func,"p_occ[]");
   gr_basis->occ[1] = MALLOC(gr_basis->maxSD[1] * gr_basis->part[1], int,func,"n_occ[]");

           /* memory to store both diagSD[] and nondiagSD[] matr elem  - unit char */

   dim                = ((gr_basis->maxSD[0] * gr_basis->maxSD[1]) * sizeof(double) 
                       + gr_basis->mem_nondiag * sizeof(STORE)) / sizeof(char);
   gr_basis->mem_pos  = MALLOC(dim, char, func,"mem_ptr[]");

       /* loop through all init group basis states |SD(Z),SD(N)> */ 

   new_file = YES;           /* create new storage files for SD matr elem */
   

/********  test output ****************
     printf("\n\n\n\nCalculation and storing of %s  SD matrix elements\n\n\n:",
                                                   (op_type ? "ANG" : "VEFF" ));
**************  end *******/


   for(init_gr_num = 0; init_gr_num < gr_basis->num_gr; init_gr_num++, init_gr++) {
     
     if((init_gr->numSD[0] * init_gr->numSD[1]) == 0)  continue;  /* no basis states */

/********  test output ****************
     printf("\n\nCalculation of SD matrix elements for init_gr = final_gr = %d", init_gr_num);
**************  end *******/

              /* find all occupied p-orbits and n_orbits  in group[k] */
 
     pn_occupied_orbits(PROTON, gr_basis, init_gr);                    /* p-orbits */
     pn_occupied_orbits(NEUTRON, gr_basis, init_gr);                   /* n-orbits */ 

          /* 
	  ** Diagonal matrix elements:
          ** The function returns the number of
	  ** bytes reserved in gr_basis->mem_pos.
	  */

     stored_diag = pn_diag_SD_elem(op_type, gr_basis, init_gr);

          /*
	  ** calculate nondiag matrix elements of <SD()|OP|SD()> for
          **                     final_gr = init_gr 
          */

     if(gr_basis->mem_nondiag > 0) {

                /* update free memory */

       gr_basis->mem_free = gr_basis->mem_pos + (stored_diag * sizeof(double))/sizeof(char);
       stored_nondiag     = pn_nondiag_SD_elem(op_type, gr_basis, init_gr);
     }
     else {
       stored_nondiag = -1; /* No storage of nondiag matrix elements */
     }
     header.init_gr        = init_gr_num;
     header.final_gr       = init_gr_num;
     header.stored_diag    = stored_diag;
     header.stored_nondiag = stored_nondiag;

     pn_write_SD_matr_elem(file_name, new_file, header, gr_basis->mem_pos);
     new_file = NO;

          /*
	  ** calculate and store nondiag matrix elements of <SD()|OP|SD()> for
          **                    final_gr > init_gr 
          */

     for(final_gr_num = init_gr_num + 1, final_gr = init_gr + 1;
                 final_gr_num  <  gr_basis->num_gr; final_gr_num++, final_gr++) {
     
       if(   ((final_gr->numSD[0] * final_gr->numSD[1]) == 0)    /* no basis states */
          || ((init_gr->m[0] - final_gr->m[0]) > max_delM_Z))   /* too large del_M_Z */
	 continue;

/********  test output ****************
     printf("\n\nCalculation of SD matrix elements for init_gr = %d and  final_gr = %d", 
                                           init_gr_num, final_gr_num);
**************  end *******/



       if(gr_basis->mem_nondiag > 0) {
	 stored_nondiag = pn_nondiag_pn_SD_elem_B(op_type, gr_basis, init_gr, final_gr);
       }
       else {
	 stored_nondiag = -1; /* No storage of nondiag matrix elements */
       }
       header.init_gr        = init_gr_num;
       header.final_gr       = final_gr_num;
       header.stored_diag    = 0;
       header.stored_nondiag = stored_nondiag;

       pn_write_SD_matr_elem(file_name, new_file, header, gr_basis->mem_pos);

     }/* end loop through final groups */
   } /* end loop through init groups */

   free(gr_basis->mem_pos);                 /* remove temporary memory */
   free(gr_basis->occ[1]);
   free(gr_basis->occ[0]);

   gr_basis->mem_pos  = NULL;
   gr_basis->mem_free = NULL;
   gr_basis->occ[0]   = NULL;
   gr_basis->occ[1]   = NULL;

/***************  test output  *************
    printf("\n\nTermination of SD matrix element storing\n\n\n\n");  
**************   end test output  ************/

} /* End: function pn_store_SD_matr_elem() */

    /* 
    ** The function
    **       pn_occupied_orbits()
    ** calculates and stores orbit numbers of all occupied orbits  in all 
    ** |SD(Z)> (part_type = 0) and |SD(N)> (part_type = 1) in present group gr[] 
    */

static void pn_occupied_orbits(int part_type, GR_BASIS *gr_basis, GROUP *gr)
{
   int
          num_part, *occ, k, par, orb;
   UL
          *sd, pos, sd_value;
  
   num_part = gr_basis->part[part_type];         /* initialization */
   occ      = gr_basis->occ[part_type];
   sd       = gr->SD[part_type];
   k        = gr->numSD[part_type];
   do {
      sd_value = *(sd++);
      for(par = 0, orb = 0, pos = UL_ONE; par < num_part; par++, orb++, pos <<= 1)  {
         for(;!(sd_value & pos); orb++, pos <<= 1);                /* particle found */ 
         *(occ++) = orb;                                         /* save orbit */ 
      }  
   } while(--k); /* end k-loop through all |SD> in sd[] */        

} /* End: function pn_occupied_orbits() */

    /*
    ** The function 
    **         pn_diag_SD_elem()
    ** calculates the diagonal SD matrix elements 
    **   <SD(Z)g,p;SD(N)g,n|OP|SD(Z)g,p;SD(N)g,n>
    ** with OP = OP(pp) + OP(nn) and OP(pn).
    ** The result is stored in memory pointed to by gr_basis->mem_free.
    ** The function returns the number of matrix elements stored.
    */

static int pn_diag_SD_elem(int op_type, GR_BASIS *gr_basis, GROUP *gr)
{
   int
             dim, k;
   double
             *ptr_diagSD;

   dim = gr->numSD[PROTON] * gr->numSD[NEUTRON];       /* dimension of |SD(Z),SD(N)> */

   ptr_diagSD = (double *)gr_basis->mem_pos;          /* temporary storage position */
   k = dim;
   while(k--) *(ptr_diagSD++) = D_ZERO;            /* initialization of diagSD[] */

   if(gr_basis->part[PROTON] >= 2) {  /* <|OP(pp)|> contr to diag matr elem */
      pn_diag_pp_matr_elem(op_type, gr_basis, gr);
   }
   if(gr_basis->part[NEUTRON] >= 2) {      /* <|OP(nn)|> contr to diag matr elem */
      pn_diag_nn_matr_elem(op_type, gr_basis, gr);
   }

   pn_diag_pn_matr_elem(op_type, gr_basis, gr); /* <|OP(pn)|> contr to diag matr elem */

   /***********  test output  **************
   ptr_diagSD = (double *)gr_basis->mem_pos;
   printf("\ndiag matr elem:\n");
   for(k = 0; k < dim; k++) {
     printf(" %8.4f ",*(ptr_diagSD++));
     if((((k + 1)/8) * 8) == ( k + 1)) printf("\n");
   }
*******************  end  ******/

   return dim ;    /* number of stored matrix elements */  

} /* End: function pn_diag_SD_elem() */

     /*
     ** The function 
     **         pn_diag_pp_matr_elem()
     ** calculates diagonal matrix elements 
     **           <SD(N),SD(Z)|OP(pp)| SD(Z),SD(N)> 
     ** where the two-particle matrix elements of OP(pp) are stored in
     ** id_diag[0] using the formula
     **        op(k,l) = ((num - k) * k)/2 -1 +l
     ** where num = 2 * gr.num_orb[type] - 3. op_type = 0 refers to VEFF and 
     ** op_type = 1 to ANG. The result is stores in diag_SD[]. 
     */

static void pn_diag_pp_matr_elem(int op_type, GR_BASIS *gr_basis, GROUP *gr)
{
   int
                  k, l, *k_ptr, *l_ptr, numSD_Z, numSD_N, part_1, num;  
   double
                   *table, *table_ptr, value, *diag_ptr; 

   switch(op_type) {
      case VEFF: table = gr_basis->op_int.id_diag[PROTON];  /* <|H(pp)|> */
                 break;
      case ANG:  table = gr_basis->op_ang.id_diag[PROTON];  /* <|J**2(pp)|>*/ 
                 break;
      default:   printf("\n\nError in function pn_diag_pp_matr_elem():");
                 printf("\nIncorrect two-particle matrix element table");
                 printf(" type = %d -- should be 1 - 2\n", op_type);
                 exit(1);
   } 
   part_1   = gr_basis->part[PROTON] - 1;                     /* initialization */
   k_ptr    = gr_basis->occ[PROTON];
   num      = (gr_basis->m_orb[PROTON] << 1) - 3;
   diag_ptr = (double *)gr_basis->mem_pos;
   numSD_Z  = gr->numSD[PROTON]; 
   numSD_N  = gr->numSD[NEUTRON];
   do {
      value = D_ZERO;                             /* calculate <SD(Z)|OP(pp)|SD(Z)> */
      k = part_1;
      do {                                                       /* k-particle loop */
         table_ptr = table + (((num - (*k_ptr)) * (*k_ptr)) >> 1) -1;
         l     = k;
         l_ptr = (++k_ptr);
         do {                                                    /* l-particle loop */ 
             value += *(table_ptr + (*(l_ptr++))); 
         } while(--l);                                           /* end l particle loop */ 
      } while(--k);                                              /* end k particle loop */
      k = numSD_N;      /* store <SD(Z)|OP(pp)|SD(Z)in all config. with same |SD(Z)> */
      while(k--) *(diag_ptr++) += value;      /* add contribution to diagSD[ ] */ 
      k_ptr++;                                               /* move to next |SD(Z)> */
   } while(--numSD_Z);     /* end  loop through all |SD(Z)> */
} /* End: function  pn_diag_pp_matr_elem() */

     /*
     ** The function 
     **         pn_diag_nn_matr_elem()
     ** calculates diagonal matrix elements 
     **           <SD(N),SD(Z)|OP(nn)| SD(Z),SD(N)> 
     ** where the two-particle matrix elements of OP(nn)
     ** are stored in id_diag[1] using the formula
     **        op(k,l) = ((num - k) * k)/2 -1 +l
     ** where num = 2 * gr.num_orb[type] - 3.
     ** op_type = 0 refers to VEFF and op_type = 1 to ANG.
     ** The result is stores in diag_SD[]. 
     */

static void  pn_diag_nn_matr_elem(int op_type, GR_BASIS *gr_basis, GROUP *gr)
{
   int
                 loop_n, k, l, *k_ptr, *l_ptr, numSD_Z, numSD_N, part_1, num;  
   double
                 *table, *table_ptr, value, *diagSD, *diag_ptr; 

   switch(op_type) {
      case VEFF: table = gr_basis->op_int.id_diag[NEUTRON];  /* <|H(nn)|> */
                 break;
      case ANG:  table = gr_basis->op_ang.id_diag[NEUTRON];  /* <|J**2(nn)|> */ 
                 break;
      default:   printf("\n\nError in function pn_diag_nn_matr_elem():");
                 printf("\nIncorrect two-particle matrix element table");
                 printf(" type = %d -- should be 1 - 3\n", op_type);
                 exit(1);
   } 
   part_1  = gr_basis->part[NEUTRON] - 1;             /* initialization */
   k_ptr   = gr_basis->occ[NEUTRON];
   num     = (gr_basis->m_orb[NEUTRON] << 1) - 3;
   diagSD  = (double *)gr_basis->mem_pos;
   numSD_Z = gr->numSD[PROTON];
   numSD_N = gr->numSD[NEUTRON];
   loop_n = numSD_N;
    do {
      value = D_ZERO;                 /* calculate <SD(N)|OP(nn)|SD(N)> */
      k = part_1;
      do {                                                       /* k-particle loop */
         table_ptr = table + (((num - (*k_ptr)) * (*k_ptr)) >> 1) -1;
         l     = k;
         l_ptr = (++k_ptr);
         do {                                                    /* l-particle loop */ 
             value += *(table_ptr + (*(l_ptr++))); 
         } while(--l);                                           /* end l particle loop */ 
      } while(--k);                                             /* end k particle loop */

      diag_ptr = diagSD;  /* store <SD(Z)|OP(nn)|SD(Z) in all config. with same |SD(N)> */
      k = numSD_Z;
      do {
         *diag_ptr += value;      /* add contribution to diagSD[ ] */ 
         diag_ptr += numSD_N;
      } while(--k);
      k_ptr++;
      diagSD++;
    } while(--loop_n);     /* end  loop through all |SD(N)> */

} /* End: function  pn_diag_nn_matr_elem() */

     /*
     ** The function 
     **         pn_diag_pn_matr_elem()
     ** calculates diagonal SD matrix elements 
     **    <SD(N),SD(Z)|OP(pn)| SD(Z),SD(N)> 
     ** The operator matrix elements are stored in pn_diag[] by the formula
     **        op(k,l) = num_orbN * k + l
     ** and the result is stores in diag_SD[]. 
     */

static void  pn_diag_pn_matr_elem(int op_type, GR_BASIS *gr_basis, GROUP *gr)
{
   int
               p_sd, n_sd, k, l, *k_ptr, *l_ptr, numSD_Z, numSD_N,
               num_partZ, num_partN, num_orbN, *occZ, *occN;  
   double 
                *table, *pn_table, value, *diag_ptr; 

   switch(op_type) {
      case VEFF: table = gr_basis->op_int.pn_diag;   /* < |H(pn)| > contribution */
                 break;
      case ANG:  table = gr_basis->op_ang.pn_diag;   /* < |J**2(pn)| > contribution */ 
                 break;
      default:   printf("\n\nError in function pn_diag_pn_matr_elem():");
                 printf("\nIncorrect two-particle matrix element table:");
                 printf(" type = %d -- should be 1 - 2\n", op_type);
                 exit(1);
   } 
   num_partZ = gr_basis->part[PROTON];                /* initalization */
   num_partN = gr_basis->part[NEUTRON];
   num_orbN  = gr_basis->m_orb[NEUTRON];
   occZ      = gr_basis->occ[PROTON];
   occN      = gr_basis->occ[NEUTRON];
   diag_ptr  = (double *)gr_basis->mem_pos;
   numSD_Z   = gr->numSD[PROTON];
   numSD_N   = gr->numSD[NEUTRON];
   for(p_sd = 0; p_sd < numSD_Z; p_sd++) {        /* proton loop */
      for(n_sd = 0; n_sd < numSD_N; n_sd++) {     /* neutron loop */
         value = D_ZERO;
         k_ptr = occZ + p_sd * num_partZ;
         k     = num_partZ;
         do {                                                 /* proton loop */
	    pn_table = table + (*(k_ptr++)) * num_orbN;
            l_ptr    = occN + n_sd * num_partN;
            l = num_partN;
            while(l--) value += *(pn_table + (*(l_ptr++))); /* neutron loop */
         } while(--k);
         *(diag_ptr++) += value;                      /* cont. from (kl) pair */
      } /* end run through all |SD(N)> */
   } /* end run through all |SD(Z)> */
} /* End: function  pn_diag_pn_matr_elem() */

    /* 
    ** The function 
    **         pn_nondiag_SD_elem()
    ** calculates  all nondiagonal SD matrix elements
    **           <SD(Z)g,p';SD(N)g,n'|OP|SD(Z)g,p;SD(N)g,n>
    ** for a single group (g) if enough memory pointed to by gr_basis->mem_free
    ** is available. The two-particle interactions is defined through
    ** the parameter op_type and includes  OP(pp), OP(nn) and OP(pn).
    ** The function returns -1 if not enough memory is available to store
    ** the matrix elements, Otherwise it returns the number of nondiagonal
    **  matrix elements stored.
    */
   
static int pn_nondiag_SD_elem(int op_type, GR_BASIS *gr_basis, GROUP *group)
{
  char
                *func = {"pn_nondiag_SD_elem(): "};
  int
                numsdZ, numsdN, k, num_pp, *pp_pos, num_nn, *nn_pos,
                stored_matr_elem, dim;
  STORE       
                *matr_pp, *matr_nn, *ptr_1, *ptr_2;

   numsdZ   = group->numSD[0];                              /* initialization */
   numsdN   = group->numSD[1];

              /* memory limit for nondiag id SD matrix elements */

   dim = MAX(numsdZ, numsdN);
   dim = (dim * (dim + 1)) / 2;
   if(dim  > gr_basis->mem_nondiag)  return -1;    /* not enough memory */

           /* 
	   ** calculate nondiagSD proton-proton matrix  
	   ** elements between states in current group
           */

   if((gr_basis ->part[PROTON] < 2) || (numsdZ < 2)) {
     pp_pos  = MALLOC(numsdZ, int, func,"pp_pos[]");  /* no pp-matrix elements */
     num_pp  = 0;
     matr_pp = MALLOC(numsdZ, STORE, func,"matr_pp[]");

     for(k = 0; k < numsdZ; k++) {                   /* necessary initialization */ 
       pp_pos[k]        =  0; 
       matr_pp[k].final = -1;
     }
   }
   else {                                            /* calculate pp-matrix elements */
     pp_pos  = MALLOC(numsdZ, int, func,"pp_pos[]");
     num_pp  = pn_nondiag_id_SD_elem(op_type, PROTON, gr_basis, group, pp_pos);

     if(num_pp == -1) return -1;    /* no more memory to store matr elem */

        /* store temporarily the proton-proton matrix elements */

     matr_pp = MALLOC(num_pp, STORE, func,"matr_pp[]");
     ptr_1   = matr_pp;
     ptr_2   = (STORE *)gr_basis->mem_free;
     k       = num_pp;
     do {
       ptr_1->value     = ptr_2->value;
       (ptr_1++)->final = (ptr_2++)->final;
     } while(--k);
   } /* end proton case */
	
           /* 
	   ** calculate nondiagSD neutron-neutron matrix  
	   ** elements between states in current group
           */

   if((gr_basis ->part[NEUTRON] < 2 ) || (numsdN < 2)) {  /* no pp-matrix elements */
     nn_pos  = MALLOC(numsdN, int, func,"nn_pos[]");
     num_nn  = 0; 
     matr_nn = MALLOC(numsdN, STORE, func,"matr_nn[]");

     for(k = 0; k < numsdN; k++) {                       /* necessary initialization */ 
       nn_pos[k]        =  0; 
       matr_nn[k].final = -1;
     }
   }
   else {                                               /* calculate pp-matrix elements */
     nn_pos  = MALLOC(numsdN, int, func,"nn_pos[]");
     num_nn  = pn_nondiag_id_SD_elem(op_type, NEUTRON, gr_basis, group, nn_pos);

     if(num_nn == -1) return -1;    /* no more memory to store matr elem */

           /* store temporarily the neutron-neutron matrix elements */

     matr_nn = MALLOC(num_nn, STORE, func,"matr_nn[]");
     ptr_1   = matr_nn;
     ptr_2   = (STORE *)gr_basis->mem_free;
     k       = num_nn;
     do {
       ptr_1->value     = ptr_2->value;
       (ptr_1++)->final = (ptr_2++)->final;
     } while(--k);
   } /* end neutron case */
	
	 /* memory check for pn matr elem from pp- and nn-nondiag matr */
   
   stored_matr_elem = -1;                                          /* initialization */
   if((num_pp * numsdN + num_nn * numsdZ) < gr_basis->mem_nondiag) {
     stored_matr_elem = pn_nondiag_pn_SD_elem_A(op_type, gr_basis, group, 
	  		                        pp_pos, matr_pp, matr_nn);
   }
   free(matr_nn);                           /* release temporary memory */
   free(nn_pos);
   free(matr_pp);
   free(pp_pos);

   return stored_matr_elem;

}  /* End: function  pn_nondiag_SD_elem() */

    /*
    ** The function
    **           pn_nondiag_id_SD_elem()
    ** calculates all nondiagonal identical particle matrix elements 
    **      <SD(N)g.n',SD(Z)g.p'|OP| SD(Z)g,p,SD(N)g,n> 
    ** within one proton/neutron group only. Pointers to all two-particle
    ** matrix elements for particle type = PROTON(NEUTRON) are stored in 
    ** structure id_diag[part_type] by the formula
    **        op(k,l) = ((num - k) * k)/2 -1 +l
    ** where num = 2 * gr.num_orb[type] - 3. The nondiag matrix elements
    ** are  saved in memory pointed to by gr_basis->mem_free.
    ** The function returns -1 if not enough memory to store matrix
    ** elements is available.
    ** Otherwise it returns number of matrix elements stored.
    */

static int pn_nondiag_id_SD_elem(int op_type, int part_type, GR_BASIS *gr_basis,
                                                      GROUP *gr, int *pos)
{
  int
                num_part_1, num, numSD, numSD_1, currSD, num_stored, number, k, l, 
                kl_phase, phase, count, *k_ptr, *l_ptr, search_num, max_elem;
  UL
                *listSD, *sd, sd_k, sd_kl, new_sd, low, high;
  STORE
                *matr_ptr, *curr_ptr;
  ID_INT
                 **table, **tab_k, *ij_ptr;
  MASK
                *mask_ptr;

   switch(op_type) {
      case VEFF: table = gr_basis->op_int.id_nondiag[part_type];  /* <|H(op_type)|>*/
                 break;
      case ANG:  table = gr_basis->op_ang.id_nondiag[part_type]; /* <|J**2(op_type)|> */ 
                 break;
      default:   printf("\n\nError in function pn_nondiag_id_elem():");
                 printf("\nIncorrect two-particle matrix element table");
                 printf(" op_type = %d -- should be VEFF or ANG\n", op_type);
                 exit(1);
   } 
  num_part_1 = gr_basis->part[part_type] - 1;
  numSD      = gr->numSD[part_type];                      /* initialization */
  numSD_1    = numSD - 1;
  num        = (gr_basis->m_orb[part_type] << 1) - 3;
  listSD     = gr->SD[part_type];
  sd         = listSD;
  matr_ptr   = (STORE *)gr_basis->mem_free;
  num_stored = 0;
  max_elem   = gr_basis->mem_nondiag;
  *(pos++)   = 0;    /* pos ptr for matr elem in matr_ptr[] */

  for(currSD = 0; currSD < numSD_1; sd++, currSD++) {
    number   = 0;                                      /* initialization */
    curr_ptr = matr_ptr;

    k_ptr  = gr_basis->occ[part_type] + currSD * gr_basis->part[part_type];
    for(k = 0; k < num_part_1; k++, k_ptr++) { 
      tab_k    = table + (((num - (*k_ptr)) * (*k_ptr)) >> 1) -1;
      sd_k     = (*sd) ^ (UL_ONE << (*k_ptr));
      l_ptr    = k_ptr + 1;
      kl_phase = +1;
      for(l = k + 1; l <= num_part_1; l++, l_ptr++, kl_phase = -kl_phase) {
	if((ij_ptr = *(tab_k+(*l_ptr))) == NULL_PTR) continue; 
	sd_kl = sd_k ^ (UL_ONE << *l_ptr); 
	search_num = currSD;
	do  {                                               /* run through all ij-pair */
	  for(; sd_kl & ij_ptr->one; ij_ptr++);  
	  if(!ij_ptr->one) break;                 /* no more ij-pairs */
	  new_sd = sd_kl ^ ij_ptr->one; /* contr found, the new SD-state */
	  
             /* 
	     ** check that |new_sd> has allowed number of particles in all
             ** j-orbits - only active if 0 < allowed number < 2*j + 1
             */
	  
	  if((count = gr_basis->part_lim[part_type].num) > 0) {    /* check masking */
	    mask_ptr = gr_basis->part_lim[part_type].table;
	    do {
	      high = new_sd & mask_ptr->mask;
	      phase = 0;
	      while(high) {
		high &= high - 1;
		phase++;
	      }
	      if((phase < mask_ptr->min) || (phase > mask_ptr->max)) break;
	      mask_ptr++;
	    } while(--count);
	    if(count > 0) continue;
	  } /* end check masking */
	  low   = new_sd & ij_ptr->two;                   /* permutation phase. */
	  phase = kl_phase;
	  while(low) {
            low &= low - 1;
            phase = -phase;
	  }
          low  = search_num +1;                               /* binary search */
          high = numSD;
	  while(1)  {
            search_num = (low + high)>> 1;
            if(new_sd < listSD[search_num])        high = search_num - 1;
            else if(new_sd > listSD[search_num]) low  = search_num + 1;
            else                                   break;
	  }

          if((max_elem--) == 0) return -1;                    /* no more memory */

	  curr_ptr->value     = phase * ij_ptr->val;
	  (curr_ptr++)->final = search_num;
	  number++;
	} while((++ij_ptr)->one); /* all two-body matrix elements for given (kl) */
      } /* end l-loop */
    } /* end k-loop */

    max_elem += number;   /* update limit for stored matrix elements */

         /* 
	 ** sort and compress nondiag matrix elements for 
	 ** each |init_SD> after increasing |final_SD>
	 */

    if(number > 1) {
      qsort(matr_ptr, (UL) number, sizeof(STORE),
	    (int(*)(const void *, const void *)) sort_nondiag);
      number = compress_nondiag_SD_elem(number, matr_ptr); 
    }
    if((max_elem -= (number + 1)) < 1) return -1;                 /* no more memory */

	 /* Each sub-block of matrix elements is terminated by an extra element */

    matr_ptr += number;            /* update pointer to next free element */
    matr_ptr->value     = 0.0;
    (matr_ptr++)->final = -1;
    number++;

    *(pos++)    = number;
    num_stored += number;
  }  /* end loop through all |SD(part_type)> */

      /*
      ** no nondiagonal matrix elements from the last
      ** |SD(part_type>, onlu a termination elements 
      */
 
  matr_ptr->value     = 0.0;
  (matr_ptr++)->final = -1;
  num_stored++;

  return    num_stored;

} /* End: function pn_nondiag_id_contribution() */

    /*
    ** The function                         
    **        int sort_nondiag()                  
    ** is a utility function for the library function qsort() in order to
    ** sort nondiagSD matrix elements of type STORE store[]
    ** after increasing store[].final.
    */
static int sort_nondiag(const STORE *one, const STORE *two)
{
  if(one->final > two->final)       return +1;
  else  if(one->final < two->final) return -1;
  else                              return  0;
} /* End: function sort_nondiag() */

      /*
      ** The function 
      **      compress_nondiag_SD_elem()
      ** takes num calculated nondiag SD elem which are sorted
      ** after increasing store[].final and add together all
      ** elements  with same store[].final.
      ** Note: The function must have num > 1;
      ** The final number of elements is returned.
      */

static int compress_nondiag_SD_elem(int num, STORE *store)
{
   int     k, l;
   STORE   *ptr1, *ptr2, *ptr3;

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
   
} /* End: function compress_nondiag_SD_elem() */

     /*
     ** The function
     **           pn_nondiag_pn_SD_elem_A()
     ** calculates all nondiagonal pn matrix elements 
     **  <SD(Z)g,p';SD(N)g,n'|OP(pp) + OP(nn) + OP(pn)|SD(Z)g,p,SD(N)g,n>
     ** within a single group g. The two-particle matrix elements are
     ** specified through the parameter op_type.
     ** The calculation is performed in following steps. For each 
     ** init |SD(Z),SD(N)>:
     **   1. Convert previous calculated nondiagonal pp-matrix elements
     **      into pn format and store them in memory pointed to by 
     **      gr_basis->mem_free. Update gr_basis->mem_free.
     **   2. Convert previous calculated nondiagonal nn-matrix elements
     **      into pn format and store them in memory pointed to by 
     **      gr_basis->mem_free. Update gr_basis->mem_free.
     **   3. Calculate nondiagonal pn-matrix elements and store them
     **      in memory pointed to by gr_basis->mem_free.
     **   4. If the remaining free memory is too small for the matrix
     **      elements the function terminates and return -1.
     **   5. Sort after final |SD(Z),SD(N)> and compress all nondiag
     **      matrix elements.
     **   6. The function returns -1 if not enough memory is available 
     **      to store matrix elements.
     **      Otherwise it returns number of matrix elements stored.
     */   

static int pn_nondiag_pn_SD_elem_A(int op_type, GR_BASIS *gr_basis, GROUP *gr,
                             int *pp_pos, STORE *matr_pp, STORE *matr_nn) 
{
  int       
             p_num, n_num, numSD_Z, numSD_N, max_elem, stored_matr, n_lim, number,   
             k, m, *k_ptr, l, *l_ptr;
  UL
             *sdZ, *sdN, sdZ_k, sdN_l, high, low;
   PN_INT
              **table, **tab_p, *ij_ptr;
   MASK
             *mask_ptr;
   RESULT
             result[2];
   STORE
             *stored_ptr, *pn_ptr, *pp_ptr, *nn_ptr;

  switch(op_type) {
      case VEFF: table = gr_basis->op_int.pn_nondiag[0];         /* <|H(pn)|> */
                 break;
      case ANG:  table = gr_basis->op_ang.pn_nondiag[0];         /* <|J**2(pn)|>*/ 
                 break;
      default:   printf("\n\nError in function pn_nondiag_pn_SD_elem_A():");
                 printf("\nIncorrect two-particle matrix element table");
                 printf(" type = %d -- should be 1 - 2\n", op_type);
                 exit(1);
  } 
  numSD_Z = gr->numSD[0];
  numSD_N = gr->numSD[1];
   
              /* run through all init |SD(Z), SD(N)> */

  max_elem    = gr_basis->mem_nondiag;           /* maximum nondiag elements to store */
  stored_ptr  = (STORE *)gr_basis->mem_free;
  pn_ptr      = stored_ptr;
  stored_matr = 0;

             /* run through all init |SD(Z), SD(N)> */

  sdZ = gr->SD[0];

  for(p_num = 0; p_num < numSD_Z; p_num++, sdZ++) {
    matr_pp += pp_pos[p_num]; 
    sdN      = gr->SD[1];

    nn_ptr   = matr_nn;
    n_lim    = (p_num < (numSD_Z -1)) ? numSD_N : numSD_N -1;
    for(n_num = 0; n_num < n_lim; n_num++, sdN++) { 
      pn_ptr = stored_ptr;              /* initialization for each init |SD(Z),SD(N)> */
      number = 0;
      k      = p_num * numSD_N;
      while(nn_ptr->final != -1) {        /* add in and modify non-diag  nn-matr elem */

	if((--max_elem) == 0) return -1;        /* no more memory */

	pn_ptr->value = nn_ptr->value;
	(pn_ptr++)->final = (nn_ptr++)->final + k;
        number++;
      }
      nn_ptr++;              /* points to next nn-group */

      pp_ptr = matr_pp;              /* add in and modify non-diag  pp-matr elem */
      while(pp_ptr->final != -1) {

	if((--max_elem) == 0) return -1;        /* no more memory */

	pn_ptr->value     = pp_ptr->value;
	(pn_ptr++)->final = (pp_ptr++)->final * numSD_N + n_num;
        number++;
      }
      k     = gr_basis->part[0];
      k_ptr = gr_basis->occ[0] + p_num * k;
      do {                                      /* proton-particle loop */
	tab_p  = table + (*(k_ptr) *  gr_basis->m_orb[1]);
	sdZ_k = (*sdZ) ^ (UL_ONE << (*(k_ptr++)));
	l      = gr_basis->part[1];
	l_ptr  = gr_basis->occ[1] + n_num * l;
	do   {                                 /* neutron particle loop */
	  sdN_l = (*sdN) ^ (UL_ONE << (*l_ptr));
          if((ij_ptr = *(tab_p + (*(l_ptr++)))) == NULL) continue;

	  result[0].sd     = 0;                         /* initialization */
	  result[1].sd     = 0;
	  
	  do   {                    /* loop through all two-particle matrix elements */
                  /* 
                  ** Pauli principle: 
                  **  proton orbit i and  neutron orbit j must both be empty. 
                  */ 

	    while(   ij_ptr->one[0] 
                  && (   (ij_ptr->one[0] & sdZ_k)
                      || (ij_ptr->one[1] & sdN_l))) ij_ptr++;

            if(!ij_ptr->one[0]) break;                 /* end of ij-pairs */

            if(result[0].sd   != (sdZ_k ^ ij_ptr->one[0])) {    /* proton case */
               result[0].sd    = (sdZ_k ^ ij_ptr->one[0]);   /* new proton SD-state */
               result[0].phase = 0;
               result[0].mask  = YES;

                  /* mask analysis:
	          ** check that |new_sdZ> has allowed number of protons  in all
                  ** j-orbits - only active if 0 < allowed number < 2*j + 1
                  */

               if((low = gr_basis->part_lim[0].num) > 0) {    /* check masking */
                  mask_ptr = gr_basis->part_lim[0].table;
                  do {
                     high = result[0].sd & mask_ptr->mask;
                     m = 0;
                     while(high) {
                        high &= high - 1;
                        m++;
                     }
                     if((m < mask_ptr->min) || (m > mask_ptr->max)) break;
                     mask_ptr++;
                  } while(--low);
                  result[0].mask = (low == 0) ? YES : NO;
               }                                /* end check proton masking */
            }   /* end check proton mask */
            if(result[0].mask == NO) continue;

	    if(result[1].sd   != (sdN_l ^ ij_ptr->one[1])) { /* neutron case */
	      result[1].sd    = (sdN_l ^ ij_ptr->one[1]);   /* new neutron SD-state */
	      result[1].phase = 0;
	      result[1].mask  = YES;  

                  /* mask analysis:
	          ** check that |new_sdN> has allowed number of neutrons  in all
                  ** j-orbits - only active if 0 < allowed number < 2*j + 1
                  */

	      if((low = gr_basis->part_lim[1].num) > 0) {
		mask_ptr = gr_basis->part_lim[1].table;
		do {
		  high = result[1].sd & mask_ptr->mask;
		  m = 0;
		  while(high) {
		    high &= high - 1;
		    m++;
		  }
		  if((m < mask_ptr->min) || (m > mask_ptr->max)) break;
		  mask_ptr++;
		} while(--low);
		result[1].mask = (low == 0) ? YES : NO;
	      } /* end check neutron mask */
            } 
            if(result[1].mask == NO) continue;

            if(result[0].phase == 0) {
	      result[0].phase = +1;              /* proton  permutation phase. */
	      low = result[0].sd & ij_ptr->two[0];
	      while(low) {
		low &= low - 1;
		result[0].phase = -result[0].phase;
	      }
	      low  =  0;           
	      high =  numSD_Z; 
	      while (1) {
		result[0].sd_num = (low + high) >> 1;
		if(result[0].sd < gr->SD[0][result[0].sd_num])
		                                        high = result[0].sd_num-1;
		else if(result[0].sd >  gr->SD[0][result[0].sd_num])
                                                        low  = result[0].sd_num+1;
		else                                  break;
	      }
            }
	    if(result[1].phase == 0) {
	      result[1].phase = +1;              /* neutron  permutation phase. */
	      low = result[1].sd & ij_ptr->two[1];
	      while(low) {
		low &= low - 1;
		result[1].phase = -result[1].phase;
	      }
	      low  = 0;
	      high = numSD_N; 
	      while (1) {
		result[1].sd_num= (low + high) >> 1;
		if(result[1].sd < gr->SD[1][result[1].sd_num]) 
                                                     high = result[1].sd_num - 1;
		else if(result[1].sd >gr->SD[1][result[1].sd_num]) 
                                                     low  = result[1].sd_num + 1;
		else                               break;
	      }
            }
                	    /* store pn-matrix elements */

	    if((--max_elem) == 0) return -1;               /* no more memory */

	    pn_ptr->value     = result[0].phase * result[1].phase * ij_ptr->val;
            (pn_ptr++)->final = result[0].sd_num * numSD_N  + result[1].sd_num;
            number++;
	  } while((++ij_ptr)->one[0]); /* loop through all matr.elem for fixed kl-pairs */
	} while(--l);   /* end l-loop for neutron particles */
      } while(--k);      /* end k-loop for proton particles */

      max_elem += number;   /* update limit for stored matrix elements */

         /* 
	 ** sort and compress nondiag matrix elements for 
	 ** each init |SD(Z),SD(N)> after increasing final |SD(Z),SD(N)>
	 */

      if(number > 1) {
	qsort(stored_ptr, (UL) number, sizeof(STORE),
	      (int(*)(const void *, const void *)) sort_nondiag);
	number = compress_nondiag_SD_elem(number, stored_ptr); 
      }

     if((max_elem -= (number + 1)) <= 1) return -1 ;  /* no more memory */

      stored_ptr           += number++;    /* update pointer to next free element */
      stored_ptr->value     = 0.0;         /* each sub-block terminates by an extra element */
      (stored_ptr++)->final = -1;
      stored_matr          += number;      /* number of elements stored */

    } /* end n_num */
  } /* end p_num */

	     /* add an extra stop element for the last init |SD(Z),SD(N)> */

  stored_ptr->value     = 0.0;  
  (stored_ptr++)->final = -1;
  stored_matr++;

   return stored_matr;               /* successfull storage of matr elem */

}  /* End: function pn_nondiag_pn_SD_elem_A() */

   /*
   ** The function 
   **         pn_write_SD_matr_elem()
   ** opens the file "file_name" and write/append diagonal
   ** and nondiagonal matrix elements 
   **       <SD(Z)g',p: SD(N)g',n|OP|SD(Z)g,p: SD(N)g,n>
   ** for (g,p), (g,n) --> (g',p'),(g',n')
   ** if(header.store_nondiag == -1) no nondiagonal matrix
   ** elements are calculated due to lack of memory.
   ** If(new_file == YES) a new file is created.
   */

static void pn_write_SD_matr_elem(char *file_name, int new_file,HEADER header,char * mem_ptr) 
{
/******/
  int
            k;
  double
            *double_ptr;
  STORE
            *store_ptr;
/****/
  UL
            num_char;
  FILE
           *file_ptr;

  if((file_ptr = fopen(file_name,(new_file == YES) ? "wb":"ab")) == NULL) {
    printf("\n\nError in function  pn_write_SD_matr_elem():");
    printf("\nNot allowed to %s the file %s\n", (new_file == YES) ? 
                                               "create":"append", file_name);
    exit(1);
  }
               /* write header elements */

  if( fwrite((const void *)&header,(size_t) sizeof(HEADER), 1, file_ptr) != 1) { 
    printf("\n\nError in function pn_write_SD_matr_elem():");
    printf("\nIn writing the header element between");
    printf("\ninit_gr = %d to final_gr = %d\n", header.init_gr, header.final_gr);
    exit(1);
  }

  /**************  test output  *******
  printf("\n\nFile_name = %s", file_name);
  printf("\nHeader: init_gr = %d  final_gr = %d  stored_diag = %d  stored_nondiag = %d\n",
       header.init_gr, header.final_gr, header.stored_diag, header.stored_nondiag);

  if(header.stored_diag > 0) {

    double_ptr = (double *)mem_ptr;

    printf("\ndiagonal matrix elements:\n");
    for(k = 0; k < header.stored_diag; k++) {
      printf(" %8.4f ",*(double_ptr++));
	   if(((k + 1)/5)*5 == (k + 1)) printf("\n");
    }
  }
    store_ptr = (STORE *)(mem_ptr + header.stored_diag * sizeof(double));

  printf("\nnondiagonal matrix elements:\n");
  for(k = 0; k < header.stored_nondiag; k++, store_ptr++) {
    printf(" %3d  %8.4f ",store_ptr->final, store_ptr->value);
	   if(((k + 1)/5)*5 == (k + 1)) printf("\n");
  }
****************  end output  ***********/

    /* write diag and nondiagSD matrix elements as a string of bytes */

  num_char =  (header.stored_diag * sizeof(double))/sizeof(char);
  if(header.stored_nondiag > 0) { 
    num_char += (header.stored_nondiag * sizeof(STORE))/sizeof(char);
  }
  if(num_char > 0) {
    if( fwrite((const void *)mem_ptr,(size_t) sizeof(char), (size_t) num_char, file_ptr)
                                                                          != num_char) { 
      printf("\n\nError in function pn_write_SD_matr_elem():");
      printf("\nIn writing %d bytes  for matrix elements between", num_char);
      printf("\ninit_gr = %d to final_gr = %d\n", header.init_gr, header.final_gr);
      exit(1);
    }
  }
  fclose(file_ptr);       /* close file */     

} /* End: function pn_write_pn_SD_matr_elem() */

     /*
     ** The function
     **           pn_nondiag_pn_SD_elem_B()
     ** calculates all nondiagonal pn-matrix elements 
     **  <SD(Z)g',p';SD(N)g',n'|OP(pp) + OP(nn) + OP(pn)|SD(Z)g,p,SD(N)g,n>
     ** between two groups g -> g'. The two-particle matrix elements are
     ** specified through the parameter op_type.
     ** The calculation is performed in following steps. For each 
     ** init |SD(Z),SD(N)>:
     **   1. Calculate nondiagonal pn-matrix elements and store them
     **      in memory pointed to by gr_basis->mem_free.
     **   2. If the free memory is too small for the matrix
     **      elements the function terminates and return -1.
     **   3. Sort after final |SD(Z),SD(N)> and compress all nondiag
     **      matrix elements.
     **   4. The function returns number of matrix elements stored.
     */   

static int pn_nondiag_pn_SD_elem_B(int op_type, GR_BASIS *gr_basis, GROUP *init_gr,
                                                                    GROUP *final_gr)
{
  int       
             p_num, n_num, init_numSD_Z, init_numSD_N, final_numSD_Z, final_numSD_N, 
             stored_matr_elem, max_elem, number, k, m, *k_ptr, l,
             *l_ptr, step_m_p;
  UL
             *sdZ, *sdN, sdZ_k, sdN_l, high, low;
   PN_INT
             **table, **tab_p, *ij_ptr;
   MASK
             *mask_ptr;
   RESULT
             result[2];
   STORE
             *stored_ptr, *pn_ptr;

      /* calculate changes in M_Z and PAR_Z from init to final group */

   if(gr_basis->parZ) {
     step_m_p =  (init_gr->m[0] - final_gr->m[0]) >> 1;
   }
   else {
     step_m_p =  (init_gr->m[0] - final_gr->m[0]) 
            + (((init_gr->par[0] * final_gr->par[0]) == +1) ? 0 : 1);
   }
  switch(op_type) {
      case VEFF: table = gr_basis->op_int.pn_nondiag[step_m_p];   /* <|H(pn)|> */
                 break;
      case ANG:  table = gr_basis->op_ang.pn_nondiag[step_m_p];   /* <|J**2(pn)|>*/ 
                 break;
      default:   printf("\n\nError in function pn_nondiag_pn_SD_elem_B():");
                 printf("\nIncorrect two-particle matrix element table");
                 printf(" type = %d -- should be 1 - 2\n", op_type);
                 exit(1);
  } 

  init_numSD_Z  = init_gr->numSD[0];
  init_numSD_N  = init_gr->numSD[1];
  final_numSD_Z = final_gr->numSD[0];
  final_numSD_N = final_gr->numSD[1];
   
  max_elem         = gr_basis->mem_nondiag;                     /* upper limit */
  stored_ptr       = (STORE *)gr_basis->mem_pos;
  pn_ptr           = stored_ptr;
  stored_matr_elem = 0;

             /* run through all init |SD(Z), SD(N)> */

  sdZ = init_gr->SD[0];

  for(p_num = 0; p_num < init_numSD_Z; p_num++, sdZ++) {
    sdN      = init_gr->SD[1];
    for(n_num = 0; n_num < init_numSD_N; n_num++, sdN++) { 
      pn_ptr = stored_ptr;     /* initialization for each init |SD(Z),SD(N)> */
      number = 0;
      k     = gr_basis->part[0];
      k_ptr = gr_basis->occ[0] + p_num * k;
      do {                                                  /* proton-particle loop */
	tab_p  = table + (*(k_ptr) *  gr_basis->m_orb[1]);
	sdZ_k = (*sdZ) ^ (UL_ONE << (*(k_ptr++)));
	l      = gr_basis->part[1];
	l_ptr  = gr_basis->occ[1] + n_num * l;
	do   {                                 /* neutron particle loop */
	  sdN_l = (*sdN) ^ (UL_ONE << (*l_ptr));
          if((ij_ptr = *(tab_p + (*(l_ptr++)))) == NULL) continue;

	  result[0].sd     = 0;                         /* initialization */
	  result[1].sd     = 0;

	  do   {                    /* loop through all two-particle matrix elements */
                  /* 
                  ** Pauli principle: 
                  **  proton orbit i and  neutron orbit j must both be empty. 
                  */ 

            while(   ij_ptr->one[0] 
                  && (   (ij_ptr->one[0] & sdZ_k)
                      || (ij_ptr->one[1] & sdN_l))) ij_ptr++;

            if(!ij_ptr->one[0]) break;                 /* end of ij-pairs */

            if(result[0].sd   != (sdZ_k ^ ij_ptr->one[0])) {    /* proton case */
               result[0].sd    = (sdZ_k ^ ij_ptr->one[0]);   /* new proton SD-state */
               result[0].phase = 0;          /* initialization of new phase */
               result[0].mask  = YES;

                  /* mask analysis:
	          ** check that |new_sdZ> has allowed number of protons  in all
                  ** j-orbits - only active if 0 < allowed number < 2*j + 1
                  */

               if((low = gr_basis->part_lim[0].num) > 0) {    /* check masking */
                  mask_ptr = gr_basis->part_lim[0].table;
                  do {
                     high = result[0].sd & mask_ptr->mask;
                     m = 0;
                     while(high) {
                        high &= high - 1;
                        m++;
                     }
                     if((m < mask_ptr->min) || (m > mask_ptr->max)) break;
                     mask_ptr++;
                  } while(--low);
                  result[0].mask = (low == 0) ? YES : NO;
               }                                /* end check proton masking */
            }   /* end check proton mask */
            if(result[0].mask == NO) continue;
    
             if(result[1].sd   != (sdN_l ^ ij_ptr->one[1])) { /* neutron case */
               result[1].sd    = (sdN_l ^ ij_ptr->one[1]);   /* new neutron SD-state */
               result[1].phase = 0;             /* initialization of new phase */
               result[1].mask  = YES;  

                  /* mask analysis:
	          ** check that |new_sdN> has allowed number of neutrons  in all
                  ** j-orbits - only active if 0 < allowed number < 2*j + 1
                  */

               if((low = gr_basis->part_lim[1].num) > 0) {
                  mask_ptr = gr_basis->part_lim[1].table;
                  do {
                     high = result[1].sd & mask_ptr->mask;
                     m = 0;
                     while(high) {
                        high &= high - 1;
                        m++;
                     }
                     if((m < mask_ptr->min) || (m > mask_ptr->max)) break;
                     mask_ptr++;
                  } while(--low);
                  result[1].mask = (low == 0) ? YES : NO;
               } /* end check neutron mask */
            } 
            if(result[1].mask == NO) continue;

            if(result[0].phase == 0) {
	      result[0].phase = +1;              /* proton  permutation phase. */
	      low = result[0].sd & ij_ptr->two[0];
	      while(low) {
		low &= low - 1;
		result[0].phase = -result[0].phase;
	      }
	      low  = 0;           
	      high =  final_numSD_Z; 
	      while (1) {
		result[0].sd_num = (low + high) >> 1;
		if(result[0].sd < final_gr->SD[0][result[0].sd_num])
		                                        high = result[0].sd_num-1;
		else if(result[0].sd > final_gr->SD[0][result[0].sd_num])
                                                        low  = result[0].sd_num+1;
		else                                  break;
	      }
            }
	    if(result[1].phase == 0) {
	      result[1].phase = +1;              /* neutron  permutation phase. */
	      low = result[1].sd & ij_ptr->two[1];
	      while(low) {
		low &= low - 1;
		result[1].phase = -result[1].phase;
	      }
	      low  = 0;
	      high = final_numSD_N; 
	      while (1) {
		result[1].sd_num= (low + high) >> 1;
		if(result[1].sd < final_gr->SD[1][result[1].sd_num]) 
                                                     high = result[1].sd_num - 1;
		else if(result[1].sd > final_gr->SD[1][result[1].sd_num]) 
                                                     low  = result[1].sd_num + 1;
		else                               break;
	      }
            }

            if((--max_elem) == 0) return -1;                    /* no nore memory */

	    pn_ptr->value     = result[0].phase * result[1].phase * ij_ptr->val;
            (pn_ptr++)->final = result[0].sd_num * final_numSD_N  + result[1].sd_num;
            number++;
            if((--max_elem) <= 0) return -1;   /* not enough memory, stop and return */

         } while((++ij_ptr)->one[0]); /* loop through all matr.elem for fixed kl-pairs */
	} while(--l);   /* end l-loop for neutron particles */
      } while(--k);      /* end k-loop for proton particles */

      max_elem += number;   /* update limit for stored matrix elements */

         /* 
	 ** sort and compress nondiag matrix elements for 
	 ** each init |SD(Z),SD(N)> after increasing final |SD(Z),SD(N)>
	 */

      if(number > 1) {
	qsort(stored_ptr, (UL) number, sizeof(STORE),
	      (int(*)(const void *, const void *)) sort_nondiag);
	number = compress_nondiag_SD_elem(number, stored_ptr); 
      }

      if((max_elem -= (number + 1)) < 1) return -1;    /* no more memory */

         /*
	 ** Each sub-block of matrix elements
	 ** is terminated by an extra element
	 */

      stored_ptr           += number++;            /* update pointer to next free element */
      stored_ptr->value     = 0.0;
      (stored_ptr++)->final = -1;
      stored_matr_elem     += number; /* number of elements stored */

    } /* end n_num */
  } /* end p_num */

   return stored_matr_elem;               /* successfull storage of matr elem */

}  /* End: function pn_nondiag_pn_SD_elem_B() */
