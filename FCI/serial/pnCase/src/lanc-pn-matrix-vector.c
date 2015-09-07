
/*******************  The module lanc-pn-matrix-vector.c  ******************/

#include "shell.h"

   /*
   ** The entrance function
   **         matrix-vector-multiplication()     
   ** takes as input an orthonormalized Lanczos vector
   ** and perform the process
   **            OP * |vec.one[]> = |vec_two[]>
   ** OP is a two-particle Hermetian operator and its matrix
   ** elements is stored in the table MATR_OP op[].
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

static HEADER  pn_read_SD_matr_elem(char *, FILE *, int, int, char *);
      /*
      ** reads from open file precalculated diagonal and nondiagonal
      ** matrix elements 
      **       <SD(Z)g',p: SD(N)g',n|OP|SD(Z)g,p: SD(N)g,n>
      ** for (g,p), (g,n) --> (g',p'),(g',n')
      */

static void pn_diag_SD_contr(int, double *, double *, double *);
      /*
      ** calculates and store the contribution to 
      **       H * |init_vec> ---> |final_vec>
      ** from the diagonal terms in a single group
      */

static void pn_nondiag_SD_contr(int, double *, double *, double *, double *, STORE *);
      /*
      ** calculates and store the contribution to 
      **       H * |init_vec> ---> |final_vec>
      ** from the nondiagonal terms where |init_vec> and |final_vec>
      ** refer to the same group or to a pair of different groups.
      */

static void pn_occupied_orbits(int, GR_BASIS *, GROUP *);
      /*
      ** calculates and stores orbit numbers of all occupied orbits  in all 
      ** |SD(Z)> (part_type = 0) and |SD(N)> (part_type = 1) in present group gr[] 
      */

static void pn_nondiag_pp_SD_elem(int, GR_BASIS *, GROUP *, double *, double *);
      /*
      ** calculates all nondiagonal proton-proton matrix elements 
      **      <SD(N)g.n',SD(Z)g.p'|OP(pp)| SD(Z)g,p,SD(N)g,n> 
      ** within one proton/neutron group only. Pointers to all two-particle
      ** matrix elements are stored in structure gr_basis->op_int.id_table[0]
      ** by the formula
      **        op(k,l) = ((num - k) * k)/2 -1 +l
      ** where num = 2 * gr.num_orb[t0] - 3.
      */

static void pn_nondiag_nn_SD_elem(int, GR_BASIS *, GROUP *, double *, double *);
      /*
      ** calculates all nondiagonal neutron-neutronmatrix elements 
      **      <SD(N)g.n',SD(Z)g.p'|OP(nn)| SD(Z)g,p,SD(N)g,n> 
      ** within one proton/neutron group only. Pointers to all two-particle
      ** matrix elements are stored in structure gr_basis->op_int.id_table[1]
      ** by the formula
      **        op(k,l) = ((num - k) * k)/2 -1 +l
      ** where num = 2 * gr.num_orb[1] - 3.
      */

static void pn_nondiag_pn_SD_elem_A(int, GR_BASIS *, GROUP *, double *, double *);
      /*
      ** calculates all nondiagonal pn-matrix elements 
      **  <SD(Z)g,p';SD(N)g,n'|OP(pn)|SD(Z)g,p,SD(N)g,n>
      ** within a single group The two-particle matrix elements are
      ** specified through the parameter op_type.
      */   

static void pn_nondiag_pn_SD_elem_B(int, GR_BASIS *, GROUP *, GROUP *, double *,
                                                    double *, double *, double *);
      /*
      ** calculates all nondiagonal pn-matrix elements 
      **  <SD(Z)g',p';SD(N)g',n'|OP(pp) + OP(nn) + OP(pn)|SD(Z)g,p,SD(N)g,n>
      ** between two groups g -> g'. The two-particle matrix elements are
      ** specified through the parameter op_type.
      */   

               /****   end local function declarations ****/

               /**** The function definitions  ****/ 
 
      /*
      ** The entrance function 
      **        matric_vector_multiplication()
      ** takes as input an orthonormalized Lanczos vector
      ** and perform the process
      **            OP * |vec.one[]> = |vec_two[]>
      ** OP is a two-particle Hermetian operator and its matrix
      ** elements is specified by the parameter op_type and stored in gr_basis
      */

void matrix_vector_multiplication(GR_BASIS *gr_basis, GROUP *init_gr, int op_type, LANC vec)
{
  char
               *func = {"matrix_vector_multiplication: "},
               *file_name;
  int
               dim, init_gr_num, final_gr_num, max_delM_Z;
  int
                k; 
  double 
               *double_ptr;
  STORE
               *store_ptr;
  GROUP
               *final_gr;
  HEADER
               header;
  FILE
               *file_ptr;

/***************  test output  ***************
printf("\n\n\nMatrix-vector calculation for %s interaction: \n\n\n",
                                         (op_type ? "ANG" : "VEFF" ));
printf("\n Initial Lanczo vector:\n\n");
    for(k = 0; k <gr_basis->tot_numSD ; k++) {
	printf("% 7.4f ", vec.one[k]);
        if((((k + 1)/10)*10) == ( k + 1)) printf("\n");
      }
************   end test output  ***/


           /* maximum change in M_Z between initial and final groups */

   max_delM_Z = 2 * gr_basis->mbas[0]->m;

        /* memory to store occupied orbits in all |SD(Z)> and SD(N)> in a group */

  gr_basis->occ[0] = MALLOC(gr_basis->maxSD[0] * gr_basis->part[0], int,func,"p_occ[]");
  gr_basis->occ[1] = MALLOC(gr_basis->maxSD[1] * gr_basis->part[1], int,func,"n_occ[]");

           /* memory to store both diagSD[] and nondiagSD[] matr elem  - unit char */

  dim                = ((gr_basis->maxSD[0] * gr_basis->maxSD[1]) * sizeof(double) 
                       + gr_basis->mem_nondiag * sizeof(STORE))/sizeof(char);
  gr_basis->mem_pos  = MALLOC(dim, char, func,"mem_ptr[]");

          /* open file containing precalculated matrix elements */

  file_name = op_type ? gr_basis->store_ang : gr_basis->store_veff;

  if((file_ptr = fopen(file_name,"rb")) == NULL) {
    printf("\n\nError in function matrix_vector_multiplication():");
    printf("\nNot allowed to open the file %s\n", file_name);
    exit(1);
  }
  rewind(file_ptr);

       /* loop through all init group basis states |SD(Z),SD(N)> */ 

  for(init_gr_num = 0; init_gr_num < gr_basis->num_gr; init_gr_num++, init_gr++) {

    if((init_gr->numSD[0] * init_gr->numSD[1]) == 0)   continue;  /* no basis states */

           /* read precalculated matrix elements from file for init_gr */

    header = pn_read_SD_matr_elem(file_name, file_ptr, init_gr_num, init_gr_num,
                                                            gr_basis->mem_pos);
               /* diagonal contribution to |lanc.vec.two[]> */
    
    double_ptr = (double *)gr_basis->mem_pos;
    pn_diag_SD_contr(init_gr->numSD[0] * init_gr->numSD[1], vec.one + init_gr->start_amp, 
		                     vec.two + init_gr->start_amp, double_ptr);

                /* nondiag contribution for final_gr = init_gr */

    if(header.stored_nondiag > 0) {

             /* add contribution using precalculated SD matrix elements */

      store_ptr = (STORE *)(gr_basis->mem_pos 
                    + (header.stored_diag * sizeof(double))/sizeof(char)); 

/************ test output  ******
      printf("\n\ninit_group = %d  final_group = %d", init_gr_num, init_gr_num);
      printf("\n Header: init_gr = %d final_gr = %d stored_diag = %d store_nondiag = %d",
	     header.init_gr, header.final_gr, header.stored_diag, header.stored_nondiag);
      printf("\ndiag_elem:\n");
      for(k = 0; k < header.stored_diag; k++) {
	printf("%10.4f", double_ptr[k]);
        if((((k + 1)/6)*6) == ( k + 1)) printf("\n");
      }
      printf("\nnondiagdiag_elem:\n");
      for(k = 0; k < header.stored_nondiag; k++) {
	printf("%4d %10.4f", store_ptr[k].final, store_ptr[k].value);
        if((((k + 1)/6)*6) == ( k + 1)) printf("\n");
      }
*****************************************/

      pn_nondiag_SD_contr(init_gr->numSD[0] * init_gr->numSD[1], 
                          vec.one + init_gr->start_amp, vec.one + init_gr->start_amp,
                          vec.two + init_gr->start_amp, vec.two + init_gr->start_amp,
                          store_ptr);
    }
    else  {
             /* 
	     ** calculate explicitly nondiagSD proton-proton
	     ** matrix elements between states in init group
             */
 
             /* find all occupied p-orbits and n_orbits  in group[k] */
 
      pn_occupied_orbits(PROTON, gr_basis, init_gr);                    /* p-orbits */
      pn_occupied_orbits(NEUTRON, gr_basis, init_gr);                   /* n-orbits */ 

      if((gr_basis ->part[PROTON] >= 2) && (init_gr->numSD[0] >= 2)) {
	pn_nondiag_pp_SD_elem(op_type, gr_basis, init_gr,
                           vec.one + init_gr->start_amp,
                           vec.two + init_gr->start_amp);
      } /* end proton case */

             /* 
	     ** calculate nondiagSD neutron-neutron matrix  
	     ** elements between states in current group
             */

      if((gr_basis ->part[NEUTRON] >= 2) && (init_gr->numSD[1] >= 2)) {
	pn_nondiag_nn_SD_elem(op_type, gr_basis, init_gr,
                           vec.one + init_gr->start_amp,
                           vec.two + init_gr->start_amp);
      } /* end neutron case */

	     /* 
	     ** calculate nondiagSD protin-neutron matrix  
	     ** elements between states in current group
             */

      pn_nondiag_pn_SD_elem_A(op_type, gr_basis, init_gr, vec.one + init_gr->start_amp,
			                                  vec.two + init_gr->start_amp);
    }

/***************  test output  ***************
printf("\nMatrix-vector calculation: \n");
printf("\n init_gr = %d  --  final_gr = %d\n", init_gr_num, final_gr_num);
printf("\n New Lanczo vector:\n\n");
    for(k = 0; k <gr_basis->tot_numSD ; k++) {
	printf("% 7.4f ", vec.two[k]);
        if((((k + 1)/10)*10) == ( k + 1)) printf("\n");
      }
************   end test output  ***/




          /*
	  ** calculate and store nondiag matrix elements of <SD()|OP|SD()> for
          **                     final_gr > init_gr 
          */

    for(final_gr_num = init_gr_num + 1, final_gr = init_gr + 1;
                 final_gr_num  <  gr_basis->num_gr; final_gr_num++, final_gr++) {

       if(   ((final_gr->numSD[0] * final_gr->numSD[1]) == 0)    /* no basis states */
          || ((init_gr->m[0] - final_gr->m[0]) > max_delM_Z))   /* too large del_M_Z */
	 continue;

           /* read precalculated matrix elements from file for init_gr --> final_gr */

      header = pn_read_SD_matr_elem(file_name, file_ptr, init_gr_num, final_gr_num,
                                                            gr_basis->mem_pos);
  
                /* nondiag contribution for final_gr > init_gr */

      if(header.stored_nondiag > 0) {

             /* add contribution using precalculated SD matrix elements */

	store_ptr = (STORE *)gr_basis->mem_pos; 

	pn_nondiag_SD_contr(init_gr->numSD[0] * init_gr->numSD[1], 
                          vec.one + init_gr->start_amp, vec.one + final_gr->start_amp,
			  vec.two + init_gr->start_amp, vec.two + final_gr->start_amp,
			  store_ptr); 

/************ test output  ******
      printf("\n\ninit_group = %d  final_group = %d", init_gr_num, final_gr_num);
      printf("\n Header: init_gr = %d final_gr = %d stored_diag = %d store_nondiag = %d",
	     header.init_gr, header.final_gr, header.stored_diag, header.stored_nondiag);
      printf("\ndiag_elem:\n");
      for(k = 0; k < header.stored_diag; k++) {
	printf("%10.4f", double_ptr[k]);
        if((((k + 1)/6)*6) == ( k + 1)) printf("\n");
      }
      printf("\nnondiagdiag_elem:\n");
      for(k = 0; k < header.stored_nondiag; k++) {
	printf("%4d %10.4f", store_ptr[k].final, store_ptr[k].value);
        if((((k + 1)/6)*6) == ( k + 1)) printf("\n");
      }
*****************************************/
      }
      else  {
	     /* 
	     ** calculate nondiagSD protin-neutron matrix  
	     ** elements for init_gr --> final_gr
             */

             /* find all occupied p-orbits and n_orbits  in init group  */
 
	pn_occupied_orbits(PROTON, gr_basis, init_gr);                    /* p-orbits */
	pn_occupied_orbits(NEUTRON, gr_basis, init_gr);                   /* n-orbits */ 

	pn_nondiag_pn_SD_elem_B(op_type, gr_basis, init_gr, final_gr, 
                          vec.one + init_gr->start_amp, vec.one + final_gr->start_amp,
			  vec.two + init_gr->start_amp, vec.two + final_gr->start_amp);
      }

/***************  test output  ***************
printf("\nMatrix-vector calculation: \n");
printf("\n init_gr = %d  --  final_gr = %d\n", init_gr_num, final_gr_num);
printf("\n New Lanczo vector:\n\n");
    for(k = 0; k <gr_basis->tot_numSD ; k++) {
	printf("% 7.4f ", vec.two[k]);
        if((((k + 1)/10)*10) == ( k + 1)) printf("\n");
      }
************   end test output  ***/




    } /* end loop through nondiagonal groups */

  } /* end loop through diagonal groups */

  fclose(file_ptr);
  free(gr_basis->mem_pos);                 /* release local memory */ 
  free(gr_basis->occ[1]); 
  free(gr_basis->occ[0]); 

  gr_basis->mem_pos  = NULL;

/***************  test output  ***************
printf("\n\n\nEnd Matrix-vector calculation: \n\n\n");

printf("\n New Lanczo vector:\n\n");
    for(k = 0; k <gr_basis->tot_numSD ; k++) {
	printf("% 7.4f ", vec.two[k]);
        if((((k + 1)/10)*10) == ( k + 1)) printf("\n");
      }
************   end test output  ***/




} /* End: function matrix_vector_multiplication() */

   /*
   ** The function 
   **         pn_read_SD_matr_elem()
   ** reads from open file precalculated diagonal and nondiagonal
   ** matrix elements 
   **       <SD(Z)g',p: SD(N)g',n|OP|SD(Z)g,p: SD(N)g,n>
   ** for (g,p), (g,n) --> (g',p'),(g',n')
   */

static HEADER  pn_read_SD_matr_elem(char *file_name, FILE *file_ptr, int init_gr, 
                                                      int final_gr, char *mem_ptr) 
{
  int
             num;
  HEADER
             header;
  double
            *double_ptr;
  STORE
            *store_ptr;

               /* read header elements */

  if(fread((void *) &header,(size_t) sizeof(HEADER),1,file_ptr) != 1)  {
    printf("\n\nError in function pn_read_SD_matr_elem():");
    printf("\nin reading from file %s", file_name);
    printf("\nthe header element between");
    printf("\ninit_gr = %d to final_gr = %d\n", init_gr, final_gr);
    exit(1);
  }
  if((header.init_gr != init_gr) || (header.final_gr != final_gr)) {
    printf("\n\nError in function pn_read_SD_matr_elem():");
    printf("\nin reading from file %s", file_name);
    printf("\nthe wrong block of precalculated SD matrix elements");
    printf("\nshould be (init_gr, final_gr) = (%d, %d)", init_gr, final_gr);
    printf("\nread values (init_gr, final_gr) = (%d, %d)\n", 
                	                   header.init_gr, header.final_gr);
    exit(1);
  }
               /* read diagonal SD matrix elements */

  if(header.stored_diag > 0) {
    double_ptr = (double *)mem_ptr;
    if(fread((void *) double_ptr,(size_t) sizeof(double),header.stored_diag,
	                                    file_ptr) != header.stored_diag)  {
      printf("\n\nError in function pn_read_SD_matr_elem():");
      printf("\nin reading from file %s",file_name);
      printf("\nthe precalculated diagonal SD matrix elements");
      printf("\n for init group %d\n",init_gr);
      exit(1);
    }
  } /* end reading precalculated diagonal Sd matrix elements */

  if(header.stored_nondiag > 0) {
    store_ptr = (STORE *)(mem_ptr + (header.stored_diag * sizeof(double))/sizeof(char));

    if((num = fread((void *) store_ptr,(size_t) sizeof(STORE),header.stored_nondiag,
	                                    file_ptr)) != header.stored_nondiag)  {
      printf("\n\nError in function pn_read_SD_matr_elem():");
      printf("\nin reading from file %s",file_name);
      printf("\nthe precalculated nondiagonal SD matrix elements");
      exit(1);
    }
  } /* end reading precalculated nondiagonal Sd matrix elements */

  return header;

} /* End: function pn_read_pn_SD_matr_elem() */

        /*
	** The function 
	**        pn_diag_SD_contr()
	** calculates and store the contribution to 
	**       H * |init_vec> ---> |final_vec>
	** from the diagonal terms in a single group
	*/

static void pn_diag_SD_contr(int num_elem, double *in_amp, double *out_amp, double *matr)
{
  do {
    *(out_amp++) += *(matr++) * (*(in_amp++));
  } while(--num_elem);
} /* End: function pn_diag_SD_contr() */


        /*
	** The function 
	**        pn_nondiag_SD_contr()
	** calculates and store the contribution to 
	**       H * |init_vec> ---> |final_vec>
	** from the nondiagonal terms where |init_vec> and |final_vec>
	** refer to the same group or to a pair of different groups.
	*/

static void pn_nondiag_SD_contr(int init_num_elem, double *in_amp1, double *in_amp2,
                              double *out_amp1, double *out_amp2, STORE *matr)
{
  int
             k;
  double  
             *ptr;
  STORE
             *matr_ptr;

  ptr      = in_amp1;               /* initialization */
  matr_ptr = matr;
  
  for(k = 0; k < init_num_elem; k++, ptr++, matr_ptr++) {
    if(matr_ptr->final != -1) {
      do {
	out_amp2[matr_ptr->final] += (*ptr) * matr_ptr->value;
      } while((++matr_ptr)->final != -1);
    }
  } /* end loop through all |init_SD> */

  ptr      = out_amp1;               /* initialization */
  matr_ptr = matr;

  for(k = 0; k < init_num_elem; k++, ptr++,matr_ptr++) {
    if(matr_ptr->final != -1) {
      do {
	*ptr += in_amp2[matr_ptr->final] * matr_ptr->value;
      } while((++matr_ptr)->final != -1);
    }
  } /* end loop through all |init_SD> */

} /* End: function pn_nondiag_SD_contr() */

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
    **           pn_nondiag_pp_SD_elem()
    ** calculates all nondiagonal proton-proton matrix elements 
    **      <SD(N)g.n',SD(Z)g.p'|OP(pp)| SD(Z)g,p,SD(N)g,n> 
    ** within one proton/neutron group only. Pointers to all two-particle
    ** matrix elements are stored in structure gr_basis->op_int.id_diag[0]
    ** by the formula
    **        op(k,l) = ((num - k) * k)/2 -1 +l
    ** where num = 2 * gr.num_orb[t0] - 3.
    */

static void pn_nondiag_pp_SD_elem(int op_type, GR_BASIS *gr_basis, GROUP *gr,
                                             double *in_amp, double *out_amp)
{
  int
                num_part_1, num, numSD, numSD_1, numSD_N, currSD, k, l, 
                kl_phase, phase, count, *k_ptr, *l_ptr, search_num;
  UL
                *listSD, *sd, sd_k, sd_kl, new_sd, low, high;
  double
                value, *in_ptr, *out_ptr;
  ID_INT
                **table, **tab_k, *ij_ptr;
  MASK
                *mask_ptr;

   switch(op_type) {
      case VEFF: table = gr_basis->op_int.id_nondiag[0];  /* <|H(op_type)|>*/
                 break;
      case ANG:  table = gr_basis->op_ang.id_nondiag[0]; /* <|J**2(op_type)|> */ 
                 break;
      default:   printf("\n\nError in function pn_nondiag_pp_elem():");
                 printf("\nIncorrect two-particle matrix element table");
                 printf(" op_type = %d -- should be VEFF or ANG\n", op_type);
                 exit(1);
   } 
  num_part_1 = gr_basis->part[0] - 1;             /* initialization */
  numSD      = gr->numSD[0];
  numSD_1    = numSD - 1;
  numSD_N    = gr->numSD[1]; 
  num        = (gr_basis->m_orb[0] << 1) - 3;
  listSD     = gr->SD[0];
  sd         = listSD;
  for(currSD = 0; currSD < numSD_1; sd++, currSD++) {
    k_ptr  = gr_basis->occ[0] + currSD * gr_basis->part[0];
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
	  
	  if((count = gr_basis->part_lim[0].num) > 0) {    /* check masking */
	    mask_ptr = gr_basis->part_lim[0].table;
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
            else if(new_sd > listSD[search_num])   low  = search_num + 1;
            else                                   break;
	  }
            	  /* add nondiag pp contribution */

	  value   = phase * ij_ptr->val;
          in_ptr  = in_amp + currSD * numSD_N;
          out_ptr = out_amp + search_num * numSD_N;
          count   =  numSD_N;	  
	  do {
	    *(out_ptr++) += *(in_ptr++) * value;
	  } while(--count);

	          /* ad the hermetian conjugate pp contribution */

	  in_ptr  = in_amp + search_num * numSD_N;
	  out_ptr = out_amp + currSD * numSD_N;
	  count   = numSD_N;
          do {
	    *(out_ptr++) += *(in_ptr++) * value;
	  } while(--count);
	} while((++ij_ptr)->one); /* all two-body matrix elements for given (kl) */
      } /* end l-loop */
    } /* end k-loop */
  }  /* end loop through all |SD(part_type)> */
} /* End: function pn_nondiag_pp_SD_elem() */

    /*
    ** The function
    **           pn_nondiag_nn_SD_elem()
    ** calculates all nondiagonal neutron-neutronmatrix elements 
    **      <SD(N)g.n',SD(Z)g.p'|OP(nn)| SD(Z)g,p,SD(N)g,n> 
    ** within one proton/neutron group only. Pointers to all two-particle
    ** matrix elements are stored in structure gr_basis->op_int.id_table[1]
    ** by the formula
    **        op(k,l) = ((num - k) * k)/2 -1 +l
    ** where num = 2 * gr.num_orb[1] - 3.
    */

static void pn_nondiag_nn_SD_elem(int op_type, GR_BASIS *gr_basis, GROUP *gr,
                                             double *in_amp, double *out_amp)
{
  int
                num_part_1, num, numSD, numSD_1, numSD_Z, currSD, k, l, 
                kl_phase, phase, count, *k_ptr, *l_ptr, search_num;
  UL
                *listSD, *sd, sd_k, sd_kl, new_sd, low, high;
  double 
                value, *in_ptr, *out_ptr;
  ID_INT
                **table, **tab_k, *ij_ptr;
  MASK
                *mask_ptr;

   switch(op_type) {
      case VEFF: table = gr_basis->op_int.id_nondiag[1];  /* <|H(op_type)|>*/
                 break;
      case ANG:  table = gr_basis->op_ang.id_nondiag[1]; /* <|J**2(op_type)|> */ 
                 break;
      default:   printf("\n\nError in function pn_nondiag_nn_elem():");
                 printf("\nIncorrect two-particle matrix element table");
                 printf(" op_type = %d -- should be VEFF or ANG\n", op_type);
                 exit(1);
   } 
  num_part_1 = gr_basis->part[1] - 1;            /* initialization */
  numSD      = gr->numSD[1];
  numSD_1    = numSD - 1;
  numSD_Z    = gr->numSD[0];
  num        = (gr_basis->m_orb[1] << 1) - 3;
  listSD     = gr->SD[1];
  sd         = listSD;
  for(currSD = 0; currSD < numSD_1; sd++, currSD++) {
    k_ptr  = gr_basis->occ[1] + currSD * gr_basis->part[1];
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
	  
	  if((count = gr_basis->part_lim[1].num) > 0) {    /* check masking */
	    mask_ptr = gr_basis->part_lim[1].table;
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
            else if(new_sd > listSD[search_num])   low  = search_num + 1;
            else                                   break;
	  }
	          /* add nondiag nn contribution */

	  value     = phase * ij_ptr->val;
          in_ptr  = in_amp + currSD;
          out_ptr = out_amp + search_num;
          count   =  numSD_Z;	  
	  do {
	    *out_ptr += *in_ptr * value;
	    in_ptr   += numSD;
	    out_ptr  += numSD;
	  } while(--count);

	          /* ad the hermetian conjugate nn contribution */

	  in_ptr  = in_amp + search_num;
	  out_ptr = out_amp + currSD;
	  count   = numSD_Z;	 
          do {
	    *out_ptr += *in_ptr * value;
	    in_ptr   += numSD;
	    out_ptr  += numSD;
	  } while(--count);
	} while((++ij_ptr)->one); /* all two-body matrix elements for given (kl) */
      } /* end l-loop */
    } /* end k-loop */
  }  /* end loop through all |SD(part_type)> */
} /* End: function pn_nondiag_nn_SD_elem() */

     /*
     ** The function
     **           pn_nondiag_pn_SD_elem_A()
     ** calculates all nondiagonal pn-matrix elements 
     **  <SD(Z)g,p';SD(N)g,n'|OP(pn)|SD(Z)g,p,SD(N)g,n>
     ** within a single group The two-particle matrix elements are
     ** specified through the parameter op_type.
     */   

static void pn_nondiag_pn_SD_elem_A(int op_type, GR_BASIS *gr_basis, GROUP *gr,
                                                   double *in_amp, double *out_amp)
{
  int       
             p_num, n_num, num1, num2, numSD_Z, numSD_N,  n_lim, k, m,
             *k_ptr, l, *l_ptr;
  UL
             *listSD_Z, *listSD_N, *sdZ, *sdN, sdZ_k, sdN_l, high, low;
  double
             value;
   PN_INT
             **table, **tab_p, *ij_ptr;
   MASK
             *mask_ptr;
   RESULT
             result[2];

  switch(op_type) {
      case VEFF: table = gr_basis->op_int.pn_nondiag[0];            /* <|H(pn)|> */
                 break;
      case ANG:  table = gr_basis->op_ang.pn_nondiag[0];            /* <|J**2(pn)|>*/ 
                 break;
      default:   printf("\n\nError in function pn_nondiag_pn_SD_elem_A():");
                 printf("\nIncorrect two-particle matrix element table");
                 printf(" type = %d -- should be 1 - 2\n", op_type);
                 exit(1);
  } 
  numSD_Z  = gr->numSD[0];
  numSD_N  = gr->numSD[1];
  listSD_Z = gr->SD[0];
  listSD_N = gr->SD[1];

             /* run through all init |SD(Z), SD(N)> */

  sdZ = gr->SD[0];

  for(p_num = 0; p_num < numSD_Z; p_num++, sdZ++) {
    sdN      = gr->SD[1];
    n_lim    = (p_num < (numSD_Z -1)) ? numSD_N : numSD_N -1;
    for(n_num = 0; n_num < n_lim; n_num++, sdN++) { 
      num1  = p_num * numSD_N + n_num;
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
	      high = numSD_Z; 
	      while (1) {
		result[0].sd_num = (low + high) >> 1;
		if(result[0].sd < listSD_Z[result[0].sd_num])
		                                        high = result[0].sd_num-1;
		else if(result[0].sd > listSD_Z[result[0].sd_num])
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
		if(result[1].sd < listSD_N[result[1].sd_num]) 
                                                     high = result[1].sd_num - 1;
		else if(result[1].sd > listSD_N[result[1].sd_num]) 
                                                     low  = result[1].sd_num + 1;
		else                               break;
	      }
            }
	    value = result[0].phase * result[1].phase * ij_ptr->val;
	    num2  = result[0].sd_num * numSD_N  + result[1].sd_num;
            out_amp[num2] += in_amp[num1] * value;

  	         /* hermetian conjugate contribution */

	    out_amp[num1] += in_amp[num2] * value; 

	  } while((++ij_ptr)->one[0]); /* loop through all matr.elem for fixed kl-pairs */
	} while(--l);   /* end l-loop for neutron particles */
      } while(--k);      /* end k-loop for proton particles */
    } /* end n_num */
  } /* end p_num */
}  /* End: function pn_nondiag_pn_SD_elem_A() */

     /*
     ** The function
     **           pn_nondiag_pn_SD_elem_B()
     ** calculates all nondiagonal pn-matrix elements 
     **  <SD(Z)g',p';SD(N)g',n'|OP(pp) + OP(nn) + OP(pn)|SD(Z)g,p,SD(N)g,n>
     ** between two groups g -> g'. The two-particle matrix elements are
     ** specified through the parameter op_type.
     */   

static void pn_nondiag_pn_SD_elem_B(int op_type, GR_BASIS *gr_basis, GROUP *init_gr,
                                GROUP *final_gr, double *in_amp1, double *in_amp2,
				                 double *out_amp1, double *out_amp2)
{
  int       
             p_num, n_num, num1, num2, init_numSD_Z, init_numSD_N, final_numSD_Z,
             final_numSD_N, k, m, *k_ptr, l, *l_ptr, step_m_p;
  UL
             *listSD_Z, *listSD_N, *sdZ, *sdN, sdZ_k, sdN_l, high, low;
  double 
             value;
   PN_INT
             **table, **tab_p, *ij_ptr;
   MASK
             *mask_ptr;
   RESULT
             result[2];

      /* calculate changes in M_Z and PAR_Z from init to final group */

   if(gr_basis->parZ) {
     step_m_p =  (init_gr->m[0] - final_gr->m[0]) >> 1;
   }
   else {
     step_m_p =  (init_gr->m[0] - final_gr->m[0]) 
            + (((init_gr->par[0] * final_gr->par[0]) == +1) ? 0 : 1);
   }
  switch(op_type) {
      case VEFF: table = gr_basis->op_int.pn_nondiag[step_m_p];       /* <|H(pn)|> */
                 break;
      case ANG:  table = gr_basis->op_ang.pn_nondiag[step_m_p];       /* <|J**2(pn)|>*/ 
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
  listSD_Z      = final_gr->SD[0];
  listSD_N      = final_gr->SD[1];

             /* run through all init |SD(Z), SD(N)> */

  sdZ = init_gr->SD[0];

  for(p_num = 0; p_num < init_numSD_Z; p_num++, sdZ++) {
    sdN      = init_gr->SD[1];
    for(n_num = 0; n_num < init_numSD_N; n_num++, sdN++) { 
      num1  = p_num * init_numSD_N + n_num;
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
		if(result[0].sd < listSD_Z[result[0].sd_num])
		                                        high = result[0].sd_num-1;
		else if(result[0].sd > listSD_Z[result[0].sd_num])
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
		if(result[1].sd < listSD_N[result[1].sd_num]) 
                                                     high = result[1].sd_num - 1;
		else if(result[1].sd > listSD_N[result[1].sd_num]) 
                                                     low  = result[1].sd_num + 1;
		else                               break;
	      }
            }
	    value = result[0].phase * result[1].phase * ij_ptr->val;
	    num2  = result[0].sd_num * final_numSD_N  + result[1].sd_num;
            out_amp2[num2] += in_amp1[num1] * value;

  	         /* hermetian conjugate contribution */

	    out_amp1[num1] += in_amp2[num2] * value; 
 
        } while((++ij_ptr)->one[0]); /* loop through all matr.elem for fixed kl-pairs */
	} while(--l);   /* end l-loop for neutron particles */
      } while(--k);      /* end k-loop for proton particles */
    } /* end n_num */
  } /* end p_num */
}  /* End: function pn_nondiag_pn_SD_elem_B() */
