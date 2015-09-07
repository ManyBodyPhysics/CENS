/*******************  The module lanc-process.c  ******************/

#include "shell.h"

     /*
     ** The entrance function
     **         pn_eigenvalue_Lanczos_process()
     ** performs a Lanczos iteration to a convergence criterium is reached
     ** or a maximum iteration limit. A temporary file is needed to store the 
     ** basis Lanczos vectors. The corresponding energy matrix is diagonalized.
     ** The function returns the set of eigenvalues wanted together with its
     ** its eigenvectors in SD basis stored in file file_name->eigen_vec.
     ** Necessary data from the calculation is stored in file file_name->h_final
     ** if a restarting Lanczos process is wanted.
     */

            /**** local data definition ****/

  typedef   struct {
     double
               diag_elem,
            nondiag_elem;
  } RESTART_MATR;   

	    /**** local function declarations ****/

static int num_J_values(GR_BASIS *gr_basis);
     /*
     ** calculates and return the number different
     ** J-values for the eigenstates in the model
     */

static int pn_max_m_value(int type, int parity, GR_BASIS const *gr_basis);
     /*
     ** calculates for a given parity maximum M-values of a set
     ** |SD(type)M P> for a given number of identical particles
     ** in a single-particle orbits specified in mbas[]. 
     ** The function returns maximum M-value.
     */

static void initialization_process(GR_BASIS *gr_basis, ITERATION *iteration,
                                                          double *start_vec);
     /*
     ** returns a start vector for a Lanczos iteration process
     */

static void random_start_vector(int numSD, double *start_vec);
     /*
     ** calculates and return a normalized start vector where
     **  amplitudes are selected through a random procedure
     */

static void append_Lanczos_vector_to_file(char const *file_name, int const n, 
                                        int const dim, double const *vec_ptr);
     /*
     ** opens the file "file_name" and store a vector with 
     **          int n             -   as local vector number on current file and
     **          int dim           -   as number of components
     **          double vec_ptr[]  -   dim amplitudes 
     ** If n > 0 the the wave functions are appended to the previous ones.
     ** If n = 0 a new file is created if "file_name" does  not exists.
     ** If n = 0 and "file_name" does exist the vector overwrites the previous ones.
     ** Note that the first vector on file is numbered n = ZERO.
     */

static int orthogonalization_process(char *file_name, double diag_elem, 
                       double *nondiag_elem,int n, int vec_dim, LANC vec);
     /*
     ** orthogonalizes a new vector to all previous Lanczos vectors, calculates 
     ** a new off diagonal matrix element and finally normalize the new vector.
     ** if(norm > matrix_limit) the function return TRUE. Otherwise FALSE
     ** indicating that no more Lanczos vectors are available 
     */

static int eigen_value_convergence(int num_eigenstates, double *delta_eigen, 
                                                 int h_dim, double *h_matrix);
     /*
     ** calculates and stores the changes in eigenvalues of the model.states lowest
     ** eigenstates. The changes for the last four iterations are saved. If any 
     ** changes are larger than accuracy the function returns FALSE,
     ** otherwise TRUE.
     */

static void save_final_data(ITERATION *iteration, int dim_h, double *h_matrix, int dim_SD, LANC vec);
    /*
    ** transforms the eigenvectors into a |SD> basis
    ** and save them in file_name->eigen_vec.save_eigenvectors.
    ** Next, calculate the energy matrix suited for a restarting 
    ** Lanczos process and save in file iteration->h_final.
    ** Finally, store in the same file the last Lanczos vector 
    ** produced in the calculateion.
    */

static void lanczo_to_SD_basis(char *file_name, int lanc_dim, double *lanc_amp,
                                                          int sd_dim, LANC vec);
     /*
     ** transform an eigenvector from Lanczos vector basis 
     ** into SD basis and store the result in vec.two[].   
     */

static void save_final_energy_matrix(char *file_name,int num_of_elements, RESTART_MATR *matrix);
    /*
    ** writes to file the energy matrx suitable
    ** for a restarting Lanczos process 
    */

static void angular_projection_process(GR_BASIS *gr_basis, GROUP *group,
                                              ITERATION *iteration ,LANC vec);
     /*
     ** takes a new normalized Lanczos vector in |vec.two> generated through 
     ** the Lanczos process H*|> ---> |vec.two>, and 
     ** calculates ang = <vec.two|J**2|vec.two>.
     **   1. If(fabs(ang - ang_ref) < ANG_LIMIT) 
     **      the function returns without any modification of |vec.two>. 
     **   2. If(fabs(ang - ang_ref) > ANG_LIMIT) the function starts a Lanczos
     **      procedure: J**2|vec.two> ---> |final> to find the first eigenstate 
     **      with <eigen|J**2|eigen> = ang_ref. This vector is saved and returned
     **      in |vec.two> as the new corrected lanczo vector.
     **   3. If the last process is not successfull the function terminates with
     **      an error message.
     */

static double ang_mom_convergence(char *file_name, double ang_ref, int dim_j, 
                                       double *j_matrix, int dim_SD,LANC vec);
     /*
     ** diagonalizes the current j_matrix[] and search the eigenvalues
     ** to find the lowest with correct value = ang_ref.
     ** If the search is successfull the corresponding eigenvector is 
     ** transformed into |SD> basis and stored and returned in |vec.two>. 
     ** Otherwise the function terminates the process.
     */
             /**** End: function declarations ****/

               /**** The function definitions  ****/ 
     /*
     ** The entrance function 
     **          pn_eigenvalue_Lanczos_process()
     ** performs a Lanczos iteration to a convergence criterium is reached
     ** or a maximum iteration limit. A temporary file is needed to store the 
     ** basis Lanczos vectors. The corresponding energy matrix is diagonalized.
     ** The function returns the set of eigenvalues wanted together with its
     ** its eigenvectors in SD basis stored in file file_name->eigen_vec.
     ** Necessary data from the calculation is stored in file file_name->h_final
     ** if a restarting Lanczos process is wanted.
     */

void pn_eigenvalue_Lanczos_process(GR_BASIS *gr_basis, GROUP *group,ITERATION *iteration)
{
  static char
               *func      = {"pn_eigenvalue_Lanczos_process(): "};
  int          
               stored_lanc_vec, run_code, dim_h, iterate, k;
  double 
               *temp, diag_elem, nondiag_elem, *h_matrix, *delta_eigen;
  LANC
               ref_vec, curr_vec;
  TID
              ex_time;

  strcpy(iteration->lanc_store,"lanc-store.dat");/* temp file for Lanczos vectors */

         /* memory for two lanczos vectors */

  ref_vec.one = MALLOC(gr_basis->tot_numSD, double, func, "ref_vec.one[]");
  ref_vec.two = MALLOC(gr_basis->tot_numSD, double, func, "ref_vec.two[]");

         /* memory for the Lanczos energy matrix */

  dim_h    = iteration->max_iterations + 1;
  h_matrix = CALLOC((dim_h * (dim_h + 1))/2, 
                                                double, func, "h_matrix[]");
 
         /* memory to save dynamical changes in eigenvalues */

  delta_eigen = MALLOC(5 * iteration->states, double, func, "delta_eigen[]");

  if(strstr(iteration->type_calc,"fixed-J")) { /* basic data for fixed_J calculation */
    iteration->num_totJ = num_J_values(gr_basis); 
  }
           /* normalized start vector for Lanczos iteration process */

  initialization_process(gr_basis, iteration, ref_vec.one);

         /* an initial lanczo vector has been found  - append to file */

  stored_lanc_vec = 0;
  append_Lanczos_vector_to_file(iteration->lanc_store, stored_lanc_vec++, 
                                           gr_basis->tot_numSD, ref_vec.one);
  printf("Stored %d energy lanczo vectors\n",stored_lanc_vec);    /* running test */
  run_time(2,&ex_time);
  printf("run time: %d hour %d min %d sec\n",ex_time.hour, ex_time.min, ex_time.sec);

         /* 
	 **  run_code  = 0 - reach (model->max_iterations) iterations
         **            = 1 - Lanczos process has converged
	 **            = 2 - no more Lanczos vectors
	 */ 

  run_code     = 0;                                /* initialization */
  dim_h        = 0;
  curr_vec.one = ref_vec.one;                 /* dynamic addresses for lanczo vectors */
  curr_vec.two = ref_vec.two;

  for(iterate = 0; iterate < iteration->max_iterations; iterate++) {

           /* initialization of a new Lanczos vector */

    for(k = 0; k < gr_basis->tot_numSD; k++) curr_vec.two[k] = D_ZERO; 

         /*
         ** performs the operation
         **     H|init_vec> = |final_vec>
         ** where the energy operator H specified through
         ** its matrix elements in m-scheme.
         */

    matrix_vector_multiplication(gr_basis, group, VEFF, curr_vec);

    diag_elem = scalar_product_of_vectors(gr_basis->tot_numSD, curr_vec.one,curr_vec.two);

             /* diagonal element in energy matrix */

    h_matrix[((stored_lanc_vec - 1) * (stored_lanc_vec + 2))/2] = diag_elem;

    dim_h++;                      /* current dimesion of h_matrix[] */

         /*
         ** orthogonalize |cur_vec.two>  against all previous
         ** Lanczos vectors stored on file. Then normalize
         ** and return a new Lanzcos vector in vec.two[]
         */

    if(!orthogonalization_process(iteration->lanc_store, diag_elem, &nondiag_elem, 
                                 stored_lanc_vec - 1, gr_basis->tot_numSD,curr_vec)) {
      run_code = 2;                       /* no more lanczo vector available */
      break;
    }
                /* nondiag element in energy matrix */ 

    h_matrix[(stored_lanc_vec * (stored_lanc_vec + 3))/2 - 1] = nondiag_elem;
 
           /* New Lanczos vector has been found.
	   ** If single-J calculation is required, check angular momentum.
           **   1. if(fabs(<final|J**2|final> - ang_ref) < ANG_LIMIT) 
           **                        continue:
           **   2. if(fabs(<final|J**2|final> - ang_ref) > ANG_LIMIT)
           **      modify |final> to fullfill the angular momentum
	   **      requirement. Then continue.
           **      If not possible the whole procedure stops.
           */

    if(strstr(iteration->type_calc,"fixed-J")) {
      angular_projection_process(gr_basis, group, iteration,curr_vec);
    }
         /* new lanczo vector has been found  - append to file */

    append_Lanczos_vector_to_file(iteration->lanc_store, stored_lanc_vec++,
                                             gr_basis->tot_numSD,  curr_vec.two);

    printf("Stored %d energy lanczo vectors\n",stored_lanc_vec);    /* running test */
    run_time(2,&ex_time);
    printf("run time: %d hour %d min %d sec\n",ex_time.hour, ex_time.min, 
                                                               ex_time.sec);
         /*
         ** The energy spectrum is tested against a given criterium. If 
         ** function returns TRUE energy convergence has been obtained.
         */

    if(eigen_value_convergence(iteration->states, delta_eigen, dim_h, h_matrix)) {
      run_code = 1;
      break;
    }
    temp         = curr_vec.one;   /* interchange the pointers */
    curr_vec.one = curr_vec.two;
    curr_vec.two = temp;
  } /* end lanczo iteration loop */

  free(delta_eigen);                           /* release local memory */

  iteration->run_code       = run_code; /* save run_code for output process */
  iteration->num_iterations = iterate;
  iteration->states         = MIN(iteration->states, dim_h);

           /* calculate and save final data */

  save_final_data(iteration, dim_h, h_matrix, gr_basis->tot_numSD, ref_vec);

  remove(iteration->lanc_store);                   /* remove temp file */
  strcpy(iteration->lanc_store,"");

  free(h_matrix);                                 /* remove local mempry */
  free(ref_vec.two);
  free(ref_vec.one);
  
} /* End: function  eigenvalue_Lanczos_process()*/

      /*
      ** The function
      **     num_J_values()
      ** calculates and return the number different
      ** J-values for the eigenstates in the model
      */

static int num_J_values(GR_BASIS *gr_basis)
{
  int     max_MZ1, max_MN1, max_MZ2, max_MN2, max_M;
 
  if(gr_basis->parZ == 0) {
    max_MZ1 = pn_max_m_value(PROTON, +1, gr_basis);
    max_MN1 = pn_max_m_value(NEUTRON, gr_basis->P, gr_basis);
    max_MZ2 = pn_max_m_value(PROTON, -1, gr_basis);
    max_MZ1 = pn_max_m_value(PROTON, -gr_basis->P, gr_basis);

    max_M = MAX(max_MZ1 + max_MN1, max_MZ2 + max_MN2);
  }
  else {
    max_M =  pn_max_m_value(PROTON, gr_basis->parZ, gr_basis);
          + pn_max_m_value(NEUTRON, gr_basis->parZ * gr_basis->P, gr_basis);
  }
  return ((max_M - gr_basis->MJ)/2 + 1);

} /* End: function num_J_values() */

     /*
     ** The function 
     **           pn_max_m_value()
     ** calculates for a given parity maximum M-values of a set
     ** |SD(type)M P> for a given number of identical particles
     ** in a single-particle orbits specified in mbas[]. 
     ** The function returns maximum M-value.
     */

static int pn_max_m_value(int type, int parity, GR_BASIS  const *gr_basis)
{ 
   char          *func = {"pn_max_m_value(): "};
   register int  loop, m, par, num_part, num_orb;
   int           *matr, *ptr;
   MBAS          *mbas;


   num_part = gr_basis->part[type];                         /* initialization */ 
   num_orb  = gr_basis->m_orb[type];
   mbas     = gr_basis->mbas[type];

   matr = MALLOC(num_part + 1, int, func, "matr[]");     /* temporary memory */

   for(loop = 0; loop < num_part; loop++) matr[loop] = loop; /* lowest config. */
   matr[num_part] = num_orb;

   do  {                                                  /* particle loop */
      for(m = 0, par = +1, ptr = matr, loop = 0; loop < num_part; loop++) {
         m   += mbas[*ptr].m; 
         par *= mbas[*(ptr++)].par; 
      }
      if(par == parity) break;

              /* new particle configuration */

      for(loop = 0; (loop < num_part && ((matr[loop]+1) >= matr[loop+1])); loop++);
      matr[loop]++;
      for(loop--; loop >= 0; loop--)  matr[loop] = loop;

   } while(matr[num_part] == num_orb);  /* end of particle loop */

   if(matr[num_part] != num_orb) {
      printf("\n\nError in function pn_max_m_value():");
      printf("\nNo %s %d particle configuration with parity %c is possible\n",
             ((type == 0) ? "PROTON" : "NEUTRON"),num_part,((parity == +1) ? '+' : '-')),
      exit(1);
   }
   free(matr);  /* release temporary memory */

   return m;                 /* return max m-value */

} /* End: function pn_max_m_value() */

      /*
      ** The function 
      **   initialization_process();
      ** returns a start vector for a Lanczos iteration process
      */

static void initialization_process(GR_BASIS *gr_basis, ITERATION *iteration,
                                                           double *start_vec)
{
  switch(strlen(iteration->type_calc)) {
      case 12 :                                   /* random-start */
               random_start_vector(gr_basis->tot_numSD, start_vec);
               break;
      case 15 :                                   /* random-continue */
               printf("\n\nError in function initialization_process():");
               printf("\ntype_calc = %s",iteration->type_calc);
               printf("\n Not implemented!!!!!!!!!\n");
               exit(1);
               break;
      case 13 :                                   /* fixed-J-start */
               printf("\n\nError in function initialization_process():");
               printf("\ntype_calc = %s\n",iteration->type_calc);
               printf("\n Not implemented!!!!!!!!!\n");
               exit(1);
               break;
      case 16 :                                   /* fixed-J-continue */
               printf("\n\nError in function initialization_process():");
               printf("\ntype_calc = %s\n",iteration->type_calc);
               printf("\n Not implemented!!!!!!!!!\n");
               exit(1);
               break;
      default :
               printf("\n\nError in function initialization_process():");
               printf("\nWrong type_calc = %s\n",iteration->type_calc);
               exit(1);
               break;
  } /* end switch() */

} /* End: function initialization_process() */

      /*
      ** The function                           
      **        random_start_vector()                
      ** calculates and return a normalized start vector where
      **  amplitudes are selected through a random procedure
      */

static void random_start_vector(int numSD, double *start_vec)
{
  int
            phase, loop;
  double
           *ptr, value, norm;

  ptr = start_vec;                                                /* initialization */

/**********************************
  srand((unsigned long) time(NULL)); 
  phase = (rand() & TRUE) ? +1 : -1;
**************************************/

  phase = +1;

  norm = D_ZERO;
  for(loop = 0; loop < numSD; loop++) {
    *(ptr++) = (value = 1.0+((double)loop)/numSD,value*phase);
    norm +=  value * value;

/****************************
    phase = (rand() & TRUE) ? +1 : -1;
******************************/

    phase = -phase;

  }
      /* Normalization of the  lanczo vector */

  if(fabs(norm) < ZERO_LIMIT)  {
    printf("\n\nError in function  random_start_vector():");
    printf("\nNorm of the initial random lanczos vector is below");
    printf(" accepted limit: norm = %E",norm);
    exit(1);
  }
  norm = 1.0 / sqrt(norm);
  ptr  = start_vec;
  for(loop = 0; loop < numSD; loop++)  {
    *(ptr++) *= norm;
  }
} /* End: function random_start_vector() */

     /*
     ** The function                                           
     **            append_Lanczos_vector_to_file( )                    
     ** opens the file "file_name" and store a vector with 
     **          int n             -   as local vector number on current file and
     **          int dim           -   as number of components
     **          double vec_ptr[]  -   dim amplitudes 
     ** If n > 0 the the wave functions are appended to the previous ones.
     ** If n = 0 a new file is created if "file_name" does  not exists.
     ** If n = 0 and "file_name" does exist the vector overwrites the previous ones.
     ** Note that the first vector on file is numbered n = ZERO.
     */

static void append_Lanczos_vector_to_file(char const *file_name, int const n, 
                                         int const dim,  double const *vec_ptr)
{
/****
  int 
            k;
****/
  FILE
           *file_ptr;

/*******************   test output ************
  printf("\n\nStorage of Lanczos vector no %d in file %s\n", n, file_name);
  for(k = 0; k < dim; k++) {
    printf("  %10.6f",vec_ptr[k]);
    if((((k + 1)/6) * 6) == (k + 1)) printf("\n");
  }
  printf("\n");
******************    end  ******************/

  if((file_ptr = fopen(file_name,(n == 0) ? "wb": "ab+")) == NULL) {
    printf("\nError in function append_Lanczos_vector_to_file();");
    printf("\nWrong file = %s to store a vector no %d", file_name,n);
    exit(1);
  }
  if(fwrite((const void *)&n,(size_t) sizeof(int), 1, file_ptr) != 1) {/* vector number */
    printf("\nError in function append_Lanczos_vector_to_file()");
    printf("\nIn writing vector number =  %d", n);
    printf("\nto file %s\n",file_name);
    exit(1);
  }   /* end of if-test  */

  if(fwrite((const void *)&dim,(size_t) sizeof(int),1,file_ptr)!=1) {/* write dimension */
    printf("\nError in function append_Lanczos_vector_to_file()");
    printf("\nIn writing numSD =  %d", dim);
    printf("\nto file %s\n",file_name);
    exit(1);
  }   /* end of if-test  */

  if(fwrite((const void *)vec_ptr,(size_t) sizeof(double),  /* write vector */
	           (size_t) dim, file_ptr) != (size_t) dim)  {
    printf("\nError in function append_Lanczos_vector_to_file()");
    printf("\nIn writing a vector to file %s\n",file_name);
    exit(1);
  } /* end of if-test  */

  fclose(file_ptr);

} /* End: function append_Lanczos_vector_to_file() */

    /*
    ** The function
    **                  orthogonalization_process()                        
    ** orthogonalizes a new vector to all previous Lanczos vectors, calculates 
    ** a new off diagonal matrix element and finally normalize the new vector.
    ** if(norm > matrix_limit) the function return TRUE. Otherwise FALSE
    ** indicating that no more Lanczos vectors are available 
    */

static int orthogonalization_process(char *file_name, double diag_elem, 
                       double *nondiag_elem, int n, int vec_dim, LANC vec)
{
  int
               k;
  double
               factor;
  FILE
              *file_ptr;

       /*
       ** orthogonalize to last stored Lanczos vector
       **   (|vec.two[]> - diag_elem * |vec.one[]>)  ---> |vec.two[]>  
       */

  add_a_vector(vec_dim, -diag_elem, vec.one, vec.two);

      /* Orthogonalize to all previous Lanczo vectors.*/

  if( (file_ptr = fopen(file_name,"rb+")) == NULL) {
    printf("\n\nError in function  orthogonalization_process():"); 
    printf("\nWrong file = %s", file_name);
    printf("\nto store lanczo vectors.\n");
    exit(1);
  }
  for(k = 0; k < n ; k++)  {

           /* Read |lanc_vec(k)> into |vec.one[]> */

    read_Lanczos_vector_from_open_file(file_name, file_ptr, k, vec_dim, vec.one);

           /* factor = <vec.one[] | vec.two[]>  */

    if(fabs(factor = scalar_product_of_vectors(vec_dim, vec.one, vec.two)) 
                                                       < ZERO_LIMIT) continue;

          /*
          ** Orthogonalize the new vector to lanczo vector k  
          **       |vec.two[]> <--- |vec.two[]> - factor * |vec.one[]>
          */

    add_a_vector(vec_dim, -factor, vec.one, vec.two);

  } /* end of loop k */

  fclose(file_ptr);  /* close the storage file for lanczo vectors */
  
          /*
	  ** Normalize |vec.two[]> and produce a new lanc_vec. 
	  ** Calculate the off-diagonal matrix element to the 
	  ** new Lanczos vector.
	  ** The new lanc_vec must  have  norm > MATRIX_LIMIT.    
          */

  if((factor =  norm_of_a_vector(vec_dim, vec.two)) < MATRIX_LIMIT)  return FALSE;

  *nondiag_elem = factor;         /* off-diagonal matrix element */
  factor = 1.0 / factor;

       /* Normalize the new lanczo vector |vec.two[]> */

  scale_a_vector(vec_dim, factor, vec.two);

  return TRUE;

} /* End: function orthogonalization_process() */ 

      /*
      ** The function                                     
      **        eigen_value_convergence()                 
      ** calculates and stores the changes in eigenvalues of the model.states lowest
      ** eigenstates. The changes for the last four iterations are saved. If any 
      ** changes are larger than accuracy the function returns FALSE,
      ** otherwise TRUE.
      */

static int eigen_value_convergence(int num_eigenstates, double *delta_eigen, int h_dim, double *h_matrix)
{
  char           
                 *func = {"eigen_value_convergence(): "};
  register int
                 loop, num;
  static   int
                 delta_num = -1;
  double
                 *list_eigenvalue;

  list_eigenvalue = MALLOC(h_dim, double, func, "list_eigenvalue[]");

           /* Only eigenvalues are calculated */

  eigenvalues(h_dim, h_matrix, list_eigenvalue);

  if(delta_num < 0) { /* save eigenvalues and return */
    for(loop = 0; loop < num_eigenstates; loop++)  { /* save eigenvalues and return */
      delta_eigen[4 * num_eigenstates + loop] =  list_eigenvalue[loop];
    }
    delta_num++;
    free(list_eigenvalue);
    return FALSE;
  }
           /* save changes in eigenvalues */

  for(loop = 0, num = delta_num * num_eigenstates; loop < num_eigenstates; loop++, num++)  {
    delta_eigen[num] = fabs(delta_eigen[4 * num_eigenstates + loop] - list_eigenvalue[loop]);
  }
           /*
           ** Increase delta_num as a pointer to a circular
           ** save-buffer containing  maximum eight sets.                  
           */

  delta_num = (3 & (++delta_num));

  for(loop = 0; loop < num_eigenstates; loop++)  { /* save eigenvalues */
    delta_eigen[4 * num_eigenstates + loop] = list_eigenvalue[loop];
  }
  if((h_dim - num_eigenstates) < 4) {
    free(list_eigenvalue);
    return FALSE;
  }
           /*
           ** If any of the saved energy differences are larger than 
           ** energy accuracy return FALSE, else return TRUE.
           */

  for(loop = 0;loop < 4 * num_eigenstates; loop++) {
    if(delta_eigen[loop] > ENERGY_LIMIT)  {
      free(list_eigenvalue);
      return FALSE;
    }  
  }
  
  free(list_eigenvalue);      /* release local memory */
  
  return TRUE;     /* The process has converged */

} /* End: function eigen_value_convergence() */

    /*
    ** The function                                           
    **              save_final_data()                       
    ** transforms the eigenvectors into a |SD> basis
    ** and save them in file_name->eigen_vec.save_eigenvectors.
    ** Next, calculate the energy matrix suited for a restarting 
    ** Lanczos process and save in file iteration->h_final.
    ** Finally, store in the same file the last Lanczos vector 
    ** produced in the calculateion.
    */

static void save_final_data(ITERATION *iteration, int dim_h, double *h_matrix, int dim_SD, LANC vec)
{
  char           
                 *func ={"save_final_data():"};
  int
                 loop, num;
  double
                 last_nondiag_elem; 
  EIGEN
                 *list_eigen;
  RESTART_MATR
                 *restart_matr;
  FILE
                 *file_ptr;

           /* local memory final energy diagonalization */

  list_eigen = MALLOC(dim_h, EIGEN, func,"list_eigen[]");
  for(loop = 0; loop < dim_h; loop++) {
    list_eigen[loop].vector = MALLOC(dim_h, double, func, "list_eigen[].vector[]");
  }
         /* final eigenvalues and eigenvectors are calculated and sorted */

  eigenvectors(dim_h, h_matrix, list_eigen);

          /* save iteration->states final eigenvalues for output */

  for(loop = 0; loop < iteration->states; loop++) {
    iteration->eigenvalues[loop] = list_eigen[loop].value;
  }
	    /* transfer iteration->states eigenvectors into SD basis and save them on file */

  for(loop = 0; loop < iteration->states; loop++) {

      /* transform eigenvector no loop to SD basis and store in vec_two[]. */

    lanczo_to_SD_basis(iteration->lanc_store, dim_h, list_eigen[loop].vector, dim_SD, vec);

       /* save the eigenvector on file */

    append_Lanczos_vector_to_file(iteration->eigen_vec, loop, dim_SD, vec.two);

  } /* end loop:  all eigenvector written to file */

        /* if(iteration->run_code != 2)
	** calculate and save on file all necessary data 
	** for a possible restart Lanczos process
        */

  if(iteration->run_code != 2) { 
    restart_matr = MALLOC(iteration->states, RESTART_MATR, func, "restart_matr[]");
 
    last_nondiag_elem = h_matrix[(dim_h * (dim_h + 3))/2 - 1];

    for(loop = 0, num = dim_h - 1; loop < iteration->states; loop++) {
      restart_matr[loop].diag_elem    = list_eigen[loop].value;
      restart_matr[loop].nondiag_elem = list_eigen[loop].vector[num] * last_nondiag_elem;
    }
         /* save the restarting energy matrix on file */

    save_final_energy_matrix(iteration->h_final, iteration->states, restart_matr);

            /* open the Lanczos vector file */ 

    if( (file_ptr = fopen(iteration->lanc_store,"rb"))==NULL) {
      printf("\n\nError in function save_final_data():");
      printf("\nWrong file = %s", iteration->lanc_store);
      printf("\nfor the Lanczos vectors.");
      exit(1);
    }
    rewind(file_ptr);

      /* read the last Lanczos  vector into vec.one[] */

    num = (2 * sizeof(int) + dim_SD * sizeof(double)) * dim_h;
 
    fseek(file_ptr,(long) num, SEEK_SET);
    read_Lanczos_vector_from_open_file(iteration->lanc_store, file_ptr, dim_h, dim_SD, vec.one);
    fclose(file_ptr); 

        /* append Lanczos vector to the restarting file */
   
    append_Lanczos_vector_to_file(iteration->h_final, iteration->states, dim_SD, vec.one);

    free(restart_matr);                           /* release local memory */
    for(loop = dim_h - 1; loop >= 0; loop--) {
      free(list_eigen[loop].vector);
    }
    free(list_eigen);

  } /* end restarting procedure */

} /* End: function save_final_data() */

      /*
      ** The function                                        
      **          lanczo_to_SD_basis()                       
      ** transform an eigenvector from Lanczos vector basis 
      ** into SD basis and store the result in vec.two[].   
      */

static void lanczo_to_SD_basis(char *file_name, int lanc_dim, double *lanc_amp,
                                                          int sd_dim, LANC vec)
{
  register int
                loop;
  FILE
                *file_ptr;

  for(loop = 0; loop < sd_dim; loop++) vec.two[loop] = D_ZERO; /* initialization */

      /* Open the lanczo storage file for reading of lanczo vectors. */

  if( (file_ptr = fopen(file_name,"rb")) == NULL) {
    printf("\n\nError in function lanczo_to_SD_basis():");
    printf("\nWrong file = %s, for Lanczos vectors.",file_name);
    exit(1);
  }
  rewind(file_ptr);
  
          /* run through all Lanczos  vectors */

  for(loop = 0; loop < lanc_dim; loop++) {

    read_Lanczos_vector_from_open_file(file_name, file_ptr, loop, sd_dim, vec.one);

              /* calculate the contribution from vector(loop) */

    add_a_vector(sd_dim, *(lanc_amp++), vec.one, vec.two);

  } /* end of loop through all lanczo_vectors */

  fclose(file_ptr);   /* close the file for lanczo vectors */

} /* End: function lanczo_to_SD_basis() */

    /*
    ** The function                               
    **     void save_final_energy_matrix()        
    ** writes to file the energy matrx suitable
    ** for a restarting Lanczos process 
    */

static void save_final_energy_matrix(char *file_name,int num_of_elements, RESTART_MATR *matrix)
{
  FILE       *file_ptr;

  if((file_ptr = fopen(file_name,"wb")) == NULL)  {
    printf("\n\nError in function save_final_energy_matrix():");
    printf("\nWrong file = %s", file_name);
    printf("\nto store final energy matrix[].\n");
    exit(1);
  }
           /* write the number of elements */

  if( fwrite((const void *)&num_of_elements, (size_t) sizeof(int), 1, file_ptr) != 1) {
    printf("\n\nError in function save_final_energy_matrix():");
    printf("\nin write num_of_elements = %d", num_of_elements);
    printf("\nto file %s\n", file_name);
    exit(1);
  }   /* end of if-test  */

    /* write the energy matrix */

  if( fwrite((const void  *)matrix, (size_t) sizeof(RESTART_MATR),
                 (size_t) num_of_elements, file_ptr) != (size_t) num_of_elements) {
    printf("\n\nError in function save_final_energy_matrix():");
    printf("\nin writing the finale energy matrix to file %s\n", file_name);
    exit(1);
  }   /* end of if-test  */

  fclose(file_ptr);

}  /* End: function save_final_energy_matrix() */

      /*
      ** The function
      **       angular_projection_process()
      ** takes a new normalized Lanczos vector in |vec.two> generated through the 
      ** Lanczos process H*|> ---> |vec.two>, and calculates 
      ** ang = <vec.two|J**2|vec.two>.
      **   1. If(fabs(ang - ang_ref) < ANG_LIMIT) 
      **      the function returns without any modification of |vec.two>. 
      **   2. If(fabs(ang - ang_ref) > ANG_LIMIT) the function starts a Lanczos
      **      procedure: J**2|vec.two> ---> |final> to find the first eigenstate 
      **      with <eigen|J**2|eigen> = ang_ref. This vector is saved and returned
      **      in |vec.two> as the new corrected lanczo vector.
      **   3. If the last process is not successfull the function terminates with
      **      an error message.
      */

static void  angular_projection_process(GR_BASIS *gr_basis, GROUP *group, 
                                               ITERATION *iteration, LANC vec) 
{
  char
            *func = {"angular_projection_process(): "},
            *temp_file = {"lanc-ang.dat"};    /* temp file for J**2| > vectors */
  int
            loop, dim_j, ang_ref, k, stored_lanc_vec;
  double
            *temp, ang, *j_matrix, diag_elem, nondiag_elem;
  LANC
            curr_vec;
  TID 
            ex_time;
  FILE
            *file_ptr;

  curr_vec.one = vec.two;     /* the new Lanczos vector are store in |ref_vec.two> */
  curr_vec.two = vec.one;         

  ang_ref = 0.25 * (double)(iteration->tot_2J * (iteration->tot_2J + 2)); 

  for(loop = 0; loop < iteration->tot_2J; loop++) {

           /* initialization of a new Lanczos vector */

    for(k = 0; k < gr_basis->tot_numSD; k++) curr_vec.two[k] = D_ZERO; 

         /*
         ** performs the operation
         **     J**2|init_vec> = |final_vec>
         ** where the angular momentum operator J**2 is  specified 
	 ** through its matrix elements in m-scheme.
         */

    matrix_vector_multiplication(gr_basis, group, ANG, curr_vec);

    diag_elem = scalar_product_of_vectors(gr_basis->tot_numSD, curr_vec.one, 
                                                               curr_vec.two);
    if(loop == 0) {
      if(fabs(diag_elem - ang_ref) < ANG_LIMIT) return;  /* Lanczos vec ok */

           /* 
	   ** lanczo procedure on |init> to find a lanczo
	   ** vector with correct angular momentum 
	   **     <mod-init|J**2|mod-init> = ang_ref
           */

      j_matrix = CALLOC(((iteration->num_totJ + 1)*(iteration->num_totJ + 2))/2, 
                                                       double,func, "j_matrix[]");

      dim_j = 0;                               /* initialization */
      stored_lanc_vec = 0;
      append_Lanczos_vector_to_file(temp_file, stored_lanc_vec++, 
                                           gr_basis->tot_numSD, curr_vec.one);
    } /* end initial loop */

             /* diagonal element in j_matrix */

    j_matrix[((stored_lanc_vec - 1) * (stored_lanc_vec + 2))/2] = diag_elem;

    dim_j++;                      /* current dimesion of h_matrix[] */

         /*
         ** orthogonalize |cur_vec.two>  against all previous
         ** Lanczos vectors stored on file. Then normalize
         ** and return a new Lanzcos vector in vec.two[]
         */

    if(!orthogonalization_process("temp_ang.dat", diag_elem, &nondiag_elem, 
                                 stored_lanc_vec - 1, gr_basis->tot_numSD,curr_vec)) {
      break;
    }
                /* nondiag element in energy matrix */ 

    j_matrix[(stored_lanc_vec * (stored_lanc_vec + 3))/2 - 1] = nondiag_elem;

         /* new lanczo vector has been found  - append to file */

    append_Lanczos_vector_to_file("temp_ang.dat", stored_lanc_vec++,
                                             gr_basis->tot_numSD, curr_vec.two);

    printf("Stored J**2 %d lanczo vectors\n",stored_lanc_vec);    /* running test */
    run_time(2,&ex_time);
    printf("run time: %d hour %d min %d sec\n",ex_time.hour, ex_time.min, 
                                                               ex_time.sec);
    temp         = curr_vec.one;   /* interchange the pointers */
    curr_vec.one = curr_vec.two;
    curr_vec.two = temp;

  } /* end lanczo iteration loop */

  ang = ang_mom_convergence(temp_file, ang_ref, dim_j, j_matrix, 
                                                    gr_basis->tot_numSD, curr_vec);

  remove(temp_file);           /* temp file to store lanco vectors J**2| > */

  free(j_matrix);               /* release local memory */

  if((file_ptr = fopen(iteration->out_data,"a")) == NULL) {
    printf("\n\nError in functin angular_projection_process():");
    printf("\nWrong file = %s for the output data\n", iteration->out_data);
    exit(1);
  }
  fprintf(file_ptr,"\nLanczo vector no %d: delta(<J**2>) = %10.2E",
		                                     stored_lanc_vec + 1, ang);
  fclose(file_ptr);

}  /* End: function angular_projection_process() */

     /*
     ** The function                                     
     **            ang_mom_convergence()                 
     ** diagonalizes the current j_matrix[] and search the eigenvalues
     ** to find the lowest with correct value = ang_ref.
     ** If the search is successfull the corresponding eigenvector is 
     ** transformed into |SD> basis and stored and returned in |vec.two>. 
     ** Otherwise the function terminates the process.
     */
     
static double ang_mom_convergence(char *file_name, double ang_ref, int dim_j, 
                                            double *j_matrix, int dim_SD, LANC vec)
{
  char
               *func = {"ang_mom_convergence(): "};
  int
               loop, k,l;
  double
               ang;
  EIGEN
               *list_eigen;

  list_eigen = MALLOC(dim_j, EIGEN, func, "list_eigen[]");

  for(loop = 0; loop < dim_j; loop++) {
    list_eigen[loop].vector = MALLOC(dim_j, double, func, "list_eigen[].vector[]");
  }
         /* Both eigenvalues and eigenvectors are calculated */

  eigenvectors(dim_j, j_matrix, list_eigen);

      /* search for correct angular momentum value = ang_ref */

  for(k = 0; k < dim_j; k++) if(fabs(list_eigen[k].value - ang_ref) < ANG_LIMIT) break;

  if(k == dim_j)  {
    printf("\n\nError in function ang_mom_convergence():");
    printf("\nNo vector with correct angular momentum was found.");
    printf("\nRequired value = %f\n", ang_ref); 
    exit(1);
  }
  ang = list_eigen[k].value;        /* return value */

        /* 
        ** The new lanczo vector will be almost the same as the one which 
        ** started the angular momentum  correction process.
        ** The sign of amplitude eigen[k].vector[0] must be positive.
        */

  if(list_eigen[k].vector[0] < 0.0) {
    for(l = 0; l < dim_j; l++) list_eigen[k].vector[l] *= -1;
  }
      /* transform eigenvector no k to |SD> basis and store the result in |vec.two> */

  lanczo_to_SD_basis(file_name, dim_j, list_eigen[k].vector, dim_SD, vec);

        /* release local memory */

  for(loop = dim_j - 1; loop >= 0; loop--)  free(list_eigen[loop].vector); 
  free(list_eigen);

  return ang;

} /* End: function ang_mom_convergence() */
