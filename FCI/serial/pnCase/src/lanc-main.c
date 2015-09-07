  /*
  **          Version Sept 2001
  ** Shell model program for a system of protons and neutrons using Lanczos 
  ** iteration method.
  ** Two possibilities are implemented:
  **      1. Lanczos iteration starts with a random initial vector with mixed
  **         angular momentum components.
  **      2. Lanczos iteration starts with an initial vector with fixed angular
  **         momentum J. At every intermediate step the new lanzcos vector is
  **         checked with respect to angular momentum admixture. If necessary 
  **         the new lanczo vector is projected onto a state with the same 
  **         angular momentum J.
  ** The program calculates and stores in memory all diagonal and as many as 
  ** possible nondiagonal matrix elements of the two-particle effective
  ** interaction and J**2.
  **
  ** It calculates
  **   1. Energy levels.                             
  **   2. Angular momentum,i.e. expectation value of  J**2                           -
  **      for a set of energy eigenstates with angular momentum J.           
  **   3. The particle distribution over the single-particle orbits.
  ** The module lanc-main.c contains the following functions:
  **
  **   1. int main(int argc, char *argv[])
  **   2. void id_shell_model_calc(FILE_NAME *file_name, INPUT const *data)
  **   3. void pn_shell_model_calc(FILE_NAME *file_name, INPUT const *data)
  **   4. void pn_data_structure(FILE_NAME *file_name, INPUT const *data, PN_SHELL *model)
  **   5. int check_pn_orbits(INPUT const *data)
  **   6. void pn_print_group_data(char const *file_name, GR_BASIS const *gr_basis, 
  **                                                                  GROUP const *group)
  **   7. void pn_result(PN_SHELL *model, GR_BASIS *gr_basis, GROUP *group, 
  **                                                                     ITERATION *iteration)
  **   8. void proton_j_occupation(PN_SHELL *model, GR_BASIS *gr_basis, GROUP *group,
  **                                                    double *j_occ, double *vector)
  **   9. void neutron_j_occupation(PN_SHELL *model, GR_BASIS *gr_basis, GROUP *group,
  **                                                   double *j_occ, double *vector)
  */

#include "shell.h"

            /**** local data definitions ****/

            /**** local function declarations ****/

  void id_shell_model_calc(INPUT  *data);
      /*
      ** perform  a shell model  Lanczos iteration process  
      ** for a system of identical particles
      **         Not ready 
      */

  void pn_shell_model_calc(INPUT  *data);
      /*
      ** perform  a shell model Lanczos iteration process  
      ** for a system of both protons and neutrons
      */

  void pn_data_structure(INPUT *, PN_SHELL *, GR_BASIS *, ITERATION *);
      /*
      ** sets up the PN_MODEL data structure 
      ** for a proton/neutron calculation
      */

  int check_pn_orbits(INPUT  *data);
      /*
      ** analyzes proton and neutron spherical single-particle orbits
      ** and return YES if they are identical, otherwise NO.
      */

  void pn_print_group_data(char *file_name, GR_BASIS  *gr_basis, GROUP  *group);
      /*
      ** writes to output file all calculated group data.
      */

  void pn_result(PN_SHELL *model, GR_BASIS *gr_basis, GROUP *group, ITERATION *iteration);
      /*
      ** takes the calculated result from the proton-neutron
      ** Lanczos process and write it to output file
      */

  void proton_j_occupation(PN_SHELL *model, GR_BASIS *gr_basis, GROUP *group,
			              	            double *j_occ, double *vector);
      /*
      ** calculate the j-occupation numbers  of protons for a shell model eigenstate 
      ** calculated through the Lanczo procedure and stored  in vector[]. 
      */

   void neutron_j_occupation(PN_SHELL *model, GR_BASIS *gr_basis, GROUP *group,
				                       double *j_occ, double *vector);
      /*
      ** calculate the j-occupation numbers  of neutrons for a shell model eigenstate 
      ** calculated through the lanczo procedure an stored  in Lanczo_one[]. 
      */

            /****   end local function declarations ****/

           /**** The function definitions  ****/ 
 
int main(int argc, char *argv[])
{
  INPUT
             data;

      /* Read file names and open file containing basic shell model input data.*/

  if(argc >= 2) {                              /* as program parameters */
    strcpy(data.in_data, argv[1]);       /* save in_data filename */
  }
  else  {                             /* or as separate typed in filename */
    printf("\n Type in file name for basic shell model data =");
    scanf("%s",data.in_data);
  } /* end of data-file-name input */

    /*
    ** All data for the shell model problem and the lanczos procedure are read 
    ** and stored in INPUT data, checked and written to output file.
    */
   
  input_process(&data);

      /* type of shell model calculation */

  if(data.part_type  <= 1) {

/**********************
   id_shell_model_calc(&file_name, &data);        * identical particle calculation *
**************************/

  }
  else {                              /* proton/neutron calculation */
    pn_shell_model_calc(&data);      /* proton-neutron calculation */
  }
  return 0;          /* sucsessfull termination */
} /* End: function main()  */

   /*
   ** The function
   **        id_shell_model_calc()
   ** perform  a shell model  Lanczos iteration process  
   ** for a system of identical particles
   **         Not ready 
   */

  void id_shell_model_calc(INPUT  *data)
{

} /* End: function id_shell_model_calc() */

   /*
   ** The function
   **        pn_shell_model_calc()
   ** perform  a shell model Lanczos iteration process  
   ** for a system of both protons and neutrons
   */

  void pn_shell_model_calc(INPUT  *data)
{
  char          *func = {"pn_shell_model_calc(): "};
  PN_SHELL
                model;
  GR_BASIS
                gr_basis;
  GROUP
                *group = NULL;
  ITERATION 
                iteration;
  FILE
                *file_ptr;
  TID
                ex_time;

         /*
         ** transfer input data for proton-neutron shell model calculation
         ** to struct PN_SHELL model struct GR_BASIS gr_basis and struct
	 **  ITERATION iteration.
         */

  pn_data_structure(data, &model, &gr_basis, &iteration);

  run_time(1,&ex_time);               /* start calculation time */

   	 /* calculate proton/neutron m-scheme single-particle basis */

  pn_single_particle_basis(&model, &gr_basis);

         /* Effective interaction:
	 **    1. read angular momentum coupled effective two-particle
         **       matrix elements from file 
         **    2. calculate and store all m-scheme two-particle
         **       effective interaction matrix elements
         **    3. For single-J shell model process calculate and 
         **       store all m-scheme matrix elements of J**2
         */

  proton_neutron_effective_interaction(data, &model, &gr_basis);

     /* calculated proton neutron m-scheme slater determinant basis */

  if((file_ptr = fopen(model.out_data,"a")) == NULL) {
    printf("\n\nError in function pn_shell_model_calc()");
    printf("\nWrong file = %s for output data\n", model.out_data);
    exit(1);
  } 
  if(pn_slater_determinant(&model, &gr_basis, &group) == FALSE)  {
    fprintf(file_ptr,"\n\nNumber of basis |SD(Z), SD(N)> states are ZERO!!!\n\n");
    fclose(file_ptr);
    exit(1);
  }
               /* Total number of basis Slater determinant */

  fprintf(file_ptr,"\n\nNumber of proton |SD> configurations = %d", gr_basis.tot_numSD_Z); 
  fprintf(file_ptr,"\nNumber of neutron |SD> configurations = %d", gr_basis.tot_numSD_N); 
  fprintf(file_ptr,"\nTotal number of proton-neutron |SD> configurations = %d",
                                                                     gr_basis.tot_numSD); 
  fclose(file_ptr);

/**************    test output *************
  pn_print_group_data(model.out_data, &gr_basis, group);
***************     end test output  ***********/

         /*
         ** Calculate and store diagonal and nondiagonal slater
	 ** determinant matrix elements of V_eff.
         */

   pn_store_SD_matr_elem(VEFF, &gr_basis, group);

         /*
         ** Calculate and store diagonal and nondiagonal slater
	 ** determinant matrix elements of V_ang.
         */

   pn_store_SD_matr_elem(ANG, &gr_basis, group);

         /* 
	 ** Lanczos iteration process to find eigenvalues 
         ** and eigenvectors for a proton-neutron system
         */

   run_time(2,&ex_time);                               /* check initialization time */
   printf("\n\nThe initialization process:\n");
   printf("run time: %d hour %d min %d sec\n\n",
               ex_time.hour, ex_time.min, ex_time.sec);

            /* memory to store the final eigenvalues */

   iteration.eigenvalues = MALLOC(iteration.states,double,func, "eigenvalues[]");

   pn_eigenvalue_Lanczos_process(&gr_basis, group, &iteration);  

   pn_result(&model, &gr_basis, group, &iteration);

} /* End: function pn_shell_model_calc() */

       /*
       ** The function 
       **      pn_data_structure()
       ** sets up the data structures PN_SHELL model, GR_BASIS gr_basis
       ** and ITERATION iteration  for a proton/neutron calculation
       */

void pn_data_structure(INPUT *data, PN_SHELL *model, GR_BASIS *gr_basis, ITERATION *iteration)
{
  char        string[ONE_LINE];
 
            /* copy data to model */

  strcpy(model->title, data->title);   /* file names for output files */
  strcpy(model->type_calc, data->type_calc);
  strcpy(model->out_data, data->out_data);

  model->part_type      = data->part_type;
  model->P              = data->P;
  model->MJ             = data->MJ;
  model->part[0]        = data->part[0];
  model->part[1]        = data->part[1];
  model->num_j[0]       = data->num_j[0];
  model->num_j[1]       = data->num_j[1];

  model->jbas[0] = data->jbas[0]; /* spherical proton/neutron single-particle orbits */
  if(check_pn_orbits(data) == YES) {
    model->jbas[1] = model->jbas[0]; 
    free(data->jbas[1]);
  }
  else {
    model->jbas[1] = data->jbas[1];
  }
          /* copy data to gr_basis */

  gr_basis->P              = data->P;
  gr_basis->MJ             = data->MJ;
  gr_basis->part[0]        = data->part[0];
  gr_basis->part[1]        = data->part[1];
  gr_basis->mem_nondiag    = data->mem_nondiag;
  gr_basis->file_nondiag   = data->file_nondiag;

        /* filename to store diag/nondiag |SD> V_eff and J**2 matrix elements */

  strcpy(string,"store-");
  strcat(string, data->title);
  strcat(string,"-pn-SD-veff.dat");
  strcpy(gr_basis->store_veff, string);

  strcpy(string,"store-");
  strcat(string, data->title);
  strcat(string,"-pn-SD-ang.dat");
  strcpy(gr_basis->store_ang, string);

          /* copy data to iteration */

  strcpy(iteration->type_calc, data->type_calc);
  strcpy(iteration->out_data,  data->out_data);
  strcpy(iteration->type_calc, data->type_calc);
  strcpy(iteration->eigen_vec, data->eigen_vec);
  strcpy(iteration->start_vec, data->start_vec);
  strcpy(iteration->h_final,   data->h_final);

  if(   !strcmp(iteration->type_calc,"random-start")          /* calculation process */
     || !strcmp(iteration->type_calc,"random-continue")){
    iteration->tot_2J        = -1000;                         /* dummy values */
    iteration->num_start_vec = 0;
    iteration->vec_no        = NULL_PTR;
  }
  else if(   !strcmp(data->type_calc,"fixed-J-start") 
           || !strcmp(data->type_calc,"fixed-J-continue")) {
    iteration->tot_2J        = data->tot_2J;
    iteration->num_start_vec = data->num_start_vec;
    iteration->vec_no        = data->vec_no;
  }
  else {
    printf("\n\nError in function pn_data_structure():");
    printf("\nWrong Type of calculation process = %s\n",model->type_calc);  
    exit(1);
  }
  iteration->max_iterations = data->max_iterations;
  iteration->states         = data->states;

} /* End: function pn_data_structure() */

   /*
   ** The function
   **       check_pn_orbits()
   ** analyzes proton and neutron spherical single-particle orbits
   ** and return YES if they are identical, otherwise NO.
   */

  int check_pn_orbits(INPUT *data)
{
  int
       k;

  if(data->num_j[0] != data->num_j[1])  return NO;

  for(k = 0; k < data->num_j[0]; k++) {

    if(data->jbas[0]->osc      != data->jbas[1]->osc)      break;
    if(data->jbas[0]->l        != data->jbas[1]->l)        break;
    if(data->jbas[0]->j        != data->jbas[1]->j)        break;
    if(data->jbas[0]->min_part != data->jbas[1]->min_part) break;
    if(data->jbas[0]->max_part != data->jbas[1]->max_part) break;

    if(fabs(data->jbas[0]->e - data->jbas[1]->e) <MATRIX_LIMIT) break;
  }
  return ((k == data->num_j[0]) ? YES : NO);
} /* End: function check_pn_orbits() */  

   /*
   ** The function 
   **         print group_data()
   ** writes to output file all calculated group data.
   */

  void pn_print_group_data(char *file_name, GR_BASIS  *gr_basis, GROUP  *group)
{
  int
            k; 
  FILE
            *file_ptr;

  if( (file_ptr= fopen(file_name,"a")) == NULL) { /* open out_data file */
    printf("\n\nError in function pn_print_group_data()");
    printf("\nWrong file = %s for output data\n", file_name);
    exit(1);
  }
                     /******   Group heading  ****/

  fprintf(file_ptr,"\n\n     ******* Group data *******");
  fprintf(file_ptr,"\n  k  parZ  parN  2*M_Z  2*M_N   numSD_Z   num_SD_N      tot_SD");
  fprintf(file_ptr,"    sect");

  for(k = 0; k < gr_basis->num_gr; k++)  {
    fprintf(file_ptr,"\n%3d    %c    %c    %4d   %4d    %6d     %6d    %8d",
	                             k, ((group[k].par[0] == +1) ? '+' : '-'),
                                          ((group[k].par[1] == +1) ? '+' : '-'),
            group[k].m[0], group[k].m[1], group[k].numSD[0],group[k].numSD[1],
                                 group[k].numSD[0]*group[k].numSD[1]);
  }
  fprintf(file_ptr,"\n");
  fclose(file_ptr);
} /* End: function pn_print_group_data() */

     /*
     ** The function
     **         pn_result()
     ** takes the calculated result from the proton-neutron
     ** Lanczos process and write it to output file
     */

  void pn_result(PN_SHELL *model, GR_BASIS *gr_basis, GROUP *group, ITERATION *iteration)
{
  char
                   *func = {"pn_result(). "};
  register int
                   loop, loop_1;
  FILE
                   *file_ptr1, *file_ptr2;
  double
                   ang_mom, *proton_j, *neutron_j;
  LANC
                   vec;

       /* temporary memory to store the j-occupation numbers.*/

  proton_j  = MALLOC(model->num_j[0], double, func, "proton_j[]");
  neutron_j = MALLOC(model->num_j[1], double, func, "neutron_j[]");

               /* memory to store two lanczos vectors */

  vec.one = MALLOC(gr_basis->tot_numSD, double, func, "vec.one[]");
  vec.two = MALLOC(gr_basis->tot_numSD, double, func, "vec.two[]");

      /* Start the output procedure */

  if((file_ptr1 = fopen(model->out_data,"a")) == NULL)   {
    printf("\n\nError in function pn_result():");
    printf("\nWrong file = %s to open the output data file\n",model->out_data);
    exit(1);
  }
  switch(iteration->run_code) {
       case 0:
	       fprintf(file_ptr1,"\n\nThe Lanczos iteration process has reached ");
               fprintf(file_ptr1," maximum iterartions = %d - specified in input data\n",
                                                          iteration->max_iterations);
	       break;
       case 1:  
	       fprintf(file_ptr1,"\n\nThe Lanczos iteration process has converged ");
               fprintf(file_ptr1," after %d iterations.\n", iteration->num_iterations);
               break;
       case 2:  
	       fprintf(file_ptr1,"\n\nNo more Lanczos vectors");
               fprintf(file_ptr1," after %d iterations.", iteration->num_iterations);
               fprintf(file_ptr1, "\nNo energy convergence\n");
               break;
      default: fprintf(file_ptr1,"\n\nAbnormal termination of the Lanczos iteration");
               fprintf(file_ptr1," process - run_code = %d\n", iteration->run_code);
               exit(1);

  } /* end switch() loop */

  fprintf(file_ptr1,"\n\n\nFinal Eigenvalues");

      /** The total parity of the eigenstates **/

  if(model->P == +1)
    fprintf(file_ptr1,"\nThe total parity is positive\n");
   else
    fprintf(file_ptr1,"\nThe total parity is negative\n");

       /* open the eigenvector file for reading */

  if((file_ptr2 = fopen(iteration->eigen_vec,"rb")) == NULL)   {
    printf("\n\nError in function pn_result():");
    printf("\nWrong file = %s to read the eigenvectors\n",iteration->eigen_vec);
    exit(1);
  }
  for(loop = 0; loop < iteration->states; loop++)    {
    read_Lanczos_vector_from_open_file(iteration->eigen_vec, file_ptr2, loop, 
                                                   gr_basis->tot_numSD, vec.one);

    for(loop_1 = 0; loop_1 < gr_basis->tot_numSD; loop_1++) vec.two[loop_1] = D_ZERO;
    matrix_vector_multiplication(gr_basis, group, ANG, vec);
    ang_mom = scalar_product_of_vectors(gr_basis->tot_numSD, vec.one, vec.two);

                 /* Calculate the orbital j-occupation.*/ 

    for(loop_1 = 0; loop_1 < model->num_j[0]; loop_1++) proton_j[loop_1] = D_ZERO;
    proton_j_occupation(model, gr_basis, group, proton_j, vec.one);  
    
    for(loop_1 = 0; loop_1 < model->num_j[1]; loop_1++) neutron_j[loop_1] = D_ZERO;
    neutron_j_occupation(model, gr_basis, group, neutron_j, vec.one);

    fprintf(file_ptr1,"\n\nE(%d)= %9.4f", loop, iteration->eigenvalues[loop]);
    fprintf(file_ptr1,"  <J**2> = %7.4f", ang_mom);
    
    fprintf(file_ptr1,"\nProton single-part. distrib.  :");    /* proton occupation */
    for(loop_1 = 0; loop_1 < model->num_j[0] ; loop_1++) {
      fprintf(file_ptr1,"  %2d/2 ",model->jbas[0][loop_1].j);
    }
    fprintf(file_ptr1,"\n                          N(j): ");
    for(loop_1 = 0; loop_1 < model->num_j[0]; loop_1++) {
      fprintf(file_ptr1," %6.3f", proton_j[loop_1]);
    }
    fprintf(file_ptr1,"\nNeutron single-part. distrib. :");    /* neutron occupation */
    for(loop_1 = 0; loop_1 < model->num_j[1]; loop_1++) {
      fprintf(file_ptr1,"  %2d/2 ",model->jbas[1][loop_1].j);
    }
    fprintf(file_ptr1,"\n                          N(j): ");
    for(loop_1 = 0; loop_1 < model->num_j[1]; loop_1++) {
      fprintf(file_ptr1," %6.3f", neutron_j[loop_1]);
    }
  }  /* eigenvalues printed */

  fclose(file_ptr2);
  fclose(file_ptr1);

  free(vec.two);                      /* release local memory */
  free(vec.one);
  free(neutron_j);
  free(proton_j);

} /* End: function pn_result() */

    /*
    ** The function                            
    **        proton_j_occupation()                   
    ** calculate the j-occupation numbers  of protons for a shell model eigenstate 
    ** calculated through the Lanczo procedure and stored  in vector[]. 
    */

  void proton_j_occupation(PN_SHELL *model, GR_BASIS *gr_basis, GROUP *group,
                                               double *j_occ, double *vector)
{
  register int
                   orb, loop, part, curr_gr;
  register UL
                   pos, sd_val;
  UL
                   *proton_SD;
  double
                   probability, *amp_ptr;

  amp_ptr  = vector;               /* points to eigenvector amplitude */

     /* run through all proton-neutron groups */

  for(curr_gr = 0; curr_gr < gr_basis->num_gr; curr_gr++,  group++)  {

    if((group->numSD[0] * group->numSD[1]) == 0)   continue;  /* no basis states */

           /* run through all |SD(Z)> in current group */

    for(loop = 0, proton_SD = group->SD[0]; loop < group->numSD[0]; proton_SD++, loop++)  {

       /* run through all |SD(N)> in the group */

       for(orb = 0, probability = 0.0; orb < group->numSD[1]; amp_ptr++, orb++) { 
          probability += (*amp_ptr) * (*amp_ptr);
       }
       sd_val = *proton_SD;                     /* orbital initialization */
       pos    = UL_ONE;
       orb    = 0;
       part   = 0;
       do   {
           for( ;!(sd_val & pos); orb++, pos <<= 1);    /* occ. orbit found.*/
           j_occ[gr_basis->mbas[0][orb].orb] += probability;
           orb++;
           pos <<= 1;
       } while(++part < gr_basis->part[0]);      
    } /* loop through all proton SD in a group */
  }  /* loop through all groups */

} /* End: function proton_j_occupation() */

   /*
   ** The function                            
   **        neutron_j_occupation()                   
   ** calculate the j-occupation numbers  of neutrons for a shell model eigenstate 
   ** calculated through the lanczo procedure an stored  in Lanczo_one[]. 
   */

  void neutron_j_occupation(PN_SHELL *model, GR_BASIS *gr_basis, GROUP *group,
                                               double *j_occ, double *vector)
{
  register int
                  orb, loop, part, curr_gr;
  register UL
                  pos, sd_val;
  UL
                  *neutron_SD;
  double
                  probability, *group_ampl_ptr, *ampl_p_ptr, *ampl_n_ptr;

  group_ampl_ptr = vector;               /* points to the eigenvector */

     /* run through all proton-neutron groups */

  for(curr_gr = 0; curr_gr < gr_basis->num_gr; curr_gr++, group++)  {

    if((group->numSD[0] * group->numSD[1]) == 0)   continue;  /* no basis states */

    ampl_n_ptr = group_ampl_ptr;
    
    for(loop = 0, neutron_SD = group->SD[1]; loop < group->numSD[1]; 
                                            ampl_n_ptr++, neutron_SD++, loop++)   {
        ampl_p_ptr = ampl_n_ptr;  
        for(orb = 0, probability = 0.0; orb < group->numSD[0]; orb++,
                                  ampl_p_ptr += group->numSD[1])  {
          probability += (*ampl_p_ptr) * (*ampl_p_ptr);
        }
        sd_val = *neutron_SD;                     /* orbital initialization */
        pos    = UL_ONE;
        orb    = 0;
        part   = 0;
        do   {
           for( ;!(sd_val & pos); orb++, pos <<= 1);    /* occ. orbit found.*/
           j_occ[gr_basis->mbas[1][orb].orb] += probability;
           orb++;
           pos <<= 1; 
        } while(++part <  gr_basis->part[1]);      
     } /* loop through all neutron SD in a group */
     group_ampl_ptr += group->numSD[0] * group->numSD[1];
  }  /* loop through all groups */

} /* End: function neutron_j_occupation() */
