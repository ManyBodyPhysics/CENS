  /*
   *          Version  pr. April 2005        
   * Program module
   *                     transition-main.c
   * to calculate matrix elements of a one-body operator between shell 
   * model eigenstates generated in a lanczo procedure for both identical particles
   * and the proton-neutron case.
   * The initial and final states refer to the same set of single-particle states
   * with the same proton and neutron number, but the total parity may be different
   * The program module contains the following functions:
   *     1.   void main()
   *     2.   void basic_input_data()
   *     3.   void type_of_calculation()
   *     4.   void memory_allocation() 
   *     5.   int max_non_diag_elem()
   *     6.   void matrix_elements_in_m_scheme()
   *     7.   double m_scheme_e2_diag()
   *     8.   double m_scheme_e2_nondiag()
   *     9.   double rad_int2()
   *    10.   double m_scheme_m1_diag()
   *    11.   double m_scheme_m1_nondiag()
   *    12.   double m_scheme_e3_nondiag()
   *    13.   double rad_int_lambda()
   */


#include    "shell.h"
#include    "sm.h"

int main(int argc, char *argv[])
{
   char          file_name[ONE_LINE], total_string[ONE_LINE];
   struct tid    ex_time;

      /* Read file names and open file containing the basic shell model data.*/

   if(argc >= 2) {                              /* as program parameters */
      if( (Files.in_data = fopen(argv[1],"r")) == NULL) {
         printf("\nargc = %d Wrong file = %s for input data\n",argc, argv[1]);
         exit(1);
      }
   }
   else  {                             /* or as separate typed in filename */
      printf("\n Type in file name for basic shell model data =");
      scanf("%s",file_name);
      if( (Files.in_data = fopen(file_name,"r")) == NULL) {
         printf("\nWrong file = %s for the input data\n", file_name);
         exit(1);
      }
   } /* end of data-file-name input */

    /*
     * All data to specify the shell model problem and the transition calculations
     * in struct shell Model.
     */

   basic_input_data();

   if((Files.out_data= fopen(File_output,"w"))==NULL) {/* Open output data file*/
      printf("\nWrong file = %s for the output data\n", File_output);
      exit(1);
   }
    /* Copy the input file Files.in_data to the output file Files.out_data.*/

   rewind(Files.in_data);
   while(1)  {
      fgets(total_string,100, Files.in_data); /* read a line */
      if(!strncmp(total_string,"END_DATA",8)) break;
      fputs(total_string,Files.out_data);    /* print a line */
   } /* All data are copied */
   fprintf(Files.out_data,"\n                   ---    RESULTS   ----\n");

   if(fclose(Files.in_data)) { /* Close the input data file */

      fprintf(Files.out_data,"\n\n Error in closing the input data file;");
      fprintf(Files.out_data,"\n  %s\n",argv[1]);
      exit(1);
   }
   run_time(1,&ex_time);      /* start run time clock */

   /* 
    * Decide the type of calculation to be performed
    * and calculate the m-scheme single-particle orbits
    */

   type_of_calculation();

   if(Model[0].Z * Model[0].N)  { 
      pn_calculation();  /* proton - neutron calculation */
   }
   else   {
      id_calculation();      /* identical particle calculation */      
   }

   run_time(2,&ex_time);  /* stop and print run time clock */

    /* print the execution time */

   fprintf(Files.out_data,"\n run time: %d hour %d min %d sec",
          ex_time.hour, ex_time.min, ex_time.sec);

   fclose(Files.out_data);

   return 0;

} /* End: function main()  */

   /*
     * The function                              
     *         basic_input_data()                
     * reads in the set of total quantum numbers from a 
     * file and store the data in the struct shell Model. 
     */

void basic_input_data(void)
{  
  int x;   /* test */
  register int i;                                             /* loop index */
  register UL  test = UL_ZERO;                                /* control that all */
  register UL  control = UL_ZERO;                             /* data are read.   */
  int          number_of_data = 24;                           /* number of data lines  */
  static char  text_string[ONE_LINE], data_string[ONE_LINE], *end;
  char *ptr,   *identity[] =
               {"Initial proton number",
                "Final proton number",
                "Initial neutron number",
                "Final neutron number",
                "Initial total angular momentum J is (even, odd)",
                "Final total angular momentum J is (even, odd)",
                "Initial total projection of angular momentum (2*M)",
                "Final total projection of angular momentum (2*M)",
                "Initial total parity(+, -)",
                "Final total parity(+, -)",
                "The number of proton j-orbits",
                "The number of neutron j-orbits",
                "Storage file for initial eigenvectors",
                "Storage file for final eigenvectors",
                 "Output file for the results",
                "Type of calculation (e2, e3, m1)",
                "Number of transitions",
                "Init(no 2J) Final(no 2J) (ground state no = 0)",
		"E2 Neutron effective charge (unit e)",
		"E2 Proton effective charge (unit e)",
		"M1 Neutron effective gl (unit e)",
		"M1 Proton  effective gl (unit e)",
		"M1 Neutron effective gs (unit)",
		"M1 Proton  effective gs (unit)"};

    /* Read line from data-file and return it as a text-string and a data-string.*/

  while(read_line(text_string, data_string) != EOF)   {

      /* Identify text string and translate the corresponding data-string.*/

    for(i = 0; i < number_of_data; i++)  {
       if(strcmp(text_string,identity[i])) continue;
       x = i;
       break;   
    }
    switch (i)  {
      case 0:                                        /*  Initial proton number */
        Model[0].Z = (int) strtol(data_string,&end,10);
        test |= UL_ONE;
        break;
      case 1:                                        /* Final proton number */
        Model[1].Z = (int) strtol(data_string,&end,10);
        test |= (UL_ONE << 1);
        break;      
      case 2:                                        /* Initial neutron number */
        Model[0].N = (int) strtol(data_string,&end,10);
        test |= (UL_ONE << 2);
        break;
      case 3:                                        /* Final neutron number */
        Model[1].N = (int) strtol(data_string,&end,10);
        test |= (UL_ONE << 3);
        break;
      case 4:
        strcpy(Model[0].J_type, data_string);
        test |= (UL_ONE << 4);         
        break;
     case 5:
        strcpy(Model[1].J_type, data_string);
        test |= (UL_ONE << 5);         
        break;
      case 6:                          /* Initial total projection of angular momentum */
        Model[0].MJ = (int) strtol(data_string,&end,10);

        /*  Test relation between particle number and total m-value */
 
        if(MOD(Model[0].Z + Model[0].N,2) != MOD(Model[0].MJ,2))  {
           printf("\n\n Wrong value of initial total MJ value !!!!!!");
           printf("\nProton number = %d  Neutron number = %d  Total MJ = %d\n\n",
                               Model[0].Z, Model[0].N, Model[0].MJ);
           exit(1);
        }
        test |= (UL_ONE << 6);
        break;
      case 7:                          /*  Final total projection of angular momentum */
        Model[1].MJ = (int) strtol(data_string,&end,10);

        /*  Test relation between particle number and total m-value */
 
        if(MOD(Model[1].Z + Model[1].N,2) != MOD(Model[1].MJ,2))  {
           printf("\n\n Wrong value of final total MJ value !!!!!!");
           printf("\nProton number = %d  Neutron number = %d  Total MJ = %d\n\n",
                               Model[1].Z, Model[1].N, Model[1].MJ);
           exit(1);
        }
        test |= (UL_ONE << 7);
        break;
      case 8:                                        /* Initial total parity */
        ptr = data_string;
        do  {
          Model[0].P= *ptr++;
        } while(Model[0].P != '+'&& Model[0].P !='-'&& Model[0].P != '\n');
        if(Model[0].P == '\n')  {
          printf("\n\nError in initial total parity quantum number!\n");
          exit(1);
        }
        test |= (UL_ONE << 8);
        break;
      case 9:                                        /*  Final total parity */ 
        ptr = data_string;
        do  {
          Model[1].P= *ptr++;
        } while(Model[1].P != '+'&& Model[1].P !='-'&& Model[1].P != '\n');
        if(Model[1].P == '\n')  {
          printf("\n\nError in final total parity quantum number !\n");
          exit(1);
        }
        test |= (UL_ONE << 9);
        break;
    case 10:                                       /* The number of proton j-orbits */
        Model[0].numj_Z =(int) strtol(data_string,&end,10);
        test |= (UL_ONE << 10);
        memory_allocation(INPUT_PROTON_JBAS,0);
	read_line(text_string, data_string); /* extra text-line before the data */
        read_orbit(Model[0].numj_Z, Model[0].jbasZ);
        break;
    case 11:                                      /* The number of neutron j-orbits */
        Model[0].numj_N =(int) strtol(data_string,&end,10);
        test |= (UL_ONE << 11);
        memory_allocation(INPUT_NEUTRON_JBAS,0);
	read_line(text_string, data_string); /* extra text-line before the data */
        read_orbit(Model[0].numj_N, Model[0].jbasN);
        break;
    case 12:                                 /* Storage file for initial eigenvectors */
        strcpy(File_eigen_init, data_string);
	removeBlanks(File_eigen_init);
        test |= (UL_ONE << 12);
        break;
    case 13:   /* Storage file for final eigenvector */
        strcpy(File_eigen_final, data_string);
	removeBlanks(File_eigen_final);
        test |= (UL_ONE << 13);
        break;
    case 14:                                            /* Output file for the results */
        strcpy(File_output, data_string);
        test |= (UL_ONE << 14);
        break;
    case 15:            /* Type of calculation (e2, e3, m1) */
          strcpy(Model[0].calc, data_string);
          test |= (UL_ONE  << 15);
          break;
    case 16:                                            /* Number of transitions */
          Model[0].num_trans = (int) strtol(data_string,&end,10);
          test |= (UL_ONE << 16);
          break;
    case 17:                 /* Init(no 2J) Final(no 2J) */
          test |= (UL_ONE << 17);
          read_transitions();
          break;
    case 18:     // Neutron effective E2 charge
          Model[0].e2_n_eff = strtod(data_string,&end);
          test |= (UL_ONE << 18);
          break;
    case 19:   // Proton effective E2 charge
          Model[0].e2_p_eff = strtod(data_string,&end);
          test |= (UL_ONE << 19);
          break;
    case 20:    // Neutron effective M1 orbital gl charge
          Model[0].m1_n_eff_gl = strtod(data_string,&end);
          test |= (UL_ONE << 20);
          break;
    case 21:    // Proton effective M1 orbital gl charge
          Model[0].m1_p_eff_gl = strtod(data_string,&end);
          test |= (UL_ONE << 21);
          break;
    case 22:   // Neutron effective M1 spin gs charge
          Model[0].m1_n_eff_gs = strtod(data_string,&end);
          test |= (UL_ONE << 22);
          break;
    case 23:   // Proton effective M1 spin gs charge
          Model[0].m1_p_eff_gs = strtod(data_string,&end);
          test |= (UL_ONE << 23);
          break;
    default:
        break;
    } /* end of switch       */

    control = 0XFFFFF;
  } /* end - while : the whole input file has been read */

        /* check if the same SD basis are used for both the initial and final states */
     
   if(   (Model[0].Z == Model[1].Z) && (Model[0].N == Model[1].N) 
      && (Model[0].MJ == Model[1].MJ) && (Model[0].P == Model[1].P))
          Model[0].same_basis = YES;
   else
          Model[0].same_basis = NO;

     /* set up the control number for the input data.*/

  if((test & control) != control) {
    printf("\n\nNot enough shell model data in input file for the transfer calculation");
    printf("\ntest = %X  control = %X\n", test, control);
    exit(1);
  }
} /* End: function basic_input_data */

    /*
     * The function
     *     type_of_calculation()
     * prepare the initial and final calculation for protons only, neutrons only
     * or for both types of particles. Type of calculation is identified through 
     *       Info_com[].type = 0 - Model.Z > 0, Model.N = 0
     *                       = 1 - Model.Z = 0, Model.N > 0
     *                       = 2 - Model.Z = 1, Model.N = 1
     *                       = 3 - Model.Z = 1, Model.N > 1
     *                       = 4 - Model.Z > 1, Model.N = 1
     *                       = 5 - Model.Z > 1, Model.N > 1
     * where Info_com[0] (Info_com[1]) refers to initial (final) proton-neutron
     * configurations.
     */

void type_of_calculation(void)
{
  register int       loop;

  if((Model[0].Z > 0) || (Model[1].Z > 0))   { 
        /*
         * From the proton spherical orbits create a m-basis for the proton SD.
         * Store the result in structure m_orbit Z_orbit and write to output 
         * file a list of proton single-particle m-states.
         */
           
    for(Z_orbit.num = 0, loop = 0; loop < Model[0].numj_Z; loop++) {/* number of states.*/
      Z_orbit.num += Model[0].jbasZ[loop].j + 1;  /* 2*j+1 */
    }
        /* Check the proton number against number of m-orbits */
    
    if((Model[0].Z > Z_orbit.num) || (Model[1].Z > Z_orbit.num)) {
      printf("\nError from function type_of_calculation():");
      printf("\nNumber of initial protons = %d  Number of m_orbits = %d",
                                                      Model[0].Z, Z_orbit.num);
      printf("\nNumber of final protons = %d  Number of m_orbits = %d",
                                                      Model[1].Z, Z_orbit.num);
      exit(1);
    }
    memory_allocation(PROTON_M_ORBITS,0);

    fprintf(Files.out_data,
      "\n\n       PROTON single-particle orbits for m >= 1/2 (symmetric around zero)\n");
    single_m_states(Model[0].numj_Z, Model[0].jbasZ, &Z_orbit);

  } /* end calculation of proton single-particle orbits */

  if((Model[0].N > 0) || (Model[1].N > 0))  { 
        /*
         * From the neutron spherical orbits create a m-basis for the neutron SD.
         * Store the result in structure m_orbit N_orbit and write to output 
         * file a list of neutron single-particle m-states.
         */
           
    for(N_orbit.num = 0, loop = 0; loop < Model[0].numj_N; loop++) {/* number of states.*/
      N_orbit.num += Model[0].jbasN[loop].j + 1;  /* 2*j+1 */
    }
       /* Check the neutron number against number of m-orbits */
    
    if((Model[0].N > N_orbit.num) || (Model[1].N > N_orbit.num)) {
      printf("\nError from function type_of_calculation():");
      printf("\nNumber of initial neutrons = %d  Number of m-orbits = %d",
                                                         Model[0].N, N_orbit.num);
      printf("\nNumber of final neutrons = %d  Number of m-orbits = %d",
                                                         Model[1].N, N_orbit.num);
      exit(1);
    }
    memory_allocation(NEUTRON_M_ORBITS,0);

    fprintf(Files.out_data,
       "\n\n       NEUTRON single-particle orbits for m >= 1/2 (symmetric around zero)\n");

    single_m_states(Model[0].numj_N, Model[0].jbasN, &N_orbit);

  } /* end calculation of neutron single-particle orbits */

    /* run through initial (loop = 0) and final (loop = 1) case */

  for(loop = 0; loop < 2; loop++)  {

        /* Set up in  Info_com.type informations for the initial states.*/

     if((Model[loop].Z != 0)&&(Model[loop].N == 0)) Info_com[loop].type = 0;/* protons*/
     else if((Model[loop].Z == 0)&&(Model[loop].N != 0)) Info_com[loop].type = 1;/* neutrons*/
     else if(Model[loop].Z * Model[loop].N)  Info_com[loop].type = 2;  /* both particles */
     else  {
        printf("\n\nSomething is wrong with the initial particle numbers.");
        printf("\n proton number = %d  neutron number = %d",Model[loop].Z, Model[loop].N);
        exit(1);
     }
         /* basic data to Info_com and Info_id or Info_pn */

     Info_com[loop].par = (Model[loop].P == '+') ? +1: -1;
     Info_com[loop].MJ  = Model[loop].MJ;
     switch(Info_com[loop].type)  {
       case 0:     /* basic data for protons only.*/      
            Info_id[loop].part = Model[loop].Z;        /* proton particles.*/
            Info_id[loop].orb  = Z_orbit.num;    /* proton single-particle m-orbits.*/
            break;
       case 1:     /* basic data for neutrons only.*/      
            Info_id[loop].part = Model[loop].N;        /* neutron particles.*/
            Info_id[loop].orb  = N_orbit.num;    /* neutron single-particle m-orbits.*/
            break;
       case 2:     /* basic data for both protons and neutrons.*/
            Info_com[loop].type =  (Model[loop].Z == 1) 
                                 ? ((Model[loop].N == 1) ? 2 : 3) : ((Model[loop].N == 1) 
                                                                                ? 4 : 5); 
            Info_pn[loop].Z_part = Model[loop].Z;      /* proton particles.*/      
            Info_pn[loop].Z_orb  = Z_orbit.num;     /* proton single-particle m-orbits.*/
            Info_pn[loop].N_part = Model[loop].N;      /* neutron particles.*/
            Info_pn[loop].N_orb  = N_orbit.num;     /* neutron single-particle m-orbits.*/

     }  /* end switch */

  }    /* end loop */

}  /* End: function type_of_calculation() */

  /*
   * The function                                                      
   *                        memory_allocation(type, func)                              
   * is a general function to allocate permanent memory during the program execution
   * for both the initial(func = 0) and the final case(func = 1). 
   */

void memory_allocation(int type, int func)
{
  register int  loop, dim;

  switch(type)  {

    case INPUT_PROTON_JBAS :

      /* Memory for the proton single-particle orbits in a struct jbas. */

      Model[0].jbasZ = (struct jbas *) vector(sizeof(struct jbas),Model[0].numj_Z + 1);
      break;
    case INPUT_NEUTRON_JBAS :

      /* Memory for the neutron single-particle orbits in a struct jbas. */

      Model[0].jbasN = (struct jbas *) vector(sizeof(struct jbas),Model[0].numj_N + 1);
      break;
    case PROTON_M_ORBITS :

       /* Allocate memory to store the the proton single-particle m-orbits.*/
 
      Z_orbit.mbas = (struct mbas *) vector(sizeof(struct mbas), Z_orbit.num);

      //  Data and memory for masking  

      Z_orbit.numMask= 0;    // number of mask proton j_orbits
      for(loop = 0; loop <  Model[0].numj_Z; loop++) {
	if(  (Model[0].jbasZ[loop].min_part > 0) 
	   ||(Model[0].jbasZ[loop].max_part < Model[0].jbasZ[loop].j + 1))
          Z_orbit.numMask++;
      }
      if( Z_orbit.numMask > 0) {
	Z_orbit.lim  =  (int *) vector(sizeof(int), 2*Z_orbit.numMask);
	Z_orbit.list =  (ULL *) vector(sizeof(ULL), Z_orbit.numMask);
      }
      break;
    case NEUTRON_M_ORBITS:

       /* Allocate memory to store the the proton single-particle m-orbits.*/
 
      N_orbit.mbas = (struct mbas *) vector(sizeof(struct mbas), N_orbit.num);

      //  Data and memory for masking  

      N_orbit.numMask= 0;    // number of mask neutron j_orbits
      for(loop = 0; loop <  Model[0].numj_N; loop++) {
	if(  (Model[0].jbasN[loop].min_part > 0) 
	   ||(Model[0].jbasN[loop].max_part < Model[0].jbasN[loop].j + 1))
          N_orbit.numMask++;
      }
      if( N_orbit.numMask > 0) {
	N_orbit.lim  =  (int *) vector(sizeof(int), 2*N_orbit.numMask);
	N_orbit.list =  (ULL *) vector(sizeof(ULL), N_orbit.numMask);
      }
      break;
    case PROTON_NEUTRON_GROUPS:

      /* Memory to store information about each group (given M and P). */

      Group[func] = (struct group *) 
                    vector(sizeof(struct group), Info_pn[func].num_gr);
      break;
   case IDENTICAL_CONFIGURATION:
      
        /* Memory to store the identical SD configuration. */

      if(Info_com[func].tr_sym)  { 

            /* even Z (or N) with MJ = 0 */ 

         if((Model[func].same_basis == YES) && (Info_com[0].J_even || Info_com[1].J_even))  {
            SD_sym[func] = (ULL *) vector(sizeof(ULL), Info_id[func].num[1]);
         }
         else if(!Model[func].same_basis && Info_com[func].J_even)  {
             SD_sym[func] = (ULL *) vector(sizeof(ULL), Info_id[func].num[1]);
         }

         SD_asym[func]    = (ULL *) vector(sizeof(ULL), Info_id[func].num[0]); 
         SD_tr_asym[func] = (ULL *) vector(sizeof(ULL), Info_id[func].num[0]);

      } /* end time-reversal symmetry */

      else if(!Info_com[func].tr_sym)  {

             /* even Z (or N)  with MJ > 0 or odd Z (or N) */
     
         SD_asym[func]  = (ULL *) vector(sizeof(ULL), Info_id[func].num[0]);

      } /* end no time-reversal symmetry */

      else  {
         printf("\n\n Error: Wrong particle number or total MJ--value for case = %d.", func);
         printf("\nin function memory_allocation()");
         printf("\nInfo_id.part = %d  Model.MJ = %d", Info_id[func].part,Info_com[func].MJ);
         exit(1);
      }
      break;
   case PROTON_CONFIGURATION:

      /* Memory to store the  proton SD configuration. */

      SD_Z[func] = (ULL *) vector(sizeof(ULL),Info_pn[func].Z_num);
      break;
   case NEUTRON_CONFIGURATION:

      /* Memory to store the neutron SD configuration. */   

      SD_N[func] = (ULL *) vector(sizeof(ULL),Info_pn[func].N_num);
      break;

   case M1_E2_E3_MATRIX_ELEMENTS :

      if((Model[0].Z != Model[1].Z) || (Model[0].N != Model[1].N))  {
         fprintf(Files.out_data,"\n\nError from Memory_allocation()");
         fprintf(Files.out_data,"\n No changes in the proton and neutron numbers");
         fprintf(Files.out_data," are allowed during M1-E2_E3 transitions.");
         fprintf(Files.out_data,"\nModel[0].Z = %d Model[0].N = %d",Model[0].Z,Model[0].N);
         fprintf(Files.out_data,"\nModel[1].Z = %d Model[1].N = %d",Model[1].Z,Model[1].N);
         exit(1);
      }
      if(Model[0].Z)  { 
       
           /* permanent memory for the one-particle index table for proton particles */
 
               /* memory for proton one-particle index table */

         Table_pn_l[0] = (struct table_l *) vector(sizeof(struct table_l),(ULL) Z_orbit.num);

           /* memory for the off-diagonal proton single-particle matrix elements. */

         Trans_pn_nondiag[0] = (struct trans_nondiag *) 
                               vector(sizeof(struct trans_nondiag), max_non_diag_elem(0));
      }/* end proton particle case  */

      if(Model[1].N)  {
                /* memory for neutron one-particle index table */

         Table_pn_l[1] = (struct table_l *) vector(sizeof(struct table_l),(ULL) N_orbit.num);

           /* memory for the off-diagonal neutron single-particle matrix elements. */

         Trans_pn_nondiag[1] = (struct trans_nondiag *) 
                              vector(sizeof(struct trans_nondiag), max_non_diag_elem(1));

      } /* end neutron  particle case  */
      break;
  case TRANSITION_MATRIX_ELEMENTS :

             /* memory for <SD|OP|SD> for same basis */ 

       if(Model[0].same_basis == YES) 
          SD_diag = (double *) vector(sizeof(double),Info_com[0].num);

              /* memory allocation to store OP*|Eigen_vec > */

      Trans_contr = (double *) vector(sizeof(double), Info_com[1].num + 10);

          /* memory to store initial/final eigenvectors.*/
  
      dim = MAX(Info_com[0].num, Info_com[1].num);

      Eigen_vec = (double *) vector(sizeof(double), dim);

          /*
           * memory to store the one-particle transition matrix 
           * elements between the initial and final states.
           */

      Transitions = (double *) vector(sizeof(double),Model[0].num_trans);
      break;
  case GROUP_CHANGE :
             /* 
             ** memory to store group changes performed by the 
             ** one-particle operator OP(p,p') and OP(n'n')
             */

      Group_change_p = (int *) vector(Info_pn[0].num_gr, sizeof(int));
      Group_change_n = (int *) vector(Info_pn[0].num_gr, sizeof(int));
      break;
  } /* end switch-loop */

} /* End:  function memory_allocation() */

    /*
     * The function 
     *        max_non_diag_elem() 
     * calculates and returns the maximum possible number of non-diagonal one-body 
     * matrix elements for protons (par_type = 0) and neutrons (par_type = 1). 
     */

int max_non_diag_elem(int par_type)
{
  register int    right, left, m_trans, parity_trans, m_value, parity, count;
  struct m_orbit  *orbit;

  count   = 0;                                         /* initialization */
  orbit   = par_type ? &N_orbit : &Z_orbit;
  m_trans = Model[1].MJ - Model[0].MJ; 
  parity_trans 
    = ((!strncmp(Model[0].calc,"e2",2)) || (!strncmp(Model[0].calc,"m1",2))) ? +1 : -1;

            /* The right single-particle index */

  for(right = 0; right < orbit->num; right++) { 
     m_value = orbit->mbas[right].m;            /* right  m_value and parity */
     parity  = orbit->mbas[right].par;

          /* The left single-particle index  */

     for(left = 0 ; left < orbit->num; left++) {

        if(left == right) continue;

            /* test m_value and parity */

        if((   (m_value + m_trans) == orbit->mbas[left].m) 
            && ((parity * parity_trans) == orbit->mbas[left].par))   count++;

     } /* end of left index */
      count++; /* add space for an extra ZERO */
  } /* end of right index  */

  return count;

} /* End; function max_non_diag_elem() */

    /*
     * The function 
     *       matrix_elements_in_m_scheme() 
     * calculates the m-scheme single-particle matrix elements of a tensor operator
     * for both protons (part_type = 0) and neutrons (part_type = 1).
     * If initial and final SD's are equal the function calculates and store the 
     * diagonal matrix elements in double Trans_diag[]. 
     * The off-diagonal ones are calculated and stored in struct trans_nondiag 
     * Trans_non_diag[]. 
     * The struct table_l Table_l[] contains pointers to the matrix elements for each 
     * value of the right single-particle index. 
     */

void matrix_elements_in_m_scheme(int part_type, int trans_type)
{
   register int           left, tr_left, right, m_value, parity, num, m_trans;
   struct table_l         *ptr_diag;
   struct m_orbit         *orbit;
   struct trans_nondiag   *ptr_nondiag;

   orbit = part_type ? &N_orbit : &Z_orbit; 

   if(  (Model[0].same_basis == YES)
      && ((trans_type == E2_TRANSITION) || (trans_type == M1_TRANSITION)))  {

              /* calculate and store the diagonal single-particle matrix elements */

      for(left = 0, ptr_diag = Table_pn_l[part_type]; left < orbit->num; left++) {

              /* store diagonale matrix elements in Transition_diag[part_type][left] */
   
         switch(trans_type)  {
            case E2_TRANSITION: 
                    ptr_diag->diag = m_scheme_e2_diag(part_type, left,orbit->mbas);
                    break;
            case M1_TRANSITION: 
                    ptr_diag->diag = m_scheme_m1_diag(part_type, left,orbit->mbas);
                    break;
         } /* end switch() */

         ptr_diag++;

      } /* end the left loop */
   } /* end of the diagonal matrix elements */

     /* Calculate the off-diagonal transition matrix elements. */

   m_trans = Model[1].MJ - Model[0].MJ; 

         /* right single-particle index */

   ptr_diag = Table_pn_l[part_type];         /* initialization of the two pointers.*/
   ptr_nondiag = Trans_pn_nondiag[part_type];

   for(right= 0; right < orbit->num; right++, ptr_diag++) { 
      m_value          = orbit->mbas[right].m + m_trans; /* transferred m-value */
      parity           = orbit->mbas[right].par;
      ptr_diag->nondiag = ptr_nondiag;/* Store the pointer reference in the index table.*/

        /*
         * The left single-particle index. 
         * If initial and final SD's are equal then only non-diagonal matrix elements
         * with left single-particle index  > right single-particle index are stored.
         */

      num = 0;   /* initialization */

      for(left = 0, tr_left = orbit->num - 1 - left; left < orbit->num; left++, tr_left--) {

         if((left == right) || (m_value != orbit->mbas[left].m)) continue; /* test mj-value */

         switch(trans_type)  {
            case E2_TRANSITION: 
                if(parity != orbit->mbas[left].par) break; /* E2 parity = +1 */
                if(fabs(ptr_nondiag->val 
                          = m_scheme_e2_nondiag(part_type, left,right,orbit->mbas)) 
                                                                       > MATRIX_LIMIT) {
                   ptr_nondiag->one    = ULL_ONE << left;
                   ptr_nondiag->tr_one = ULL_ONE << tr_left;
                   ptr_nondiag->two    
                      = (ULL_ONE<<MAX(left, right)) - (ULL_ONE<<(MIN(left, right) +1));
                   num++;
                   ptr_nondiag++;
                }
                break;
	 case M1_TRANSITION:
                if(parity != orbit->mbas[left].par) break; /* M1 parity = +1 */
                if(fabs(ptr_nondiag->val 
                          = m_scheme_m1_nondiag(part_type, left,right,orbit->mbas)) 
                                                                       > MATRIX_LIMIT) {
                   ptr_nondiag->one    = ULL_ONE << left;
                   ptr_nondiag->tr_one = ULL_ONE << tr_left;
                   ptr_nondiag->two    
                      = (ULL_ONE<<MAX(left, right)) - (ULL_ONE<<(MIN(left, right) +1));
                   num++;
                   ptr_nondiag++;
                }
                break;
         case E3_TRANSITION: 
                if(parity == orbit->mbas[left].par) break; /* E3 parity = -1 */
                if(fabs(ptr_nondiag->val 
                          = m_scheme_e3_nondiag(part_type, left,right,orbit->mbas)) 
                                                                       > MATRIX_LIMIT) {
                   ptr_nondiag->one    = ULL_ONE << left;
                   ptr_nondiag->tr_one = ULL_ONE << tr_left;
                   ptr_nondiag->two    
                      = (ULL_ONE<<MAX(left, right)) - (ULL_ONE<<(MIN(left, right) +1));
                   num++;
                   ptr_nondiag++;
                }
         } /* end switch() loop */

      } /* end of left index */

          /* Terminate each left index group with a ZERO */

      if(num)   (ptr_nondiag++)->one = ULL_ZERO;
      else      ptr_diag->nondiag    =  NULL_PTR;

   } /* end of right index*/

} /* End: function matrix_elements_in_m_scheme() */

    /*
     * The function 
     *          m_scheme_e2_diag()
     calculates and return the diagonal single particle matrix elements of the
     * electric quadrupol operator
     *     <number | e_p * r**2 *  Y_2 | number>
     *           =   <n l j m|e_p * r**2 *  Y_2| n l j m> 
     */

double m_scheme_e2_diag(int part_type, int number,struct mbas *orbit)
{
  register int  n, l, j, m;
  double        value;

  value = D_ZERO;                /* initialization */

     /*  Single particle quantum numbers l, j, m */

  n = orbit[number].osc;
  l = orbit[number].l;
  j = orbit[number].j;   /*  twice the value of j  */
  m = orbit[number].m;   /*  and m.              */

  if((l == 0) || (j == 1))  return value; /* matrix element is ZERO */

  value  =  clebsch_gordan(2*l, 4,2*l, 0, 0) * clebsch_gordan(j, 4, j, m, 0)
            * 0.63078313 * rad_int2(n, l, n, l);

  if(j == (2*l + 1))  {
     value *= sqrt(((3.0 + j) * (j - 2.0)) / (j * (j + 1.0)));
  }
  else if(j == (2*l - 1))  {
     value *= sqrt(((4.0 + j) * (j - 1.0)) / ((j + 1.0) * (j + 2.0)));
  }
  else   {
        printf("\n\nError in function m_scheme_e2_diag():");
        printf("\n 2 * j = %d  and l = %d\n", j, l);
        exit(1);
  }

  return (value *(part_type ? Model[0].e2_n_eff : Model[0].e2_p_eff));

} /* End: function m_scheme_e2_diag() */

   /*
    * The function
    *          m_scheme_e2_nondiag()
    * calculates and return the non-diagonal single particle matrix elements 
    * of the  electric quadrupol operator
    *     <number | e_p * r**2 *  Y_2 | number>
    *                       =   <n l j m|e_p * r**2 *  Y_2| n l j m>
    * The single-particle states have the following phase conventions:
    *       a)  time-reversed phase factor: i^l 
    *       b)  a special "Whitehead" phase factor: (-1)^K_j
    *           where K_j = 1/2 MOD(2j + 1, 4)   for orbits with m_j > 0
    */

double m_scheme_e2_nondiag(int part_type, int left, int right, struct mbas *orbit)
{
  register int   n_f, l_f, j_f, m_f, n_i, l_i, j_i, m_i, phase;
  double         value;

     /*  left side single particle quantum numbers  */

  n_f   = orbit[left].osc;
  l_f   = orbit[left].l;
  j_f   = orbit[left].j;      /*  twice the value of j    */
  m_f   = orbit[left].m;      /*  and m_f                 */  
  phase = orbit[left].phase;  /* "Whitehead" phase factor */

      /*  right single particle quantum numbers */

  n_i   = orbit[right].osc;
  l_i   = orbit[right].l;
  j_i   = orbit[right].j;       /* twice the value of j     */
  m_i   = orbit[right].m;       /*  and m_i                 */
  phase *= orbit[right].phase;  /* "Whitehead" phase factor */

  if(PHASE(l_i + l_f) != +1) return D_ZERO;       /* (l_i + l_f) must be even */

  if(fabs(value = clebsch_gordan(2*l_i, 4,2*l_f, 0, 0)) < ZERO_LIMIT)  return D_ZERO;

      /* check total spin conservation */

  if(fabs(value *= clebsch_gordan(j_i, 4, j_f, m_i, m_f - m_i)) < ZERO_LIMIT) return D_ZERO;

  value *=  0.63078313 * sqrt((double)((2*l_i+1)*(j_i+1))) * rad_int2(n_f, l_f, n_i, l_i);

         /* 
          * Calculation of the 6J-symbol using the eqs. (6.3.3),(6.3.4) 
          * in Edmonds including the phase = (-)^(l_i + 1/2 + j_i)
          * and the time-reversed phase factor; 
          *      i^(l_i - l_f) = (-1)^(l_i - l_f)/2
          */

  if(j_i == (2*l_i + 1))   {
     if(j_f == (2*l_f + 1))   {
        value *=  PHASE((l_i - l_f)/2)
                  * sqrt(  ((6.0 + j_i + j_f) * (j_i + j_f - 4.0))
                         / (4.0 * j_i * (j_i + 1.0) * j_f * (j_f + 1.0)));
     }
     else if(j_f == (2*l_f - 1))   {
        value *= - PHASE((l_i - l_f)/2)
                   *  sqrt(  ((4.0 + j_i - j_f) * (6.0 + j_f - j_i))
                           / (4.0 * (j_f + 1.0) * (j_f + 2.0) * j_i * (j_i + 1.0)));
     }
     else   {
        printf("\n\nError in function m_scheme_e2_nondiag():");
        printf("\n 2 * j_f = %d  and l_f = %d\n", j_f, l_f);
        exit(1);
     } 
  }
  else if(j_i == (2*l_i - 1))   {
     if(j_f == (2*l_f + 1))   {
        value *= PHASE((l_i - l_f)/2)
                 *sqrt(  ((4.0 + j_f - j_i) * (6.0 + j_i - j_f))
                       / (4.0 * (j_i + 1.0) * (j_i + 2.0) * j_f * (j_f + 1.0)));
     }
     else if(j_f == (2*l_f - 1))   {
        value *= PHASE((l_i - l_f)/2)
                 * sqrt(  ((8.0 + j_f + j_i) * (j_f + j_i - 2.0))
                        / (4.0 * (j_f + 1.0) * (j_f + 2.0) * (j_i + 1.0) * (j_i + 2.0)));
     }
  }
  else   {
        printf("\n\nError in function m_scheme_e2_nondiag():");
        printf("\n 2 * j_i = %d  and l_i = %d\n", j_i, l_i);
        exit(1);
  } 
        /* Here we add the "Whitehead" phase factor: phase = K_j_i + K_j_f */
 
  return (value * phase * (part_type ? Model[0].e2_n_eff : Model[0].e2_p_eff));

} /* End: function m_scheme_e2_nondiag() */

     /*
      * The function  
      *           rad_int2( ) 
      * calculates the radial integral : Rnl's - within the same oscillator 
      * shell:  N = 2n + l
      *    < Nl | r**2 | Nl > = N + 3/2 
      * < N l-2 | r**2 | Nl > = 2sqrt{(n+1)(n+l+1/2)}
      */

double rad_int2(int n_1, int l_1, int n_2, int l_2)
{
  if((2*n_1+l_1) != (2*n_2+l_2))  return D_ZERO;
  if(l_1 == l_2)  return (double)(2*n_1 + l_1 + 1.5);
  else if((l_1 - l_2) == 2) return (2.0*sqrt((double)((n_1+1)*(n_1+l_1+0.5))));
  else if((l_2 - l_1) == 2) return (2.0*sqrt((double)((n_2+1)*(n_2+l_2+0.5))));
  else return D_ZERO;

}  /* End: function rad_int2() */

    /*
     * The function 
     *          m_scheme_m1_diag()
     * calculates and return the diagonal single particle matrix
     * elements of the magnetic dipole operator
     *     <number|g_l l_vec + gs s_vec|number>
     *          = <n l j m|g_l l_vec + gs s_vec|| n l j m>
     * Ref: Bohr and Mottelsen: Nuclear Structure I, p.336
     */


double m_scheme_m1_diag(int part_type, int number,struct mbas *orbit)
{

  int     n, l, j, m;
  double  factor = 0.488602512; // = sqrt(3/(4*pi))
  double  value;

     /*  Single particle quantum numbers l, j, m */

  n = orbit[number].osc;
  l = orbit[number].l;
  j = orbit[number].j;   /*  twice the value of j  */
  m = orbit[number].m;   /*  and m.              */

  value =  factor 
          *(clebsch_gordan(j,2,j,m,0)/clebsch_gordan(j, 2, j, j, 0)) 
          *((double)j/((double)(2*(2*l+1))));

  if(j == (2*l + 1)) {
    value *= part_type  
	      ? (2*l * Model[0].m1_n_eff_gl + Model[0].m1_n_eff_gs) 
             : (2*l * Model[0].m1_p_eff_gl + Model[0].m1_p_eff_gs); 
  }      
  else if(j == (2*l - 1)) {
    value *= part_type  
	    ? ((2*l + 2) * Model[0].m1_n_eff_gl - Model[0].m1_n_eff_gs) 
            : ((2*l + 2) * Model[0].m1_p_eff_gl - Model[0].m1_p_eff_gs); 
  }
  else   {
        printf("\n\nError in function m_scheme_m1_diag():");
        printf("\n Wrong: 2 * j = %d  and l = %d\n", j, l);
        exit(1);
  } 
	
  return value;

} /* End; function m_scheme_m1_diag() */

   /*
    * The function
    *          m_scheme_m1_nondiag()
    * calculates and return the non-diagonal single particle matrix elements 
    * of the  magnetic dipole  operator
    *     <left |g_l l_vec + gs s_vec| right>
    *    = <n_f l_f j_f m_f|g_l l_vec + gs s_vec|| n_i l_i j_i m_i>
    * The single-particle states have the following phase conventions:
    *       a)  time-reversed phase factor: i^l 
    *       b)  a special "Whitehead" phase factor: (-1)^K_j
    *           where K_j = 1/2 MOD(2j + 1, 4)   for orbits with m_j > 0
    * Ref: Bohr and Mottelsen: Nuclear Structure I, p.336
    */

double m_scheme_m1_nondiag(int part_type, int left, int right, struct mbas *orbit)
{
  register int   n_f, l_f, j_f, m_f, n_i, l_i, j_i, m_i, phase;

  double        factor = 0.488602512; // = sqrt(3/(4*pi))

  double         value;

     /*  left side single particle quantum numbers  */

  n_f   = orbit[left].osc;
  l_f   = orbit[left].l;
  j_f   = orbit[left].j;      /*  twice the value of j    */
  m_f   = orbit[left].m;      /*  and m_f                 */  
  phase = orbit[left].phase;  /* "Whitehead" phase factor */

      /*  right single particle quantum numbers */

  n_i   = orbit[right].osc;
  l_i   = orbit[right].l;
  j_i   = orbit[right].j;       /* twice the value of j     */
  m_i   = orbit[right].m;       /*  and m_i                 */
  phase *= orbit[right].phase;  /* "Whitehead" phase factor */

  if((n_i != n_f) || (l_i != l_f)) return D_ZERO;

  if((j_i == (2*l_i - 1)) && (j_f == (2*l_f + 1)))   {
    value = - factor*phase*sqrt((double)(j_i + 1)/((double)(j_f + 1))) 
              * clebsch_gordan(j_i, 2, j_f, 1, 0)
              * clebsch_gordan(j_i, 2, j_f, m_i, m_f - m_i);
  }             
  else if((j_i == (2*l_i + 1)) && (j_f == (2*l_f - 1)))   {
    value = factor*phase*clebsch_gordan(j_f, 2, j_i, 1, 0)
              * clebsch_gordan(j_i, 2, j_f, m_i, m_f - m_i);
  }
 else   {
        printf("\n\nError in function m_scheme_m1_nondiag():");
        printf("\n Wrong; 2 * j_i = %d  and l_i = %d\n", j_i, l_i);
        exit(1);
  } 

  return value * ( (part_type)  
		   ? (Model[0].m1_n_eff_gs - Model[0].m1_n_eff_gl) 
	           : (Model[0].m1_p_eff_gs - Model[0].m1_p_eff_gl));

} /* end: function m_scheme_m1_nondiag() */

   /*
    * The function
    *          m_scheme_e3_nondiag()
    * calculates and return the non-diagonal single particle matrix elements 
    * of the  electric octupole operator
    *     <number | e_p * r**3 *  Y_3 | number>
    *                       =   <n l j m|e_p * r**3 *  Y_3| n l j m>
    * The single-particle states have the following phase conventions:
    *       a)  time-reversed phase factor: i^l 
    *       b)  a special "Whitehead" phase factor: (-1)^K_j
    *           where K_j = 1/2 MOD(2j + 1, 4)   for orbits with m_j > 0
    */

double m_scheme_e3_nondiag(int part_type, int left, int right, struct mbas *orbit)
{
  register int   n_f, l_f, j_f, m_f, n_i, l_i, j_i, m_i, phase;
  double         value;

     /*  left side single particle quantum numbers  */

  n_f   = orbit[left].osc;
  l_f   = orbit[left].l;
  j_f   = orbit[left].j;      /*  twice the value of j    */
  m_f   = orbit[left].m;      /*  and m_f                 */  
  phase = orbit[left].phase;  /* "Whitehead" phase factor */

      /*  right single particle quantum numbers */

  n_i   = orbit[right].osc;
  l_i   = orbit[right].l;
  j_i   = orbit[right].j;       /* twice the value of j     */
  m_i   = orbit[right].m;       /*  and m_i                 */
  phase *= orbit[right].phase;  /* "Whitehead" phase factor */

  if(PHASE(l_i + l_f) != -1) return D_ZERO;       /* (l_i + l_f) must be odd */

  if(fabs(value = clebsch_gordan(2*l_i, 6,2*l_f, 0, 0)) < ZERO_LIMIT)  return D_ZERO;

      /* check total spin conservation */

  if(fabs(value *= clebsch_gordan(j_i, 6, j_f, m_i, m_f - m_i)) < ZERO_LIMIT) return D_ZERO;

  value *= 0.746352665 * sqrt((double)((2*l_i+1)*(j_i+1))) 
             * rad_int_lambda(n_f, l_f, n_i, l_i, 3);

         /* 
          * Calculation of the 6J-symbol using the eqs. (6.3.3),(6.3.4) 
          * in Edmonds including the phase = (-)^(l_f - 1/2 + j_i)
          * and the time-reversed phase factor; 
          *      i^(l_i - l_f) = (-1)^(l_i - l_f)/2
          */

  if(j_i == (2*l_i + 1))   {
     if(j_f == (2*l_f + 1))   {
        value *=  PHASE((l_i - l_f)/2)
                  * sqrt(  ((8.0 + j_i + j_f) * (j_i + j_f - 6.0))
                         / (4.0 * j_i * (j_i + 1.0) * j_f * (j_f + 1.0)));
     }
     else if(j_f == (2*l_f - 1))   {
        value *= PHASE((l_i - l_f)/2)
                   *  sqrt(  ((6.0 + j_i - j_f) * (8.0 + j_f - j_i))
                           / (4.0 * (j_f + 1.0) * (j_f + 2.0) * j_i * (j_i + 1.0)));
     }
     else   {
        printf("\n\nError in function m_scheme_e3_nondiag():");
        printf("\n 2 * j_f = %d  and l_f = %d\n", j_f, l_f);
        exit(1);
     } 
  }
  else if(j_i == (2*l_i - 1))   {
     if(j_f == (2*l_f + 1))   {
        value *= - PHASE((l_i - l_f)/2)
                 *sqrt(  ((6.0 + j_f - j_i) * (8.0 + j_i - j_f))
                       / (4.0 * (j_i + 1.0) * (j_i + 2.0) * j_f * (j_f + 1.0)));
     }
     else if(j_f == (2*l_f - 1))   {
        value *=  - PHASE((l_i - l_f)/2)
                 * sqrt(  ((10.0 + j_f + j_i) * (j_f + j_i - 4.0))
                        / (4.0 * (j_f + 1.0) * (j_f + 2.0) * (j_i + 1.0) * (j_i + 2.0)));
     }
  }
  else   {
        printf("\n\nError in function m_scheme_e3_nondiag():");
        printf("\n 2 * j_i = %d  and l_i = %d\n", j_i, l_i);
        exit(1);
  } 
        /* Here we add the "Whitehead" phase factor: phase = K_j_i + K_j_f */
 
  return (value * phase * (part_type ? Model[0].e2_n_eff : Model[0].e2_p_eff));

} /* End: function m_scheme_e2_nondiag() */

     /*
      * The function  
      *           rad_int_lambda( ) 
      * calculates the radial integral 
      *         <n_2,l_2 | r**(lambda) |n_1, l_1>
      * Reference: K. Heyde: The Nuclear Shell Model,  
      * eqn.(9.18), page 320.
      */

double rad_int_lambda(int n_1, int l_1, int n_2, int l_2, int lambda)
{
   register int  nn_1, nn_2, sigma, low, high;
   double        sum;

        /* the factor 10**(lambda/2) for lambda <= 6 */

   static double factor[7] = {1.0, 3.16227766, 10.0, 31.62277660, 100.0, 316.22776602, 1000.0};

         /* fct(K)= gamma(K/2) / 10**(K/2 -1) , K = 0,....,50. */

   static double fct[51] = { 
            0.0,           5.60499122,    1.0,           2.80249561e-1, 1.0e-1, 
            4.20374341e-2, 2.0e-2,        1.05093585e-2, 6.0e-3,        3.67827549e-3,
            2.40e-3,       1.65522397e-3, 1.20e-3,       9.10373183e-4, 7.2e-4,
            5.91742569e-4, 5.04e-4,       4.43806927e-4, 4.032e-4,      3.77235888e-4,
            3.6288e-4,     3.58374093e-4, 3.6288e-4,     3.76292798e-4, 3.99168e-4,
            4.32736718e-4, 4.790016e-4,   5.40920897e-4, 6.2270208e-4,  7.30243211e-4,
            8.71782912e-4, 1.05885266e-3, 1.30767437e-3, 1.64122162e-3, 2.09227899e-3,
            2.70801567e-3, 3.55687428e-3, 4.73902742e-3, 6.40237371e-3, 8.76720072e-3,
            1.21645100e-2, 1.70960414e-2, 2.43290201e-2, 3.50468849e-2, 5.10909422e-2,
            7.53508025e-2, 1.12400073e-1, 1.69539306e-1, 2.58520167e-1, 3.98417368e-1,
            6.20448402e-1};
 

   if(lambda > 6)   {
      printf("\n\nError in function rad_int_lambda():");
      printf("  Not allowed lambda-value!!! lambda  = %d\n", lambda);   
      exit(1);
   }
   if((l_2 > (l_1 + lambda)) || (l_2 < abs(l_1 - lambda)))   {
      printf("\n\nError in function rad_int_lambda():");
      printf(" Not allowed angular momentum relation");
      printf("\n<n_2 = %d, l_2 = %d | r**%d | n_1 = %d, l_1 = %d >\n",
                                           n_2, l_2, lambda, n_1, l_1);    
      exit(1);
   }
   nn_1  = 2 * n_1;
   nn_2  = 2 * n_2;
   low   = MAX(MAX(nn_1 - l_2 + l_1 - lambda, nn_2  - l_1 + l_2 - lambda), 0);
   high  = MIN(nn_1, nn_2);
   for(sum = 0.0, sigma = low; sigma <= high; sigma += 2)   {
      sum +=    fct[l_1 + l_2 + lambda + 3 + sigma]
             / (  fct[2 + sigma] 
                * fct[nn_1 + 2 - sigma] 
                * fct[nn_2  + 2- sigma]
                * fct[l_2 - l_1 + lambda + 2 - nn_1 + sigma]
                * fct[l_1 - l_2 + lambda + 2 - nn_2 + sigma]);
   }
   sum *= sqrt(  (fct[nn_1 + 2] * fct[nn_2 + 2])
               / (fct[nn_1 + 3 + 2 * l_1] * fct[nn_2 + 3 + 2 * l_2]))
               * fct[l_2 - l_1 + lambda + 2] 
               * fct[l_1 - l_2 + lambda + 2] 
               * factor[lambda]; 

   return sum;
} /* End: function rad_int_lambda() */
