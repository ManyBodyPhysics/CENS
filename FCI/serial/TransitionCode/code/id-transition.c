      /*
      **  Program module
      **            id-transition.c
      ** calculates different one-particle transition
      ** for identical particles between eigenstates
      ** generated in a Lanczos procedure.
      ** The program module contains the following functions:
      **     1.   id_calculation()
      **     2.   void id_slater_determinant()
      **     3.   void number_of_identical_SD()
      **     4.   void number_of_SD_tr_sym()
      **     5.   void number_of_SD_no_sym()
      **     6.   void identical_SD_configuration()
      **     7.   void identical_SD_tr_sym()
      **     8.   void identical_SD_no_sym()
      **     9.   void id_particle_transitions()
      **    10.   void id_diag_trans()
      **    11.   void read_id_eigenvector()
      **    12.   void  id_trans_contribution()
      **    13.   void diag_contr()
      **    14.   void id_nondiag_asymI()
      **    15.   void id_nondiag_symI()
      **    16.   void id_nondiag_asymII()
      **    17.   void id_nondiag_symII()
      **    18.   void id_nondiag_nosymIII()
      **    19.   void id_nondiag_nosymIV()
      **    20.   double overlap()
      */

#include "shell.h"
#include "sm.h"

void id_calculation(void)
{
   int     part_type;   /* = 0 (1) for protons (neutrons) */

   part_type = Model[0].Z ? 0 : 1;

     /* Calculate and store the initial and final SD. */

   id_slater_determinant(); 

    /*
     * Choice of three different types of calculations: 
     *    1. E2 transitions between eigenstates given in a SD basis.         
     *    2. M1 transitions between eigenstates given in a SD basis.         
     *    3. E3 transitions between eigenstates given in a SD basis.
     * In all cases the proton-neutron numbers are the same in 
     * initial and final states. 
     */

   if(!strncmp(Model[0].calc,"e2",2)) {

        /*
         * Allocate memory for the index-l table, the diagonal and off-diagonal 
         * single-particle  matrix elements.
         */

      memory_allocation(M1_E2_E3_MATRIX_ELEMENTS,0); 

      matrix_elements_in_m_scheme(part_type, E2_TRANSITION);

         /* Allocate memory for the calculation of the transition matrix elements.*/

      memory_allocation(TRANSITION_MATRIX_ELEMENTS,0);

         /* Calculate the E2 transitions.*/

      id_particle_transitions(E2_TRANSITION); 

   }
   if(!strncmp(Model[0].calc,"m1",2)) {

        /*
         * Allocate memory for the index-l table, the diagonal and off-diagonal 
         * single-particle  matrix elements.
         */

      memory_allocation(M1_E2_E3_MATRIX_ELEMENTS,0); 

      matrix_elements_in_m_scheme(part_type, M1_TRANSITION);

         /* Allocate memory for the calculation of the transition matrix elements.*/

      memory_allocation(TRANSITION_MATRIX_ELEMENTS,0);

         /* Calculate the E2 transitions.*/

      id_particle_transitions(M1_TRANSITION); 

   }
   if(!strncmp(Model[0].calc,"e3",2))  {

        /*
         * Allocate memory for the index-l table, the diagonal and off-diagonal 
         * single-particle  matrix elements for both protons and neutrons. 
         */

      memory_allocation(M1_E2_E3_MATRIX_ELEMENTS,0); 

      matrix_elements_in_m_scheme(part_type, E3_TRANSITION);

         /* Allocate memory for the calculation of the transition matrix elements.*/

      memory_allocation(TRANSITION_MATRIX_ELEMENTS,0);

         /* Calculate the E3 transitions.*/

      id_particle_transitions(E3_TRANSITION); 

   }
 }  /* End: function  id_calculation() */

  /*
   * The function                                                      
   *            void id_slater_determinant(func) 
   * allocates memory, calculates and store the complete set of slater 
   * determinants (SD) for a system of identical particles for both
   * the initial and the final eigenstates.       
   */

void id_slater_determinant(void)
{
   register int   sym_conf;

   Info_com[0].J_even = !strcmp(Model[0].J_type,"even") ? YES : NO; /* initialization */
   Info_com[1].J_even = !strcmp(Model[1].J_type,"even") ? YES : NO;
   Info_com[0].tr_sym = ((Info_com[0].MJ == 0 ) && !MOD(Info_id[0].part,2)) ? YES : NO;
   Info_com[1].tr_sym = ((Info_com[1].MJ == 0 ) && !MOD(Info_id[1].part,2)) ? YES : NO;

   sym_conf = NO; 
   if(Model[0].same_basis)  {    
      if(Info_com[0].tr_sym)   {
         if(Info_com[0].J_even || Info_com[1].J_even) sym_conf = YES;

         number_of_identical_SD(0, sym_conf);       /* initial and final basis basis */
         memory_allocation(IDENTICAL_CONFIGURATION,0);
         identical_SD_configuration(0, sym_conf);

         Info_id[1].num[0]     = Info_id[0].num[0];
         Info_id[1].num[1]     = Info_id[0].num[1];
         Info_com[0].num       = Info_id[0].num[0] 
                                 + (Info_com[0].J_even ? Info_id[0].num[1] : 0);
         Info_com[1].num       = Info_id[1].num[0] 
                                 + (Info_com[1].J_even ? Info_id[1].num[1] : 0);

         SD_asym[1]    = SD_asym[0];      /* initial and final SD basis is equal */
         SD_tr_asym[1] = SD_tr_asym[0];
         if(sym_conf == YES) SD_sym[1] = SD_sym[0]; 
      }
      else  {

         number_of_identical_SD(0, sym_conf);       /* initial and final basis */
         memory_allocation(IDENTICAL_CONFIGURATION,0);
         identical_SD_configuration(0, sym_conf);

         Info_id[1].num[0]     = Info_id[0].num[0];
         Info_com[0].num       = Info_id[0].num[0];
         Info_com[1].num       = Info_id[1].num[0];

         SD_asym[1] = SD_asym[0];      /* initial and final SD basis is equal */
      }

   } /* end same initial and final basis */
       
   else if(!Model[0].same_basis)  {  
      if(Info_com[0].tr_sym)  {    /* initial basis */
         if(Info_com[0].J_even) sym_conf = YES;

         number_of_identical_SD(0, sym_conf);           /* initial basis */
         memory_allocation(IDENTICAL_CONFIGURATION,0);
         identical_SD_configuration(0, sym_conf);
    
         Info_com[0].num = Info_id[0].num[0] + Info_id[0].num[1];
      }
      else   {
         number_of_identical_SD(0, sym_conf);           /* initial basis */
         memory_allocation(IDENTICAL_CONFIGURATION,0);
         identical_SD_configuration(0, sym_conf);

         Info_com[0].num = Info_id[0].num[0];
      }  /* end initial basis */

      sym_conf = NO;
      if(Info_com[1].tr_sym)  {    /* final basis */
         if(Info_com[1].J_even) sym_conf = YES;

         number_of_identical_SD(1, sym_conf);           /* final basis */
         memory_allocation(IDENTICAL_CONFIGURATION,1);
         identical_SD_configuration(1, sym_conf);
    
         Info_com[1].num = Info_id[1].num[0] + Info_id[1].num[1];
      }
      else   {
         number_of_identical_SD(1, sym_conf);           /* final basis */
         memory_allocation(IDENTICAL_CONFIGURATION,1);
         identical_SD_configuration(1, sym_conf);

         Info_com[1].num = Info_id[1].num[0];

      }  /* end final basis */

   } /* end different initial and final basis */   

} /* End: function id_slater_determinant() */

  /*
   * The function
   *          number_of_identical_SD(func)
   * calculates the number of SD for the case of identical particles
   *         Info_com[func].type = 0  -- protons only
   *                       .type = 1  -- neutrons only
   * for the initial(func = 0) and the final case(func = 1).
   * If possible time-reversed symmetry is assumed.
   */

void number_of_identical_SD(int func, int sym_conf)
{
  if(Info_com[func].type == 0)  { /* protons only */
     if(Info_com[func].tr_sym)  {
        number_of_SD_tr_sym(func, sym_conf, Z_orbit);  /* even Z with MJ = 0 */
        if(func == 0) {
           fprintf(Files.out_data,
              "\n\nNumber of initial asymmetric proton slater determinants = %8d",
                                                               Info_id[func].num[0]);
           fprintf(Files.out_data,
              "\nNumber of  initial symmetric proton slater determinants = %8d",
                                                               Info_id[func].num[1]);
           fprintf(Files.out_data,
               "\nTotal number of initial slater determinants              = %8d",
                                   2 * Info_id[func].num[0] +  Info_id[func].num[1]);
        }
        else  {
           fprintf(Files.out_data,
              "\n\nNumber of final asymmetric proton slater determinants = %8d",
                                                               Info_id[func].num[0]);
           fprintf(Files.out_data,
              "\nNumber of  final symmetric proton slater determinants = %8d",
                                                               Info_id[func].num[1]);
           fprintf(Files.out_data,
               "\nTotal number of final slater determinants              = %8d",
                                   2 * Info_id[func].num[0] +  Info_id[func].num[1]);
        }
     }
     else if(!Info_com[func].tr_sym)   {
        number_of_SD_no_sym(func, Z_orbit);  /* even Z with MJ > 0 or odd Z */
        if(func == 0)   {
           fprintf(Files.out_data,
              "\n\nNumber of initial non-symmetric proton slater determinants = %8d",
                                                                  Info_id[func].num[0]);
        }
        else    {
            fprintf(Files.out_data,
              "\n\nNumber of final non-symmetric proton slater determinants = %8d",
                                                                 Info_id[func].num[0]);
        }
     }
     else  {
        printf("\n\nError: Wrong particle number or total MJ-value for case = %d.",func);
        printf("\nin function number_of_identical_SD()");
        printf("\nModel.Z = %d Model.MJ = %d",
                  Info_id[func].part,Info_com[func].MJ);
        exit(1);
     }
  } /* end proton case */

  if(Info_com[func].type == 1)  { /* neutrons only */
     if(Info_com[func].tr_sym)  {
        number_of_SD_tr_sym(func, sym_conf, N_orbit);  /* even N with MJ = 0 */
        if(func == 0) {
           fprintf(Files.out_data,
              "\n\nNumber of initial asymmetric neutron slater determinants = %8d",
                                                               Info_id[func].num[0]);
           fprintf(Files.out_data,
              "\nNumber of  initial symmetric neutron slater determinants = %8d",
                                                               Info_id[func].num[1]);
           fprintf(Files.out_data,
               "\nTotal number of initial slater determinants              = %8d",
                                   2 * Info_id[func].num[0] +  Info_id[func].num[1]);
        }
        else  {
           fprintf(Files.out_data,
              "\n\nNumber of final asymmetric neutron slater determinants = %8d",
                                                               Info_id[func].num[0]);
           fprintf(Files.out_data,
              "\nNumber of  final symmetric neutron slater determinants = %8d",
                                                               Info_id[func].num[1]);
           fprintf(Files.out_data,
               "\nTotal number of final slater determinants              = %8d",
                                   2 * Info_id[func].num[0] +  Info_id[func].num[1]);
        }
     }
     else if(!Info_com[func].tr_sym)   {
        number_of_SD_no_sym(func, N_orbit);  /* even N with MJ > 0 or odd N */
        if(func == 0)   {
           fprintf(Files.out_data,
              "\n\nNumber of initial non-symmetric neutron slater determinants = %8d",
                                                                  Info_id[func].num[0]);
        }
        else    {
            fprintf(Files.out_data,
              "\n\nNumber of final non-symmetric neutron slater determinants = %8d",
                                                                 Info_id[func].num[0]);
        }
     }
     else  {
        printf("\n\nError: Wrong particle number or total MJ-value for case = %d.",func);
        printf("\nin function number_of_identical_SD()");
        printf("\nInfo_id.part = %d Model.MJ = %d",Info_id[func].part,Info_com[func].MJ);
        exit(1);
     }
  } /* end neutron case */
   
} /* End: function number_of_identical_SD() */

  /*
   * The function                                           
   *          number_of_SD_tr_sym(func)                   
   * calculates the number of slater determinants (SD) for even number of identical 
   * particles with MJ = 0 for both initial(func = 0) and finale case(func = 1).
   * It calculates the number of asymmetric (Info_id.num[0]) and the number of 
   * symmetric (Info_id.num[1] SD).
   */

void number_of_SD_tr_sym(int func, int sym_conf, struct m_orbit orbit)
{
  register int     loop,m, parity, last_orbit, num, phase;
  register ULL      config, time_reversed_config, high;
  int              *matr;

  matr = (int *) vector(sizeof(int), Info_id[func].part + 1); /* temporary memory */

      /* the lowest  configuration */

  for(loop = 0; loop < Info_id[func].part; loop++) matr[loop] = loop;
  matr[loop] = orbit.num;

       /* calculate the total number of SD. */

  Info_id[func].num[0] = 0;  /* initialization */
  Info_id[func].num[1] = 0;
  last_orbit = orbit.num - 1;

  do  {  /* run through all possible configurations */ 
    
       /* find m-value and parity */

    for(m = 0, parity = +1, loop = 0; loop < Info_id[func].part; loop++)  {
       m += orbit.mbas[matr[loop]].m;
       parity *= orbit.mbas[matr[loop]].par;
    }
    if((m == Info_com[func].MJ) && (parity == Info_com[func].par))  {
    
      // Code config  and its time-reversed

      config               = ULL_ZERO;
      time_reversed_config = ULL_ZERO;
      for(loop = 0; loop < Info_id[func].part; loop++)  {
         config               += ULL_ONE << matr[loop];
         time_reversed_config += ULL_ONE <<(last_orbit - matr[loop]);
      }

      if( config <= time_reversed_config) {
	if(config < time_reversed_config) loop = 0;
        else if((config == time_reversed_config) && (sym_conf == YES)) loop = 1;
        else loop = 2;

	if((orbit.numMask == 0) && (loop < 2)) { 
	  Info_id[func].num[loop]++;
	}
	else if((orbit.numMask > 0) && (loop < 2)) {
	  for(num = 0; num < orbit.numMask; num++) {
	    high = config & orbit.list[num];
	    phase = 0;
	    while(high) {
	      high &= high - 1;
	      phase++;
	    }
	    if(  (phase < orbit.lim[2*num])
		 ||(phase > orbit.lim[2*num + 1])) break;
	  }
	  if(num == orbit.numMask) Info_id[func].num[loop]++;
	}
      } // end config < time_reversed_config

    } // end m- and parity test

      /* new particle configuration */

      for(loop = 0; (loop < Info_id[func].part) && ((matr[loop]+1) >= matr[loop+1]); loop++);
      matr[loop]++;
      for(loop--; loop >= 0; loop--)  matr[loop] = loop;

  } while(matr[Info_id[func].part] == orbit.num);  // end all config

  free(matr);  // release temporary memory

  return;

} /* End: function number_of_SD_tr_sym() */

  /*
   * The function                                           
   *          number_of_SD_no_sym_(func)                   
   * calculates the number of slater determinants (SD) for identical particles without
   * time-reversal symmetry for both the initial(func = 0) and the final case(func = 1).
   * It calculates the number of  SD (Info_id.num[0]). (Info_num[1] = 0).
   */

void number_of_SD_no_sym(int func, struct m_orbit orbit)
{
  register int     loop, m, parity, num, phase;
  int              *matr;
  ULL              config, high;

  matr = (int *) vector(sizeof(int), Info_id[func].part + 1); /* temporary memory */

      /* the lowest neutron configuration */

  for(loop = 0; loop < Info_id[func].part; loop++) matr[loop] = loop;
  matr[loop] = orbit.num;

        /* calculate the total number of SD.*/

  Info_id[func].num[0] = 0;   /* initialization */
  Info_id[func].num[1] = 0;  

  do  {    /* run through all configurations */

        /* find m-value and parity */

    for(m = 0, parity = +1, loop = 0; loop < Info_id[func].part; loop++)  {
       m      += orbit.mbas[matr[loop]].m;
       parity *= orbit.mbas[matr[loop]].par;
    }

    if((m == Info_com[func].MJ) && (parity == Info_com[func].par)) {
      if(orbit.numMask == 0) { 
	 Info_id[func].num[0]++;
      }
      else  {   // check for particle limitations
  	        // code the configuration
	for(config = ULL_ZERO, loop = 0; loop < Info_id[func].part; loop++)  {
	  config += ULL_ONE << matr[loop];
	}
	for(num = 0; num < orbit.numMask; num++) {
	  high = config & orbit.list[num];
	  phase = 0;
	  while(high) {
	    high &= high - 1;
	    phase++;
	  }
	  if(  (phase < orbit.lim[2*num])
	       ||(phase > orbit.lim[2*num + 1])) break;
	}
	if(num == orbit.numMask) Info_id[func].num[0]++;
      }
    }
      /* new particle configuration */

    for(loop = 0; (loop < Info_id[func].part) && ((matr[loop]+1) >= matr[loop+1]); loop++);
    matr[loop]++;
    for(loop--; loop >= 0; loop--) matr[loop] = loop;

  } while(matr[Info_id[func].part] == orbit.num);  /* end running through all configurations */

  free(matr);  /* release temporary memory */

  return;

} /* End: function number_of_SD_no_sym() */

  /*
   * The function
   *          identical_SD_configuration(func)
   * calculates and stores the complete set  of SD for the case of identical particles
   *         Info_com.type = 0  -- protons only
   *                 .type = 1  -- neutrons only
   * for both the initial(func = 0) and the final case(func = 1).
   * If possible time-reversed symmetry is assumed.
   */

void identical_SD_configuration(int func, int sym_conf)
{
  if(Info_com[func].type == 0)  { /* protons only */

     if(Info_com[func].tr_sym)

        identical_SD_tr_sym(func, sym_conf, Z_orbit); /* even Z with MJ = 0 */

     else if(!Info_com[func].tr_sym)

        identical_SD_no_sym(func, Z_orbit);   /* even Z with MJ > 0 or odd Z */

     else  {
        printf("\n\nError: Wrong particle number or total MJ-value for case = %d.",func);
        printf("\nin function identical_SD_configuration()");
        printf("\nModel.Z = %d Model.MJ = %d\n",Info_id[func].part,Info_com[func].MJ);
        exit(1);
     }
  } /* end proton case */

  if(Info_com[func].type == 1)  { /* neutrons only */

     if(Info_com[func].tr_sym)

        identical_SD_tr_sym(func, sym_conf, N_orbit); /* even N with MJ = 0 */

     else if(!Info_com[func].tr_sym)

        identical_SD_no_sym(func, N_orbit);  /* even N with MJ > 0 or odd Z */

     else  {
        printf("\n\nError: Wrong particle number or total MJ-value for case = %d .",func);
        printf("\nin function identical_SD_configuration()");
        printf("\nModel.N = %d Model.MJ = %d",Info_id[func].part,Info_com[func].MJ);
        exit(1);
     }
  } /* end neutron case */
   
} /* End: function identical_SD_configuration() */

  /*
   * The function                                           
   *          identical_SD_tr_sym(func)                   
   * calculates and stores the complete set of slater determinants (SD) for identical
   * particles with MJ = 0 for both initial(func = 0) and the final case(func = 1).
   * It calculates both the asymmetric SD and the symmetric SD. 
   */

void identical_SD_tr_sym(int func, int sym_conf, struct m_orbit orbit)
{
  register int     loop, m, parity, last_orbit,num, phase,count_asym, count_sym;
  register ULL      sd_asym, sd_tr_asym, high;
  int              *matr;
     
  matr = (int *) vector(sizeof(int), Info_id[func].part + 1); /* temporary memory */

      /* the lowest neutron configuration */

  for(loop = 0; loop < Info_id[func].part; loop++) matr[loop] = loop;
  matr[loop] = orbit.num;

       /* calculate the total number of SD.*/

  last_orbit = orbit.num - 1;
  count_asym  = 0;
  count_sym   = 0;

  do  {  /* run through all configuration */
    
      /* find m-value and parity */

    for(m = 0, parity = +1, loop = 0; loop < Info_id[func].part; loop++)  {
       m      += orbit.mbas[matr[loop]].m;
       parity *= orbit.mbas[matr[loop]].par;
    }
    if((m == Info_com[func].MJ) && (parity == Info_com[func].par))  {

        /* Code the configuration  and its time-reversed one */

      sd_asym    = 0;
      sd_tr_asym = 0;
      for(loop = 0; loop < Info_id[func].part; loop++)  {
         sd_asym    += ULL_ONE << matr[loop];
         sd_tr_asym += ULL_ONE <<(last_orbit - matr[loop]);
      }
      if(sd_asym < sd_tr_asym)           loop = 0;
      else if(   (sd_asym == sd_tr_asym) 
	      && (sym_conf == YES))      loop = 1;
      else                               loop = 2;

      if(orbit.numMask == 0) {	
	if(loop == 0) {
	  SD_asym[func][count_asym]      = sd_asym;
	  SD_tr_asym[func][count_asym++] = sd_tr_asym;
	}  
	else if(loop == 1) {
	  SD_sym[func][count_sym++] = sd_asym;
	}
      }
      else {
	for(num = 0; num < orbit.numMask; num++) {
	  high = sd_asym & orbit.list[num];
	  phase = 0;
	  while(high) {
	    high &= high - 1;
	    phase++;
	  }
	  if(  (phase < orbit.lim[2*num])
	       ||(phase > orbit.lim[2*num + 1])) break;
	}
	if(num == orbit.numMask) {
	  if(loop == 0) {
	    SD_asym[func][count_asym]      = sd_asym;
	    SD_tr_asym[func][count_asym++] = sd_tr_asym;
	  }  
	  else if(loop == 1) {
	    SD_sym[func][count_sym++] = sd_asym;
	  }
	} 
      } // end  orbit.numMask > 0
    }     /* end parity and angular momentum test */

         /* new particle configuration */

    for(loop = 0; (loop < Info_id[func].part) && ((matr[loop]+1) >= matr[loop+1]); loop++);
    matr[loop]++;
    for(loop--; loop >= 0; loop--)  matr[loop] = loop;
    
  } while(matr[Info_id[func].part] == orbit.num);  /* end of particle loop */

    /* Check the number of identical asymmetric and symmetric SD */

    if((count_asym != Info_id[func].num[0]) || (count_sym != Info_id[func].num[1]))  {
      printf("\n\nError in calculating the identical particle  SDfor case = %d", func);
      printf("\nFrom number_of_identical_SD(): num_asym = %d num_sym = %d",
                                               Info_id[func].num[0],Info_id[func].num[1]);
      printf("\n From function identical_SD_tr_sym_evenJ()");
      printf("\ncount_asym = %d  count_sym = %d\n",count_asym, count_sym);
      exit(1);
    }
  free(matr);  /* release temporary memory */

} /* End: function identical_SD_tr_sym() */

  /*
   * The function                                           
   *          identical_SD_no_sym(func)                   
   * calculates and stores the complete set of slater determinants (SD) 
   * for identical particles with no time reversal symmtry both for 
   * initial(func = 0) and the final case(func = 1).
   */

void identical_SD_no_sym(int func, struct m_orbit orbit)
{
  register int      loop, m, parity, count,num, phase;
  register ULL       sd_val, high;
  int               *matr;

  matr = (int *) vector(sizeof(int), Info_id[func].part + 1); /* temporary memory */

      /* the lowest neutron configuration */

  for(loop = 0; loop < Info_id[func].part; loop++) matr[loop] = loop;
  matr[loop] = orbit.num;

        /* calculate the total number of SD.*/

  count  = 0;   /* initialization */

  do  {    /* run through all configurations */
        /* find m-value and parity */

    for(m = 0, parity = +1, loop = 0; loop < Info_id[func].part; loop++)  {
       m        += orbit.mbas[matr[loop]].m;
       parity   *= orbit.mbas[matr[loop]].par;
    }
    if((m == Info_com[func].MJ) && (parity == Info_com[func].par))  {

        /* Code the configuration */

      for(sd_val = ULL_ZERO, loop = 0; loop < Info_id[func].part; loop++)  {
         sd_val += ULL_ONE << matr[loop];
      }
      if(orbit.numMask == 0) { 
	SD_asym[func][count++] = sd_val;
      }
      else  {   // check for particle limitations
	for(num = 0; num < orbit.numMask; num++) {
	  high = sd_val & orbit.list[num];
	  phase = 0;
	  while(high) {
	    high &= high - 1;
	    phase++;
	  }
	  if(  (phase < orbit.lim[2*num])
	       ||(phase > orbit.lim[2*num + 1])) break;
	}
	if(num ==orbit.numMask)  SD_asym[func][count++] = sd_val;
      }
    }
         /* new particle configuration */

    for(loop = 0; (loop < Info_id[func].part) && ((matr[loop]+1) >= matr[loop+1]); loop++);
    matr[loop]++;
    for(loop--; loop >= 0; loop--)  matr[loop] = loop;
    
  } while(matr[Info_id[func].part] == orbit.num);  /* end running through all configurations */

    /* Check the number of identical SD */

  if(Info_id[func].num[0] != count)  {
    printf("\n\nError in calculating the identical particle SD for case = %d", func);
    printf("\nFrom number_of_identical_SD(): number = %d",Info_id[func].num[0]);
    printf("\n From function identical_SD_no_sym(): count = %d\n",count);
    exit(1);
  }
  
  free(matr);  /* release temporary memory */

} /* End: function identical_SD_no_sym() */

    /*
     * The function                                          
     *       void  id_particle_transitions()                
     * calculate for identical particles all one-particle
     * transition matrix elements  between eigenstates 
     * generated through a Lanczos procedure.
     */

void id_particle_transitions(int trans_type)
{
   register int     trans, loop, init_num, init_save, final_num, j_trans, dim;
   double           clebsch, value = D_ZERO;
   struct table_l   *table_id;
   double           xxx;
                /*
                 * If initial and final SD's are equal  calculate and 
                 * store all diagonal matrix elements between SD's 
                 */
   table_id = Info_com[0].type ? Table_pn_l[1] : Table_pn_l[0]; /* initialization */

   if(Model[0].same_basis == YES)   {
      id_diag_trans(Info_id[0].part, Info_id[0].num[0],                /* asym. or nosym.*/
                           table_id, SD_asym[0], SD_diag);
      if((Info_com[0].tr_sym) && (Info_com[0].J_even && Info_com[1].J_even))  {
         id_diag_trans(Info_id[0].part, Info_id[0].num[1],
                       table_id, SD_sym[0], (SD_diag + Info_id[0].num[0]));
      }
   } /* end calculation of diagonal matrix elements */

  j_trans = 0;       /* initialization */
   switch(trans_type)  {     /* Basic diagonal matrix elements and output heading */
       case  E2_TRANSITION: 
             j_trans = 4;  /* twice angular momentum transfer  */

             fprintf(Files.out_data,"\n\n            E2 transitions:");
             fprintf(Files.out_data,"\nInit(no 2J) Final(no 2J) (gr.st.no=0):");
             fprintf(Files.out_data," B(E2) in  (e**2) * (b**4)");
             break;
       case  M1_TRANSITION: 
             j_trans = 2;  /* twice angular momentum transfer  */
 
             fprintf(Files.out_data,"\n\n            M1 transitions:");
             fprintf(Files.out_data,"\nInit(no 2J) Final(no 2J) (gr.st.no=0):");
             fprintf(Files.out_data," B(M1) in unit (e*hbar/(2*M*c)");
             break;
       case  E3_TRANSITION: 
             j_trans = 6;  /* twice angular momentum transfer  */

             fprintf(Files.out_data,"\n\n            E3 transitions:");
             fprintf(Files.out_data,"\nInit(no 2J) Final(no 2J) (gr.st.no=0):");
             fprintf(Files.out_data,"  B(E3) in  (e**2) * (b**6)");
             break;
       default:
          printf("\n\nError in function one_particle_transitions():");
          printf("\n   transition type = %d is wrong!!!!!\n",trans_type);
          exit(1);
   } /* end switch() */

    /* Start the calculation of the transition matrix elements.*/

   for(trans = 0,init_save = -1; trans < Model[0].num_trans; trans ++) {

          /* Check angular momentum relation */

      clebsch = clebsch_gordan(Model[0].trans[trans].j_init,j_trans, 
                    Model[0].trans[trans].j_final, Model[0].MJ, Model[1].MJ - Model[0].MJ);

      if(fabs(clebsch) < ZERO_LIMIT)  {
         fprintf(Files.out_data,"\nTransition no %d can not be calculated!!!",trans);
         fprintf(Files.out_data,"\nTransition type: %s",Model[0].calc);
         fprintf(Files.out_data,"   2*J_init = %2d 2*J_final = %2d\n",
                               Model[0].trans[trans].j_init,Model[0].trans[trans].j_final);
         continue;
      }
             /* read the initial eigenvector into Eigen_vec[] */

      if((init_num = Model[0].trans[trans].num_init) != init_save) {

         read_id_eigenvector(File_eigen_init, init_num, Info_id[0].num[0],
                                                        Info_id[0].num[1], Eigen_vec); 

                   /* initialize the vector Trans_contr to ZERO. */

         dim = (Info_com[1].tr_sym && Info_com[1].J_even) 
                ? Info_id[1].num[0] + Info_id[1].num[1] : Info_id[1].num[0];

         for(loop = 0; loop < dim; loop++) Trans_contr[loop] = D_ZERO;

                   /* Calculate: |Trans_contr> = OP * |Eigen_vec> */

         switch(trans_type)  {                /* transition matrix elements */
             case  E2_TRANSITION: 
             case  M1_TRANSITION: 
                  id_trans_contribution(Info_id[0].part,Info_id[0].orb, Info_id[0].num[0], 
                                      Info_id[0].num[1], Info_id[1].num[0], Info_id[1].num[1],
                                                               table_id, Eigen_vec, Trans_contr);
                                  break;
             case  E3_TRANSITION: 
                  id_trans_contribution(Info_id[0].part,Info_id[0].orb, Info_id[0].num[0], 
                                      Info_id[0].num[1], Info_id[1].num[0], Info_id[1].num[1],
                                                               table_id, Eigen_vec, Trans_contr);

                                  break;
         } /* end switch() */
      }  /* end if(new initial eigenvector) */

      final_num = Model[0].trans[trans].num_final;

      read_id_eigenvector(File_eigen_final, final_num, Info_id[1].num[0],
                                                            Info_id[1].num[1], Eigen_vec); 

      xxx   = overlap(Info_id[1].num[0], Info_com[1].J_even*Info_id[1].num[1],
                        Trans_contr, Eigen_vec);
      value = xxx / clebsch;
   
      value *= value * (Model[0].trans[trans].j_final +1)/(Model[0].trans[trans].j_init +1);
      fprintf(Files.out_data,"\n    %2d %2d        %2d %2d            ", /* Print result */
          init_num, Model[0].trans[trans].j_init,final_num, Model[0].trans[trans].j_final);
      fprintf(Files.out_data,"               %14.8E",value);
      init_save = init_num; /* saved for next transition */
   } /* all transitions calculated */
} /* End: function id_particle_transitions() */

      /*
       * The function 
       *                id_diag_trans()
       * calculates for identical particles the diagonal transition matrix
       * elements <SD|OP|SD> of the one-particle operator OP (E2 or M1). 
       * The result is stored in a diag_ptr[]. The operator OP is defined
       * through its one-particle matrix elements found in table_ptr[].diag.
       */

void id_diag_trans(int part_num, int num_sd, struct table_l *table_ptr, 
                                                         ULL *sd_ptr, double *diag_ptr)
{
   register int         orb, part, loop;
   register ULL          pos;
   double               value;
 
   for(loop = 0; loop < num_sd; loop++, sd_ptr++) { /* run through all SD's */

      value  = 0.0;             /* initialization */
      pos    = ULL_ONE;
      orb    = 0;
      part   = 0;
      do   {            /* particle loop */
         for( ; !((*sd_ptr) & pos); orb++, pos <<= 1);
         value += table_ptr[orb++].diag;
         pos <<= 1;
      } while(++part < part_num);  /* end l particle loop */
      diag_ptr[loop]   = value;
   } /* end loop through all SD */
} /* End: function id_diag_trans() */

    /*
     * The function            
     *    void read_id_eigenvector()
     * reads the eigenvector no n generated through a Lanczos procedure
     * with num_asym asymmetric and num_sym symmetric elements. If no 
     * time-reversal symmetry num_sym = 0.
     * The function returns the amplitudes in vec_ptr[].
     */

void read_id_eigenvector(char *file_name,int eigen_num, int num_asym,
                                         int num_sym,  double *vec_ptr)
{
   int     n, totDim; 
   FILE    *file_ptr;

   printf("\n\n\nStart reading eigenvector number = %d\n\n",eigen_num);

    /* Open the file  for the  eigenvectors. */

   if((file_ptr = fopen(file_name,"rb")) ==  NULL_PTR)  {
      printf("\n\nError from function read_id_eigenvector:");
      printf("\nWrong eigenvector file_name = %s\n", file_name);
      exit(1);
   }
   totDim = num_asym + num_sym;
   fseek(file_ptr,(long)(eigen_num * (2 * sizeof(int) + totDim * sizeof(float))),SEEK_SET);

   if(fread((void *)&n,(size_t)sizeof(int),1,file_ptr) != 1) { /* read  vector number*/
      printf("\n\nError in function read_id_eigenvector():");
      printf("\nSomething is wrong in reading the vector number n = %d\n",n);
      exit(1);
   }
   if(eigen_num != n)  {
      printf("\n\nError in function read_id_eigenvector():");
      printf("\nWrong eigenvector number in file %s",file_name );
      printf("\n The correct eigenvector has not been found");
      printf("\neigenvector number = %d   number read = %d\n",eigen_num,n);
      exit(1);
   } /* end of if-test */

   if(fread((void *)&n,(size_t)sizeof(int),1,file_ptr) != 1) { /* read  dimension from file */
      printf("\n\nError in function read_id_eigenvector():");
      printf("\nSomething is wrong in reading total dimension = %d\n", n);
      printf("\n for eigenvector no %d\n", eigen_num);
      exit(1);
   }
   if(n != totDim)  {
      printf("\n\nError in function read_id_eigenvector():");
      printf("\nWrong dimension read for eigenvector no %d", n);
      printf("\n The correct eigenvector has not been found");
      printf("\nDimension should be %d -  number read = %d\n", totDim, n);
      exit(1);
   } /* end of if-test */

       /*
       ** read vector in float format 
       ** and convert to double
       */
   {
     int     k;
     float   *temp_mem, *init_ptr;
     double  *final_ptr;

     temp_mem   = (float *) malloc((ULL)(totDim * sizeof(float)));

     if(fread((void *) temp_mem, (size_t) sizeof(float),        /* read the vector */
	      (size_t) totDim, file_ptr) != (size_t) totDim) {
       printf("\n\n Error in function read_id_eigenvector():");
       printf("\nWrong number of elements in file %s",file_name );
       printf("\nvector_num = %d elements = %d\n",eigen_num, totDim);
       exit(1);
     }   /* end of if-test  */

     init_ptr  = temp_mem;
     final_ptr = vec_ptr;

     for(k = 0; k < totDim; k++)  {
       *(final_ptr++) = (double)(*(init_ptr++));
     }
     free(temp_mem);
   } 

   fclose(file_ptr);
   return;

} /* End: function read_id_eigenvector() */

     /*
      * The function
      *            id_trans_contribution()
      * operates for identical particles with a one-particle tensor
      * operator of M1, E2 or E3 type on a initial eigenfunction
      * stored in the vector init[] to produce a new vector final[]
      *               |final > = OP * |init >
      */

void  id_trans_contribution(int part_num, int orb_num, int init_asym_sd, int init_sym_sd, 
           int final_asym_sd, int final_sym_sd, struct table_l *table_id, double *init, double *final)
{
   if((Info_com[0].tr_sym) && (Info_com[1].tr_sym))   {

        /* initial and final eigenvectors  have time-reversal symmetry, */

      if(Model[0].same_basis == YES)   {   /* same basis  - diagonal contribution */

         diag_contr(init_asym_sd, 2, SD_diag, init, final);
         if(Info_com[0].J_even && Info_com[1].J_even)  {             
            diag_contr(init_sym_sd, 1, SD_diag + init_asym_sd,
                              init + init_asym_sd, final + init_asym_sd);
         }
      } /* end diagonal contributions */

                /* nondiagonal contribution */

      id_nondiag_asymI(part_num, orb_num, init_asym_sd,
                           final_asym_sd, final_sym_sd, table_id, init, final); 
      if(Info_com[0].J_even) { 
         id_nondiag_symI(part_num, orb_num, init_sym_sd,
                  final_asym_sd, final_sym_sd, table_id, init + init_asym_sd, final); 
      }
   } /* end time-reversal symmetry in both initial and final states. */

   else if(Info_com[0].tr_sym && !Info_com[1].tr_sym)   {

               /* 
                * The initial eigenvectors are calculated with time-reversal 
                * symmetry, but the final eigenvectors have no symmetry. 
                */

      id_nondiag_asymII(part_num, orb_num, init_asym_sd, final_asym_sd, table_id, init, final);
      if(Info_com[0].J_even)   { 
         id_nondiag_symII
            (part_num, orb_num, init_sym_sd, final_asym_sd, table_id, init + init_asym_sd, final);
      }   
   } /* end time-reversal symmetry in init states and no symmetry in final states */

   else if(!Info_com[0].tr_sym && Info_com[1].tr_sym)   {

               /* 
                * The final eigenvectors are calculated with time-reversal 
                * symmetry, but the initial eigenvectors have no symmetry. 
                */

      id_nondiag_nosymIII
         (part_num, orb_num, init_asym_sd, final_asym_sd, final_sym_sd, table_id, init, final);
   } /* end time-reversal symmetry in final states and no symmetry in init states */

   else if(!Info_com[0].tr_sym && !Info_com[1].tr_sym)   {

            /* Both initial and final eigenvectors are calculated without time-reversal.*/

      if(Model[0].same_basis == YES)   {   /* same basis  - diagonal contribution */

         diag_contr(init_asym_sd, 1, SD_diag, init, final);

      } /* end diagonal contributions */

                  /* non-diagonal contributions */

      id_nondiag_nosymIV(part_num, orb_num, init_asym_sd, final_asym_sd, table_id, init, final);

   } /* end  with no symmetry  in both initial and final states.*/

   else   { 
      printf("\n\nError from function m1_e2_e3transition()");
      printf("\nSomething is wrong for the time-reversal symmetry!!!!!");
      printf("Info_com[0].tr_sym = %d  Info_com[1].tr_sym = %d\n",
                                 Info_com[0].tr_sym,Info_com[1].tr_sym);
      exit(1);
   }

}  /* End: function id_trans_contribution() */

       /*
        * The function 
        *          diag_contr()
        * calculates and stores the diagonal contribution <SD|OP|SD> 
        * initial vector |init[]> to final vector |final[]>
        */

void diag_contr(int num_SD, int factor, double *sd_diag, double *init, double *final)
{
   if(num_SD == 0)  return;     
   do  {
      *(final++) = factor * (*(init++)) * (*(sd_diag++));
   } while(--num_SD);
   return;
} /* End: function diag_contr() */

  /*
   * The function                                 
   *        id_nondiag_asymI()         
   * calculates and stores for the identical particle case the contribution from the 
   * non-diagonal matrix elements <SD'|OP|SD> to the final vector from all the 
   * asym. and tr_asym.  terms in the init vector |init> which has time-reversal symmetry
   * into a final state with time-reversal symmetry.
   * The function includes the case with the same basis for initial and final states, 
   * i.e. even parity transitions as well as the case with different basis in initial
   * and final states, i.e. odd parity transitions.
   * In the first case we must have  pointer relations 
   *            SD_asym[1] = SD_asym[0] and SD_tr_asym[1] = SD_tr_asym[0].
   */

void id_nondiag_asymI(int part_num, int num_orb, int num_init_asym_sd, int num_final_asym_sd,
                   int num_final_sym_sd, struct table_l *table_id, double *init, double *final)
{
   register int          orb, part, loop, phase, search_num1, search_num2, search_num3, 
                         low, high, j_final_phase;
   register ULL           pos, l_sd, l_tr_sd, tr_one, new_sd, new_tr_sd, sd_init, sd_tr_init;
   double                *final_sym; 
   struct trans_nondiag  *kl_ptr;           /* pointer to kl-matrix element list */

   j_final_phase = Info_com[1].J_even ?  +1 : -1;  /* initialization */
   tr_one   = ULL_ONE <<(num_orb -1);
   final_sym = final + num_final_asym_sd;

      /* run through all asymmetric |SD> in the initial time-reversed eigenstate */

   for(loop = 0; loop < num_init_asym_sd; init++,loop++)  {

      if(fabs(*init) < ZERO_LIMIT) continue;

      sd_init    = SD_asym[0][loop];       /* initialization */
      sd_tr_init = SD_tr_asym[0][loop];

           /* run through all contributions for a init asym |SD> */

      for(part = 0, orb = 0, pos = ULL_ONE; part < part_num; part++, orb++, pos  <<= 1)  {
         for( ; !(sd_init & pos); orb++, pos <<= 1);               /* l part. found */
         if((kl_ptr = table_id[orb].nondiag) == NULL_PTR) continue;

         l_sd       = sd_init ^ pos;               /* remove the l-particle */
         l_tr_sd = sd_tr_init ^ (tr_one >> orb); 

         search_num1 = 0;           /* search position in final asym |SD[]> */
         search_num2 = 0;           /* search position in final sym |SD[]>  */
         do   {                /* for fixed l run through all k contributions */

                /* Pauli principle: k  position must be  allowed in |l_sd> */ 

            for( ;(l_sd & kl_ptr->one); kl_ptr++);
            if(!(kl_ptr->one)) break;  /* no more one-particle matrix elements */

            new_sd    = l_sd ^ kl_ptr->one;       /* contribution found, the new SD-state */

	    //   Mask test on |new_sd>
	    {
	      int  num, phase;
	      ULL  high;
	      struct m_orbit   orbit;

	      orbit = Model[0].Z ? Z_orbit : N_orbit;

	      if(orbit.numMask > 0) {
		for(num = 0; num < orbit.numMask; num++) {
		  high = new_sd & orbit.list[num];
		  phase = 0;
		  while(high) {
		    high &= high - 1;
		    phase++;
		  }
		  if(  (phase < orbit.lim[2*num])
		       ||(phase > orbit.lim[2*num + 1])) break;
		}
		if(num < orbit.numMask) continue;
	      }
	    }
	    // end Mask test on   |new_sd>

            new_tr_sd = l_tr_sd ^ kl_ptr->tr_one; /* the new time-reversed SD-state */

            low    = new_sd & kl_ptr->two; /* Calculate the kl  permutation phase. */
            phase = +1;
            while(low) {
               low &= low - 1;
               phase = -phase;
            }
            if(new_sd < new_tr_sd)  {      /* the new_sd in SD_asym.*/
               low  = search_num1;
               high = num_final_asym_sd;
               while(1)  {
                  search_num1 = (low + high)>> 1;
                  if(new_sd < SD_asym[1][search_num1])      high = search_num1 -1;
                  else if(new_sd > SD_asym[1][search_num1]) low  = search_num1 + 1;
                  else                                      break;
               }
                 /*
                  * Both  SD_asym[0][loop]    --> SD_asym[1][loop]
                  *       SD_tr_asym[0][loop] --> SD_tr_asym[1][loop]
                  * contributions are stored in final[search_num1].
                  */
                  
               final[search_num1] +=  2 * phase * kl_ptr->val * (*init);

            } /* end new config in final SD_asym[] */

            else  if((j_final_phase == +1) && (new_sd == new_tr_sd)) { /* new_sd in SD_sym.*/
               low  = search_num2;                       
               high = num_final_sym_sd;
               while(1)  {
                  search_num2 = (low + high)>> 1;
                  if(new_sd < SD_sym[1][search_num2])      high = search_num2 -1;
                  else if(new_sd > SD_sym[1][search_num2]) low  = search_num2 + 1;
                  else                                  break;
               }
                 /*
                  * Both  SD_asym[0][loop]    --> SD_sym[1][loop]
                  *       SD_tr_asym[0][loop] --> SD_sym[1][loop]
                  * contributions are stored in final_sym[search_num2].
                  */

               final_sym[search_num2]  +=  2 * phase * kl_ptr->val * (*init);

            } /* end new config in final SD_sym[] */

            else if(new_sd > new_tr_sd)  {

                    /*
                     * |new_sd> is in SD_tr_asym[] and |tr_new_sd> in SD_asym[]. 
                     * Thus a binary search in SD_asym[].
                     */

               low = 0;
               high = num_final_asym_sd;
               while(1)  {
                  search_num3 = (low + high) >> 1;
                  if(new_tr_sd < SD_asym[1][search_num3])      high = search_num3 - 1;
                  else if(new_tr_sd > SD_asym[1][search_num3]) low  = search_num3 + 1;
                  else                                         break;
               }
                 /*
                  * Both  SD_asym[0][loop]    --> SD_sym[1][loop]
                  *       SD_tr_asym[0][loop] --> SD_sym[1][loop]
                  * contributions are stored in final[search_num3].
                  */

               final[search_num3] +=  2 * phase * kl_ptr->val * j_final_phase * (*init);

            } /* end new config in final SD_tr_asym[] */

         } while((++kl_ptr)->one); /* end do-loop for contribution from all matr. elem */
      } /* end l particle loop */

   } /* end of for-loop through the SD components */

} /* End: function id_nondiag_asymI() */

  /*
   * The function                                 
   *        id_nondiag_symI()         
   * calculates and stores for the identical particle case the contribution from the 
   * non-diagonal matrix elements <SD'|OP|SD> to the final vector from all the 
   * sym.  terms in the init vector |init> which has time-reversal symmetry
   * into a final state with time-reversal symmetry.
   * The function includes the case with the same basis for initial and final states, 
   * i.e. even parity transitions as well as the case with different basis in initial
   * and final states, i.e. odd parity transitions.
   * In the first case we must have  pointer relations 
   *            SD_asym[1] = SD_asym[0] and SD_tr_asym[1] = SD_tr_asym[0].
   * Futhermore, both J_initial and (L + J_final) must be even. Otherwise, the function
   * should not be called.
   */

void id_nondiag_symI(int part_num, int num_orb, int num_init_sym_sd, int num_final_asym_sd,
                   int num_final_sym_sd, struct table_l *table_id, double *init, double *final)
{
   register int          orb, part, loop, phase, search_num1, search_num2, 
                         low, high, j_final_phase;
   register ULL           pos, l_sd, l_tr_sd, tr_one, new_sd, new_tr_sd, sd_init;
   double                *final_sym; 
   struct trans_nondiag  *kl_ptr;           /* pointer to kl-matrix element list */

   j_final_phase = Info_com[1].J_even ?  +1 : -1;   /* initializations */
   tr_one   = ULL_ONE <<(num_orb -1);
   final_sym = final + num_final_asym_sd;

      /* run through all asymmetric |SD> in the initial time-reversed eigenstate */

   for(loop = 0; loop < num_init_sym_sd; init++,loop++)  {

      if(fabs(*init) < ZERO_LIMIT) continue;

      sd_init    = SD_sym[0][loop];       /* initialization */

           /* run through all contributions for a init asym |SD> */

      for(part = 0, orb = 0, pos = ULL_ONE; part < part_num; part++, orb++, pos  <<= 1)  {
         for( ; !(sd_init & pos); orb++, pos <<= 1);               /* l part. found */
         if((kl_ptr = table_id[orb].nondiag) == NULL_PTR) continue;

         l_sd    = sd_init ^ pos;               /* remove the l-particle */
         l_tr_sd = sd_init ^ (tr_one >> orb); 

         search_num1 = 0;           /* search position in final asym |SD[]> */
         search_num2 = 0;           /* search position in final sym |SD[]>  */
         do   {                /* for fixed l run through all k contributions */

                /* Pauli principle: k  position must be  allowed in |l_sd> */ 

            for( ;(l_sd & kl_ptr->one); kl_ptr++);
            if(!(kl_ptr->one)) break;  /* no more one-particle matrix elements */

            new_sd = l_sd ^ kl_ptr->one;          /* contribution found, the new SD-state */


	    //   Mask test on |new_sd>
	    {
	      int  num, phase;
	      ULL  high;
	      struct m_orbit   orbit;

	      orbit = Model[0].Z ? Z_orbit : N_orbit;

	      if(orbit.numMask > 0) {
		for(num = 0; num < orbit.numMask; num++) {
		  high = new_sd & orbit.list[num];
		  phase = 0;
		  while(high) {
		    high &= high - 1;
		    phase++;
		  }
		  if(  (phase < orbit.lim[2*num])
		       ||(phase > orbit.lim[2*num + 1])) break;
		}
		if(num < orbit.numMask) continue;
	      }
	    }
	    // end Mask test on   |new_sd>

            new_tr_sd = l_tr_sd ^ kl_ptr->tr_one; /* the new time-reversed SD-state */

            low    = new_sd & kl_ptr->two; /* Calculate the kl  permutation phase. */
            phase = +1;
            while(low) {
               low &= low - 1;
               phase = -phase;
            }
            if(new_sd < new_tr_sd)  {      /* the new_sd in SD_asym.*/
               low  = search_num1;
               high = num_final_asym_sd;
               while(1)  {
                  search_num1 = (low + high)>> 1;
                  if(new_sd < SD_asym[1][search_num1])      high = search_num1 -1;
                  else if(new_sd > SD_asym[1][search_num1]) low  = search_num1 + 1;
                  else                                      break;
               }
                 /*
                  * Both  SD_sym[0][loop] --> SD_asym[1][loop]
                  *       SD_sym[0][loop] --> SD_tr_asym[1][loop]
                  * contributions are stored in final[search_num1].
                  */

               final[search_num1] +=  2 * phase * kl_ptr->val * (*init);

            } /* end new config in final SD_asym[] */

            else  if((j_final_phase == +1) && (new_sd == new_tr_sd)) { /* new_sd in SD_sym.*/
               low  = search_num2;                       
               high = num_final_sym_sd;
               while(1)  {
                  search_num2 = (low + high)>> 1;
                  if(new_sd < SD_sym[1][search_num2])      high = search_num2 -1;
                  else if(new_sd > SD_sym[1][search_num2]) low  = search_num2 + 1;
                  else                                  break;
               }
                 /*  SD_sym[0][loop] --> SD_sym[1][loop] contribution */

               final_sym[search_num2]  += phase * kl_ptr->val * (*init);

            } /* end new config in final SD_sym[] */
        } while((++kl_ptr)->one); /* end do-loop for contribution from all matr. elem */
      } /* end l particle loop */

   } /* end of for-loop through the SD components */

} /* End: function id_nondiag_symI() */

  /*
   * The function                                 
   *        id_nondiag_asymII()         
   * calculates and stores for the identical particle case the contribution from the 
   * non-diagonal matrix elements <SD'|OP|SD> to the final vector from all the 
   * asym. and tr_asym. terms for the case of time-reversal symmetry in the init vector |init>
   * into a final state with no time-reversal symmetry.
   */

void id_nondiag_asymII(int part_num, int num_orb, int num_asym_sd,
                            int num_final_sd, struct table_l *table_id, double *init, double *final)
{
   register int          orb, part, loop, phase, search_num, low, high, j_init_phase;
   register ULL           pos, tr_one, l_sd, new_sd, sd_init, sd_tr_init;
   struct trans_nondiag  *kl_ptr;           /* pointer to kl-matrix element list */

   tr_one       = ULL_ONE <<(Info_id[0].orb -1);
   j_init_phase = Info_com[0].J_even ? +1 : -1; /* = +1(-1) if init state has even (odd) J_value */

      /* run through all asymmetric |SD> in the initial time-reversed eigenstate */

   for(loop = 0; loop < num_asym_sd; init++,loop++)  {

      if(fabs(*init) < ZERO_LIMIT) continue;

      sd_init    = SD_asym[0][loop];       /* initialization */
      sd_tr_init = SD_tr_asym[0][loop];

           /* run through all contributions for a init asym |SD> */

      for(part = 0, orb = 0, pos = ULL_ONE; part < part_num; part++, orb++, pos  <<= 1)  {
         for( ; !(sd_init & pos); orb++, pos <<= 1);               /* l part. found */
         if((kl_ptr = table_id[orb].nondiag) != NULL_PTR)  {
            l_sd       = sd_init ^ pos;               /* remove the l-particle */
            search_num = 0;           /* search position in final |SD[]> */
            do   {                /* for fixed l run through all k contributions */

                /* Pauli principle: k  position must be  allowed in |l_sd> */ 

               for( ;(l_sd & kl_ptr->one); kl_ptr++);
               if(!(kl_ptr->one)) break;  /* no more one-particle matrix elements */

               new_sd = l_sd ^ kl_ptr->one; /* contribution found, the new SD-state */

	    //   Mask test on |new_sd>
	    {
	      int  num, phase;
	      ULL  high;
	      struct m_orbit   orbit;

	      orbit = Model[0].Z ? Z_orbit : N_orbit;

	      if(orbit.numMask > 0) {
		for(num = 0; num < orbit.numMask; num++) {
		  high = new_sd & orbit.list[num];
		  phase = 0;
		  while(high) {
		    high &= high - 1;
		    phase++;
		  }
		  if(  (phase < orbit.lim[2*num])
		       ||(phase > orbit.lim[2*num + 1])) break;
		}
		if(num < orbit.numMask) continue;
	      }
	    }
	    // end Mask test on   |new_sd>

               low    = new_sd & kl_ptr->two; /* Calculate the kl  permutation phase. */
               phase = +1;
               while(low) {
                  low &= low - 1;
                  phase = -phase;
               }
               low  = search_num; /* search for the |new_sd>  in no-sym. |SD_asym[1]> */
               high = num_final_sd;
               while(1)  {
                  search_num = (low + high)>> 1;
                  if(new_sd < SD_asym[1][search_num])      high = search_num -1;
                  else if(new_sd > SD_asym[1][search_num]) low  = search_num + 1;
                  else                                      break;
               }

                  /* contribution SD_asym[0][loop] --> SD_nonsym[1][search_num]  */

               final[search_num] += phase * kl_ptr->val * (*init);

             } while((++kl_ptr)->one); /* end do-loop for contr. from all ij_pairs */

         } /* end new config in non-sym |SD_asym[1]> from |SD_asym[0]> */

                /* contributions from the time-reversed asym. |SD_asym[0]> */

         if((kl_ptr = table_id[num_orb - 1 - orb].nondiag) != NULL_PTR)  {

            l_sd = sd_tr_init ^ (tr_one >> orb);   /* remove the time-reversed l-particle */

            search_num = 0;           /* search position in final SD[] */
            do   {                /* for fixed l run through all k contributions */

                /* Pauli principle: k  position must be  allowed in |l_sd> */ 

               for( ;(l_sd & kl_ptr->one); kl_ptr++);
               if(!(kl_ptr->one)) break;  /* no more one-particle matrix elements */

                new_sd = l_sd ^ kl_ptr->one;    /* the new time-reversed SD-state */
              
               low = new_sd & kl_ptr->two; /* Calculate the kl  permutation phase. */
               phase = +1;
               while(low) {
                  low &= low - 1;
                  phase = -phase;
               }
               low  = search_num; /* search for the |new_sd>  in no-sym. |SD_asym[1]> */
               high = num_final_sd;
               while(1)  {
                  search_num = (low + high)>> 1;
                  if(new_sd < SD_asym[1][search_num])      high = search_num -1;
                  else if(new_sd > SD_asym[1][search_num]) low  = search_num + 1;
                  else                                     break;
               }
                    /* contribution SD_tr_asym[0][loop] --> SD_nonsym[1][search_num]  */

               final[search_num] += phase * kl_ptr->val * j_init_phase * (*init);

            } while((++kl_ptr)->one); /* end do-loop for contr. from all ij_pairs */

         } /* end new config in final non_sym |SD[1]> from |SD_tr_asym[0]> */

      } /* end l particle loop */

   } /* end of for-loop through the asym. SD components */

} /* End: function id_nondiag_asymII() */

  /*
   * The function                                 
   *        id_nondiag_symII()         
   * calculates and stores for the identical particle case the contribution from the 
   * non-diagonal matrix elements <SD'|OP|SD> to the final vector from all the 
   * sym. terms for the case of time-reversal symmetry in the init vector |init_SD_sym>
   * into a final state with no time-reversal symmetry.
   * Note that init[] starts with the first symmetric amplitudes.
   */

void id_nondiag_symII(int part_num, int num_orb, int num_sym_sd, 
                      int num_final_sd, struct table_l *table_id, double *init, double *final)
{
   register int          orb, part, loop, phase, search_num, low, high;
   register ULL           pos, l_sd, new_sd, sd_init;
   struct trans_nondiag  *kl_ptr;           /* pointer to kl-matrix element list */

       /* run through all symmetric |SD> in the initial time-reversed eigenstate */

   for(loop = 0; loop < num_sym_sd; init++, loop++)  {

      if(fabs(*init) < ZERO_LIMIT) continue;

      sd_init = SD_sym[0][loop];

            /* run through all contributions for the given |SD> */

      for(part = 0, orb = 0, pos = ULL_ONE; part < part_num; part++, orb++, pos  <<= 1)  {
         for( ; !(sd_init & pos); orb++, pos <<= 1);               /* l part. found */
         if((kl_ptr = table_id[orb].nondiag) == NULL_PTR) continue;

         l_sd    = sd_init ^ pos;               /* remove the l-particle */
         search_num = 0;           /* search position in final |SD[]> */
         do   {                /* for fixed l run through all k contributions */

                /* Pauli principle: k  position must be  allowed in |l_sd> */ 

            for( ;(l_sd & kl_ptr->one); kl_ptr++);
            if(!(kl_ptr->one)) break;  /* no more one-particle matrix elements */

            new_sd = l_sd ^ kl_ptr->one; /* contribution found, the new SD-state */

	    //   Mask test on |new_sd>
	    {
	      int  num, phase;
	      ULL  high;
	      struct m_orbit   orbit;

	      orbit = Model[0].Z ? Z_orbit : N_orbit;

	      if(orbit.numMask > 0) {
		for(num = 0; num < orbit.numMask; num++) {
		  high = new_sd & orbit.list[num];
		  phase = 0;
		  while(high) {
		    high &= high - 1;
		    phase++;
		  }
		  if(  (phase < orbit.lim[2*num])
		       ||(phase > orbit.lim[2*num + 1])) break;
		}
		if(num < orbit.numMask) continue;
	      }
	    }
	    // end Mask test on   |new_sd>

            low    = new_sd & kl_ptr->two; /* Calculate the kl  permutation phase. */
            phase = +1;
            while(low) {
               low &= low - 1;
               phase = -phase;
            }
            low  = search_num; /* search for the |new_sd>  in no-sym. |SD_asym[1]> */
            high = num_final_sd;
            while(1)  {
               search_num = (low + high)>> 1;
               if(new_sd < SD_asym[1][search_num])      high = search_num - 1;
               else if(new_sd > SD_asym[1][search_num]) low  = search_num + 1;
               else                                     break;
            }
              /* contribution SD_sym[0][loop] --> SD_nosym[1][search_num] */

            final[search_num] += phase * kl_ptr->val * (*init);

         } while((++kl_ptr)->one); /* end do-loop for contribution from all ij_pairs */

      } /* end l particle loop */

   } /* end for-loop through all sym |SD> */

} /* End: function id_nondiag_symII() */

       /*
	* The function                                 
	*        id_nondiag_nosymIII()         
   	* calculates and stores for the identical particle case the contribution from
   	* the non-diagonal matrix elements <SD'|OP|SD> to the final vector for no
   	* time-reversal symmetry in initial into final states with time-reversal
        * symmetry.
   	*/

void id_nondiag_nosymIII(int part_num, int num_orb, int num_init_sd, int num_final_asym_sd,
                       int num_final_sym_sd, struct table_l *table_id, double *init, double *final)
{
   register int          orb, part, loop, phase, search_num1, search_num2, 
                         search_num3, high, j_final_phase;
   register ULL           low, pos, pos_tr, l_sd, new_sd, new_tr_sd,sd_init;
   double                *final_sym;
   struct trans_nondiag  *kl_ptr;            /* pointer to k-matrix element list    */

   j_final_phase = Info_com[1].J_even ? +1 : -1;     /* initialization */
   final_sym     = final + num_final_asym_sd;

                /* run through all initial nonsym. |SD> */

   for(loop = 0; loop < num_init_sd; init++, loop++)  {

      if(fabs(*init) < ZERO_LIMIT) continue;

      sd_init = SD_asym[0][loop];
     
      for(part = 0, orb = 0, pos = ULL_ONE; part < part_num; part++, orb++, pos <<= 1)   {
        for( ; !(sd_init & pos); orb++, pos <<= 1);   /* l loop */
        if((kl_ptr = table_id[orb].nondiag) == NULL_PTR) continue;

        l_sd    = sd_init ^ pos;                       /* l part. found */

        search_num1 = 0;                        /* search number for final |SD_asym> */
        search_num2 = 0;                        /* search number for final |SD_sym>  */ 

        do   {                        /* run through all matrix elements */

                /* Pauli principle: index k must be  allowed in |l_sd> */ 

           for( ;(l_sd & kl_ptr->one); kl_ptr++);
           if(!(kl_ptr->one)) break;  /* no more matrix elements */

           new_sd = l_sd ^ kl_ptr->one; /* contribution found, the new SD-state */

	    //   Mask test on |new_sd>
	    {
	      int  num, phase;
	      ULL  high;
	      struct m_orbit   orbit;

	      orbit = Model[0].Z ? Z_orbit : N_orbit;

	      if(orbit.numMask > 0) {
		for(num = 0; num < orbit.numMask; num++) {
		  high = new_sd & orbit.list[num];
		  phase = 0;
		  while(high) {
		    high &= high - 1;
		    phase++;
		  }
		  if(  (phase < orbit.lim[2*num])
		       ||(phase > orbit.lim[2*num + 1])) break;
		}
		if(num < orbit.numMask) continue;
	      }
	    }
	    // end Mask test on   |new_sd>

           low    = new_sd & kl_ptr->two; /* Calculate the kl  permutation phase. */
           phase = +1;
           while(low) {
              low &= low - 1;
              phase = -phase;
           }
                  /* calculate |new_tr_sd> - the time-reversed of |new_sd> */
  
           new_tr_sd = ULL_ZERO;
           pos_tr    = ULL_ONE;
           high      = num_orb - 1; /* orb num for a time-reversed particle */
           low       = part_num;
           do   {
              for( ; !(new_sd & pos_tr); high--, pos_tr <<= 1);
              new_tr_sd += ULL_ONE << high;
              pos_tr <<= 1;
              high--;
           } while(--low);   /* time-reversed |new_tr_sd> of |new_sd> calculated */ 

           if(new_sd < new_tr_sd)  {      /* the new_sd in SD_asym.*/
              low  = search_num1;
              high = num_final_asym_sd;
              while(1)  {
                 search_num1 = (low + high)>> 1;
                 if(new_sd < SD_asym[1][search_num1])      high = search_num1 -1;
                 else if(new_sd > SD_asym[1][search_num1]) low  = search_num1 + 1;
                 else                                      break;
              }
               /* contribution SD_non_sym[0][loop] --> SD_asym[1][search_num1] */

              final[search_num1] += phase * kl_ptr->val * (*init);

           } /* end new config in SD_asym[] */

           else  if((j_final_phase == +1) && (new_sd == new_tr_sd)) { /* new_sd in SD_sym.*/
              low  = search_num2;                       
              high = num_final_sym_sd;
              while(1)  {
                 search_num2 = (low + high)>> 1;
                 if(new_sd < SD_sym[1][search_num2])      high = search_num2 -1;
                 else if(new_sd > SD_sym[1][search_num2]) low  = search_num2 + 1;
                 else                                  break;
              }   
                 /* contribution SD_non_sym[0][loop] --> SD_sym[1][search_num2] */

              final_sym[search_num2]  += phase * kl_ptr->val * (*init);

           } /* end new config in SD_sym[] */

           else if(new_sd > new_tr_sd)  {

                    /*
                     * |new_sd> is in SD_tr_asym[] and |tr_new_sd> in SD_asym[]. 
                     * Thus a binary search in SD_asym[].
                     */

              low = 0;
              high = num_final_asym_sd;
              while(1)  {
                 search_num3 = (low + high) >> 1;
                 if(new_tr_sd < SD_asym[1][search_num3])      high = search_num3 - 1;
                 else if(new_tr_sd > SD_asym[1][search_num3]) low  = search_num3 + 1;
                 else                                         break;
              }
                 /* 
                  * contribution SD_non_sym[0][loop] --> SD_sym[1][search_num2] 
                  * stored in final[search_num3] 
                  */

              final[search_num3] +=  phase * kl_ptr->val * j_final_phase * (*init);

           } /* end new config in SD_tr_asym[] */

        } while((++kl_ptr)->one); /* end do-loop for contribution from all matr. elem */
     } /* end l particle loop */

  } /* end of for-loop through the SD components */

} /* End: function id_nondiag_nosymIII() */

       /*
	* The function                                 
	*        id_nondiag_nosymIV()         
   	* calculates and stores for the identical particle case the contribution from
   	* the non-diagonal matrix elements <SD'|OP|SD> to the final vector when both
        * the initial and final states have  no time-reversal symmetry.
   	*/

void id_nondiag_nosymIV(int part_num, int num_orb, int num_init_sd, int num_final_sd,
                                            struct table_l *table_id,double *init, double *final)
{
   register int          orb, part, loop, phase, search_num, high;
   register ULL           low, pos, l_sd, new_sd, sd_init;
   double                in_ampl;           /* current amplitude                   */
   struct trans_nondiag  *kl_ptr;            /* pointer to k-matrix element list    */

      /* run through all components in the asymmetric group */

   for(loop = 0; loop < num_init_sd; loop++)  {

      if(fabs(in_ampl = init[loop]) < ZERO_LIMIT) continue;

      sd_init    = SD_asym[0][loop];       /* initialization */
     
      for(part = 0, orb = 0, pos = ULL_ONE; part < part_num; part++, orb++, pos <<= 1)   {
        for( ; !(sd_init & pos); orb++, pos <<= 1);   /* l loop */
        if((kl_ptr = table_id[orb].nondiag) == NULL_PTR) continue;
        l_sd    = sd_init ^ pos;                       /* l part. found */
        search_num = 0;                        /* search in SD_asym[] */
        do   {                        /* run through all matrix elements */

                /* Pauli principle: index k must be  allowed in |l_sd> */ 

           for( ;(l_sd & kl_ptr->one); kl_ptr++);
           if(!(kl_ptr->one)) break;  /* no more matrix elements */

           new_sd = l_sd ^ kl_ptr->one; /* contribution found, the new SD-state */

	    //   Mask test on |new_sd>
	    {
	      int  num, phase;
	      ULL  high;
	      struct m_orbit   orbit;

	      orbit = Model[0].Z ? Z_orbit : N_orbit;

	      if(orbit.numMask > 0) {
		for(num = 0; num < orbit.numMask; num++) {
		  high = new_sd & orbit.list[num];
		  phase = 0;
		  while(high) {
		    high &= high - 1;
		    phase++;
		  }
		  if(  (phase < orbit.lim[2*num])
		       ||(phase > orbit.lim[2*num + 1])) break;
		}
		if(num < orbit.numMask) continue;
	      }
	    }
	    // end Mask test on   |new_sd>

           low    = new_sd & kl_ptr->two; /* Calculate the kl  permutation phase. */
           phase = +1;
           while(low) {
              low &= low - 1;
              phase = -phase;
           }
                   /* binary search */

           low  = search_num;
           high = num_final_sd;
           while(1)  {
              search_num = (low + high)>> 1;
              if(new_sd < SD_asym[1][search_num])      high = search_num - 1;
              else if(new_sd > SD_asym[1][search_num]) low  = search_num + 1;
              else                                   break;
           }
              /* contribution SD_nosym[0][loop] --> SD_nosym[1][search_num] */

           final[search_num] += phase * kl_ptr->val * in_ampl;

        } while((++kl_ptr)->one); /* end do-loop for contribution from all matr. elem */
     } /* end l particle loop */

  } /* end of for-loop through the SD components */

} /* End: function id_nondiag_nosymIV() */

  /*
   * The function 
   *             overlap()
   * calculates overlap between  |trans_contr> = OP * |init_vec>
   * and a given final eigenvector |final_eigen>
   */

double overlap(int num_asym_sd, int num_sym_sd, 
                       double *trans_contr, double *final_eigen)
{
   int         number;
   double      *eigen_ptr, result;

   result  = 0.0;    /* initialization */
   number    = num_asym_sd;
   if(number)   {
      eigen_ptr = final_eigen;
      do {
         result += (*(eigen_ptr++)) * (*(trans_contr++));
      } while(--number);

   } /* end contribution from the asymmetric/ no-symmetric terms. */

          /* 
           * Additional contribution from the symmetric components in the final
           * eigenvector if the final vector have time-reversal symmetry.
           */

   if( (number = num_sym_sd) ) {
      eigen_ptr = final_eigen + num_asym_sd;
      do {                /* contribution frn the sym terms */
         result += (*(eigen_ptr++)) * (*(trans_contr++));
      } while(--number);
   }  /* end symmetric contributions */

   return result;

} /* End: function overlap()  */
