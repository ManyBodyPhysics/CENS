      /*
      **  Program module
      **            pn-transition.c
      ** calculates different one-particle transition
      ** between eigenstates generated in a Lanczos procedure.
      ** for a system with both protons and neutrons
      ** The program module contains the following functions:
      **     1.    void pn_calculation()
      **     2.    void pn_slater_determinant()
      **     3.    void parity_check()
      **     4.    void number_of_groups()
      **     5.    void group_m_limit()
      **     6.    int number_of_proton_SD()
      **     7.    void proton_SD_configuration()
      **     8.    int number_of_neutron_SD()
      **     9.    void neutron_SD_configuration()
      **    10.    void pn_particle_transitions()
      **    11.    void group_change()
      **    12.    void pn_diag_matrix_elements()
      **    13.    void read_pn_eigenvector()
      **    14.    void pn_nondiag_iteration()
      **    15.    void nondiag_proton-contribution()
      **    16.    void nondiag_neutron_contribution()
      **    17.    double pn_overlap()
      */

#include "shell.h"
#include "sm.h"

void pn_calculation(void)
{

    /* Calculate and store the initial and final SD. */

   pn_slater_determinant();
 
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

      matrix_elements_in_m_scheme(0, E2_TRANSITION); /* protons */
      matrix_elements_in_m_scheme(1, E2_TRANSITION); /* neutrons */

         /* Allocate memory for the calculation of the transition matrix elements.*/

      memory_allocation(TRANSITION_MATRIX_ELEMENTS,0);

         /* Calculate the E2 transitions.*/

      pn_particle_transitions(E2_TRANSITION); 

   }
  
   if(!strncmp(Model[0].calc,"m1",2)) {

        /*
         * Allocate memory for the index-l table, the diagonal and off-diagonal 
         * single-particle  matrix elements.
         */

      memory_allocation(M1_E2_E3_MATRIX_ELEMENTS,0); 

      matrix_elements_in_m_scheme(0, M1_TRANSITION); /* protons */
      matrix_elements_in_m_scheme(1, M1_TRANSITION); /* neutrons */


         /* Allocate memory for the calculation of the transition matrix elements.*/

      memory_allocation(TRANSITION_MATRIX_ELEMENTS,0);

         /* Calculate the M1 transitions.*/

      pn_particle_transitions(M1_TRANSITION); 

   }
   if(!strncmp(Model[0].calc,"e3",2))  {

        /*
         * Allocate memory for the index-l table, the diagonal and off-diagonal 
         * single-particle  matrix elements for both protons and neutrons. 
         */

      memory_allocation(M1_E2_E3_MATRIX_ELEMENTS,0); 

      matrix_elements_in_m_scheme(0, E2_TRANSITION); /* protons */
      matrix_elements_in_m_scheme(1, E2_TRANSITION); /* neutrons */


         /* Allocate memory for the calculation of the transition matrix elements.*/

      memory_allocation(TRANSITION_MATRIX_ELEMENTS,0);

         /* Calculate the E3 transitions.*/

      pn_particle_transitions(E3_TRANSITION); 
    }

} /* End: function pn_calculation() */

  /*
   * The function                                                      
   *            void pn_slater_determinant() 
   * allocates memory, calculates and store the complete set of slater 
   * determinants (SD) for a system of Z protons and N neutrons both 
   * for the initial(func = 0) and the final case(func = 1).       
   */

void pn_slater_determinant(void)
{
   register int       func, limit, group;

   for(func = 0; func < 2; func++)   {

      parity_check(func);            /* parity for proton/neutron SD */
      number_of_groups(func);        /*  number of proton-neutron SD groups.*/ 

        /* memory to store details of the proton and neutron groups. */

      memory_allocation(PROTON_NEUTRON_GROUPS,func);

         /* Calculate the complete set of SD */

      Info_pn[func].Z_num = number_of_proton_SD(func);    /* protons */ 
      memory_allocation(PROTON_CONFIGURATION,func);
      proton_SD_configuration(func);

      Info_pn[func].N_num = number_of_neutron_SD(func);         /* neutrons */  
      memory_allocation(NEUTRON_CONFIGURATION,func);
      neutron_SD_configuration(func);          

       /* Calculate the total number of basic proton - neutron SD configurations. */ 

      limit = Info_pn[func].num_gr; 
      for(Info_com[func].num = 0, group = 0; group < limit; group++)  {
         Group[func][group].ampl_start 
             = Info_com[func].num; /* start pos. in lanczo vectors */ 
         Info_com[func].num  += Group[func][group].Z_numSD * Group[func][group].N_numSD;
      }
                   /* Basic output data */

     fprintf(Files.out_data,"\n\nNumber of proton SD configurations = %d", 
                                                               Info_pn[func].Z_num);
     fprintf(Files.out_data,"  Number of neutron SD configurations = %d", 
                                                               Info_pn[func].N_num);
     fprintf(Files.out_data,"\nTotal number of proton-neutron SD configurations = %d",
                                                                   Info_com[func].num);

   } /* end run through initial and final states, */

} /* End: function pn_slater_determinant() */

  /*
   * The function                                      
   *               parity_check()                      
   * calculates for the initial case (func = 0) and the final case ( func = 1) 
   * the parities for the proton and neutron separately and store the information
   * Info_pn[func].Z_par (proton particles) and Info_pn[func].N_par (neutron particles)
   * The information is stored as 
   *                      = -1    negative parity only                
   *                      = +1    positive parity only               
   *                      =  0    both parities                      
   */

void parity_check(int func)
{
  register int    loop, loop_max;
 
  Info_pn[func].Z_par = Z_orbit.mbas[0].par; /* proton parities */
  loop_max      = Info_pn[func].Z_orb / 2;
  for(loop = 1; (loop < loop_max) && (Z_orbit.mbas[loop].par == Info_pn[func].Z_par); loop++);
  if(loop < loop_max)  Info_pn[func].Z_par = 0; /* both parities */

    /* possible modification of parity due to proton number */

  if((Info_pn[func].Z_par == -1) && (!PARITY(Info_pn[func].Z_part))) Info_pn[func].Z_par = +1;

  Info_pn[func].N_par = N_orbit.mbas[0].par; /* neutron parities */
  loop_max      = Info_pn[func].N_orb / 2;
  for(loop = 1; (loop < loop_max) && (N_orbit.mbas[loop].par == Info_pn[func].N_par); loop++);
  if( loop < loop_max) Info_pn[func].N_par = 0; /* both parities */

    /* possible modification of parity due to neutron number */

  if((Info_pn[func].N_par == -1) && (!PARITY(Info_pn[func].N_part))) Info_pn[func].N_par = +1;

  if((Info_pn[func].Z_par * Info_pn[func].N_par) == -Info_com[func].par)  { /* Error check.*/
     printf("\n No proton-neutron SD is possible due to parity for case %d.",func);
     printf("\n From function parity_check() with total parity +\n");
     exit(1);
  }
  if((Info_pn[func].Z_par * Info_pn[func].N_par) == Info_com[func].par)  
     return;        /* no parity reduction.*/

  if((Info_pn[func].Z_par == 0) && (Info_pn[func].N_par != 0))  /* possible parity reduction.*/
       Info_pn[func].Z_par =  Info_com[func].par * Info_pn[func].N_par;
  else if((Info_pn[func].N_par == 0) && (Info_pn[func].Z_par != 0))
       Info_pn[func].N_par =  Info_com[func].par * Info_pn[func].Z_par;

} /* End: function parity_check() */

  /*
   * The function
   *        number_of_groups()
   * calculates the maximum m--values for protons and neutrons in the groups for initial 
   * case (func = 0) or final case (func = 1). Then it finds the total number of groups. 
   * It is assumed that both protons and neutrons are present i.e. Info_com[func].type > 1.
   */

void number_of_groups(int func)
{
  register int  loop_z, loop_n;

    /* Calculate maximum m-values for protons */

  group_m_limit(0, Info_pn[func].Z_par, Info_pn[func].Z_part, Z_orbit, Info_pn[func].Zm_max);

   /* Calculate maximum m-values for neutrons */

  group_m_limit(1, Info_pn[func].N_par, Info_pn[func].N_part, N_orbit, Info_pn[func].Nm_max);

    /* calculate max. and min. m-values and number of groups */

  Info_pn[func].num_gr = 0;    /* initialization */

  loop_z  = (Info_pn[func].Z_par == 0) ? 1: 0;
  do   {
     loop_n = ((Info_com[func].par == -1) && (Info_pn[func].Z_par == 0)) ? 1 - loop_z : loop_z;
     if(Info_pn[func].Zm_max[loop_z] <= Info_pn[func].Nm_max[loop_n])  {    
        Info_pn[func].Nm_min[loop_n] = Info_com[func].MJ - Info_pn[func].Zm_max[loop_z];
        Info_pn[func].Zm_min[loop_z] 
            = MAX(-Info_pn[func].Zm_max[loop_z],  Info_com[func].MJ 
                                                - Info_pn[func].Nm_max[loop_n]);
        Info_pn[func].Nm_max[loop_n] = Info_com[func].MJ - Info_pn[func].Zm_min[loop_z];
     }  
     else  {
        Info_pn[func].Zm_min[loop_z] = Info_com[func].MJ - Info_pn[func].Nm_max[loop_n];
        Info_pn[func].Nm_min[loop_n] 
           = MAX(-Info_pn[func].Nm_max[loop_n],  Info_com[func].MJ 
                                               - Info_pn[func].Zm_max[loop_z]);
        Info_pn[func].Zm_max[loop_z] = Info_com[func].MJ - Info_pn[func].Nm_min[loop_n];
     }  
     Info_pn[func].num_gr 
       = MAX(Info_pn[func].num_gr,  Info_pn[func].Zm_max[loop_z] 
                                  - Info_pn[func].Zm_min[loop_z] + 2);

  } while((--loop_z) >= 0);  /* end do - loop through the two parity cases.*/

} /* End: function number_of_groups() */

  /*
   * The function
   *        group_m_limit()
   * calculates the maximum m--values for protons and neutrons in the groups.
   */

void group_m_limit(int type, int parity, int part_num, struct m_orbit orbit, int *result)
{
  register int loop_1, loop_2,  m1, m2 = 0, par1, par2, number;

  if(abs(parity) > 1)  {
     printf("\nError message: Wrong input parity in function group_m_limit()");
     printf("\nInput parity =  %d\n",parity);
     exit(1);
   }

    /* possible particle - hole transformation */

  number = MIN(part_num, orbit.num - part_num);

    /* Calculate first maximum m-value */

  for(m1 = 0, par1 = +1, loop_1 = 0; loop_1 < number; loop_1++)  {
     m1   += orbit.mbas[loop_1].m;
     par1 *= orbit.mbas[loop_1].par;
  }
  if(abs(par1) != 1)  {
     printf("\nError: Data error from  function group_m_limit()");
     printf("\nParticle type = %d and  parity  = %d\n",type,par1);
     exit(1);
  }            
  result[0] = m1;
  result[1] = m1;
  if(parity == par1) return;  /* only one type of parity */
 
       /* remove m-value and parity for one particle */

  for(loop_2 = number - 1; loop_2 >= 0; loop_2--)  {
     m2   = m1 - orbit.mbas[loop_2].m;
     par2 = orbit.mbas[loop_2].par;
     for(loop_1 = number; loop_1 < orbit.num; loop_1++)  {
        if(orbit.mbas[loop_1].par != par2) break;
     }
     if(loop_1 < orbit.num) break;
  }

  if(loop_2 < 0)  { /* NB! never reach loop_2 < 0  */
     printf("\nError in function group_m_limit()");
     printf("\nSomething wrong with the single-particle orbits.");
     printf("\nOpposite parity not found for type %d",type);
     exit(1);
  }
  m2  += orbit.mbas[loop_1].m;
  par2 = -par1;

  if(parity == par2) {
    result[0] = m2;
    return;
  }
  switch(par2)  {
     case +1: result[0] = m2;
              break;
     case -1: result[1] = m2;
              break;
     default:  printf("\nData error in function group_m_limit()");
               printf("\nsecond parity - par2 = %d.",par2);
               exit(1);
  } /* end switch-loop */
    
} /* End: function group_m_limit() */

  /*
   * The function                                        
   *          number_of_proton_SD()                      
   * calculates and return the number of proton SD for the initial 
   * case (func = 0) or the final case (func = 1).
   * and store the result in pn->Z_num. 
   */

int number_of_proton_SD(int func)
{
   register int     loop, num, m, phase, parity;
   int              *matr;
   ULL              config, high;

   matr = (int *) vector(sizeof(int),Info_pn[func].Z_part + 1);  /* temporary memory */
   for(loop = 0; loop < Info_pn[func].Z_part; loop++) matr[loop] = loop; /* lowest config. */
   matr[Info_pn[func].Z_part] = Info_pn[func].Z_orb;
   num = 0;                   /* initialization */
   do  {   /* the proton loop */
 
     // test for Masking on current config
     

       // check for particle limitations
  	        // code the configuration

     for(config = ULL_ZERO, loop = 0; loop < Info_pn[func].Z_part; loop++)  {
       config += ULL_ONE << matr[loop];
     }
     for(loop = 0; loop < Z_orbit.numMask; loop++) {
       high = config & Z_orbit.list[loop];
       phase = 0;
       while(high) {
	 high &= high - 1;
	 phase++;
       }
       if(  (phase < Z_orbit.lim[2*loop])
	    ||(phase > Z_orbit.lim[2*loop + 1])) break;
     }
       // end  Masking on current config

     if(loop ==  Z_orbit.numMask) {
       for(parity = +1, m = 0, loop = 0; loop < Info_pn[func].Z_part; loop++)  {
         m      += Z_orbit.mbas[matr[loop]].m; 
         parity *= Z_orbit.mbas[matr[loop]].par; 
       }  
       switch(Info_pn[func].Z_par * parity)  {
          case +1: if((m <= Info_pn[func].Zm_max[0]) && (m >= Info_pn[func].Zm_min[0])) num++;
                   break;
          case  0: if(     (parity == +1) && (m <= Info_pn[func].Zm_max[0])
                                         && (m >= Info_pn[func].Zm_min[0])) num++;
                   else if((parity == -1) && (m <= Info_pn[func].Zm_max[1])
                                         && (m >= Info_pn[func].Zm_min[1])) num++;
          case -1: break;
       } /* end of switch loop */
     } 
      /* new particle configuration */

      for(loop = 0; (loop < Info_pn[func].Z_part) && ((matr[loop]+1) >= matr[loop+1]); loop++);
      matr[loop]++;
      for(loop--; loop >= 0; loop--)  {
         matr[loop] = loop;
      }
   } while(matr[Info_pn[func].Z_part] == Info_pn[func].Z_orb);  /* end of proton loop */

   if(num == 0)  {
      printf("\nError in the number of proton configurations for case = %d", func);
      printf("  number of SD = %d", num);
      printf("\n Message from function number_of_proton_SD()\n");
      exit(1);
   }
   free(matr);  /* Release temporary memory */
   return num;
} /* End: function number_of_proton_SD() */

  /*
   * The function                                                 
   *           proton_SD_configuration()                       
   * calculates all proton SD configurations for initial case (func = 0) or 
   * final case (func = 1). Then the SD are organized into groups 
   * after m-value and parity. The first group has max. m-value and + parity, 
   * next max. m-value and - parity, down to min. m-value and - parity.
   * Some of the group may not contain any SD.
   The sequence of primary SD are finally stored in SD_Z[ ]. 
   */

void proton_SD_configuration(int func)
{
  register int     loop, num, group, m, m_lower_limit, phase, start_pos;
  int              *matr;
  ULL              *st_ptr, high;
  struct sd        *sd_Z,  *sd_ptr;
  struct group     *gr;

  gr  = Group[func];    /* initialization of pointer */

    /*
     * Allocate temporary memory to store SD with parity 
     * and m-value plus particle bit configurations.       
     */

  sd_Z = (struct sd *)vector(sizeof(struct sd),Info_pn[func].Z_num + 1);

  matr = (int *) vector(sizeof(int),Info_pn[func].Z_part + 1);
  for(loop = 0; loop < Info_pn[func].Z_part; loop++) matr[loop] = loop; /* lowest config. */
  matr[Info_pn[func].Z_part] = Info_pn[func].Z_orb;
  
  sd_ptr = sd_Z;     /* initialization of pointer */
  num    = 0;        /* check the number of proton SD */
  do  {               /* the particle loop */
    sd_ptr->m_value = 0;
    sd_ptr->parity  = +1;
    sd_ptr->config  = ULL_ZERO;
    for(loop = 0; loop < Info_pn[func].Z_part; loop++)  {
      sd_ptr->m_value += Z_orbit.mbas[matr[loop]].m;
      sd_ptr->parity  *= Z_orbit.mbas[matr[loop]].par;
      sd_ptr->config  += ULL_ONE << matr[loop];
    }

    // check for particle limitations

    for(loop = 0; loop < Z_orbit.numMask; loop++) {
      high = sd_ptr->config & Z_orbit.list[loop];
      phase = 0;
      while(high) {
	high &= high - 1;
	phase++;
      }
      if(  (phase < Z_orbit.lim[2*loop])
	 ||(phase > Z_orbit.lim[2*loop + 1])) break;
    }
       // end  Masking on current config

    if(loop == Z_orbit.numMask) {
    
       /* Check if the parity and m--value is acceptable */

      switch(Info_pn[func].Z_par * sd_ptr->parity)  {
	case +1: if(   (sd_ptr->m_value <= Info_pn[func].Zm_max[0]) 
                       && (sd_ptr->m_value >= Info_pn[func].Zm_min[0]))  {
	           sd_ptr++;
                   num++;
	         } 
                   break;
        case  0: if((sd_ptr->parity == +1) && (sd_ptr->m_value <= Info_pn[func].Zm_max[0])
                                           && (sd_ptr->m_value >= Info_pn[func].Zm_min[0]))  {
                   sd_ptr++;
                   num++;
                 }    
	         else if((sd_ptr->parity == -1) && (sd_ptr->m_value <= Info_pn[func].Zm_max[1])
                                                && (sd_ptr->m_value >= Info_pn[func].Zm_min[1])){
                   sd_ptr++;
                   num++;
                 }
	case -1: break;
      } /* end of switch loop */
    }
            /* new proton configuration */

    for(loop = 0; (loop < Info_pn[func].Z_part) && ((matr[loop]+1) >= matr[loop+1]); loop++);
    matr[loop]++;
    for(loop--; loop >= 0; loop--)  {
      matr[loop] = loop;
    }
  } while(matr[Info_pn[func].Z_part]== Info_pn[func].Z_orb);  /* end of particle loop */

    /* Check the number of proton SD */

  if(num != Info_pn[func].Z_num)  {
    printf("\nError in calculating the proton SD for case = %d", func);
    printf("\n From function number_of_proton_SD(): num_sd = %d",Info_pn[func].Z_num);
    printf("\n From function proton_SD_configuration(): num_sd = %d\n",num);
    exit(1);
  }

    /* Release temporary memory for the particle bit configurations.*/

  free(matr);

    /* All proton SD are sorted into groups with given M_Z value and with max M_Z first.*/

  qsort(sd_Z, (size_t) num,(size_t) sizeof(struct sd),
	(int(*)(const void  *,const void  *)) plus_m_value);

  group     = 0;    /* initialization of group counter, etc. */
  start_pos = 0;
  sd_ptr    = sd_Z;
  
     /* Stop test for the do-loop */

  m_lower_limit = MIN(Info_pn[func].Zm_min[0], Info_pn[func].Zm_min[1]) - 2;
  sd_Z[Info_pn[func].Z_num].m_value = m_lower_limit; 
  
  do  {
    m = sd_ptr->m_value;
    num = 0;
    while( sd_ptr[++num].m_value == m);
    
         /* 
          * Sort the SD with M_Z = m into two parity groups with 
          * parity = +1 first.
          */

    if(num > 1) qsort(sd_ptr,(size_t) num,(size_t) sizeof(struct sd),
                          (int(*)(const void *,const void *)) plus_parity);
        
    if(sd_ptr[0].parity == -1)   {  
          
          /* No proton SD with current m_value has + parity */
        
       gr[group].Z_m       =  m;   /* store the necessary data */
       gr[group].Z_par     = +1;
       gr[group].Z_start   =  0;
       gr[group++].Z_numSD =  0;
    
          /* All proton SD with current m_value has - parity */
    
       gr[group].Z_m       =  m;   /* store the necessary data */
       gr[group].Z_par     = -1;
       gr[group].Z_start   =  start_pos;
       gr[group++].Z_numSD =  num;
       
         /* sort the proton SD with - parity after increasing configuration number.*/

       if(num > 1) qsort(sd_ptr,(size_t) num,(size_t) sizeof(struct sd),
                         (int(*) (const void  *,const void  *)) config_comp);
    
      sd_ptr    += num;  /* points to beginning of next SD group */
      start_pos += num;
    
    } /* end if-case with only one type of parity SD */       
    
    else   {
       
        /* find number of proton SD with + parity */
      
      for( loop = 0; ((sd_ptr[loop].parity) == +1) && (loop < num); loop++);
      
          /* Store the necessary data */

      gr[group].Z_m       = m;   /* store the necessary data */
      gr[group].Z_par     = +1;
      gr[group].Z_start   = start_pos;
      gr[group++].Z_numSD = loop;

       /* sort proton SD with + parity after increasing configuration number.*/

      if(loop > 1) qsort(sd_ptr,(size_t) loop,(size_t) sizeof(struct sd),
                         (int(*) (const void  *,const void  *)) config_comp);

         /* next proton SD group with - parity */

      sd_ptr    += loop;  /* points to beginning of - parity SD */
      start_pos += loop;
         
         /* Store the necessary data */

      gr[group].Z_m       =  m;   /* store the necessary data */
      gr[group].Z_par     = -1;
      gr[group].Z_start   =  (loop == num) ? 0 : start_pos;
      gr[group++].Z_numSD =  num - loop;

       /* sort proton SD with - parity after increasing configuration number.*/

      if((num - loop > 1)) qsort(sd_ptr,(size_t)(num - loop),(size_t) sizeof(struct sd),
                                (int(*)(const void  *,const void  *)) config_comp);
      
      start_pos += num - loop;   /* points to beginning of next SD */
      sd_ptr    += num - loop;
    } /* end if-case with both type of parity SD */     
  } while(sd_ptr->m_value > m_lower_limit);     

  if((Info_pn[func].num_gr) != group)  {
     printf("\n\nError from function proton_SD_configuration: Group number is wrong!!");
     printf("\nFrom function number_of_groups() - group number = %d",Info_pn[func].num_gr);
     printf("\nFrom function proton_SD_configuration() - group number = %d",group);
     printf("\n for case = %d\n",func);
     exit(1);
   }
    /*
     * Transform all SD configurations to permanent storage in  
     * SD_Z[] and release temporary memory allocation.
     */

  sd_ptr = sd_Z;                 /* initialization */
  st_ptr = SD_Z[func];
  for(loop = 0; loop < Info_pn[func].Z_num; loop++) {
    *(st_ptr++) = (sd_ptr++)->config;
  }
  free(sd_Z);

} /* End: function proton_SD_configuration() */  

  /*
   * The function                                         
   *          number_of_neutron_SD()                      
   * calculates and return the number of neutron SD for the initial 
   * case (func = 0) or the final case (func = 1)
   * and store the result in Info_pn[func].N_num.       
   */

int number_of_neutron_SD(int func)
{
   register int     loop, num, m, phase, parity;
   int              *matr;
   ULL              config, high;

   matr = (int *)vector(sizeof(int),Info_pn[func].N_part + 1); /* temporary memory */

      /* the lowest configuration */

   for(loop = 0; loop < Info_pn[func].N_part; loop++) matr[loop] = loop; /* lowest config. */
   matr[Info_pn[func].N_part] = Info_pn[func].N_orb;
   num = 0;   /* initialization */
   do  {   /* the particle  loop */

    // test for Masking on current config
     
     // check for particle limitations
  	        // code the configuration

     for(config = ULL_ZERO, loop = 0; loop < Info_pn[func].N_part; loop++)  {
       config += ULL_ONE << matr[loop];
     }
     for(loop = 0; loop < N_orbit.numMask; loop++) {
       high = config & N_orbit.list[loop];
       phase = 0;
       while(high) {
	 high &= high - 1;
	 phase++;
       }
       if(  (phase < N_orbit.lim[2*loop])
	    ||(phase > N_orbit.lim[2*loop + 1])) break;
     }
 
    // end  Masking on current config

     if(loop == N_orbit.numMask) {
       for(parity = +1, m = 0, loop = 0; loop < Info_pn[func].N_part; loop++)  {
	 m      += N_orbit.mbas[matr[loop]].m; 
	 parity *= N_orbit.mbas[matr[loop]].par; 
       }  
       switch(Info_pn[func].N_par * parity)  {
         case +1: if((m <= Info_pn[func].Nm_max[0]) && (m >= Info_pn[func].Nm_min[0])) num++;
                  break;
         case  0: if(     (parity == +1) && (m <= Info_pn[func].Nm_max[0])
                                        && (m >= Info_pn[func].Nm_min[0])) num++;
                  else if((parity == -1) && (m <= Info_pn[func].Nm_max[1])
                                        && (m >= Info_pn[func].Nm_min[1])) num++;
         case -1: break;
       } /* end of switch loop */
     }
      /* new particle configuration */

      for(loop = 0; (loop < Info_pn[func].N_part) && ((matr[loop]+1) >= matr[loop+1]); loop++);
      matr[loop]++;
      for(loop--; loop >= 0; loop--)  {
         matr[loop] = loop;
      }
   } while(matr[Info_pn[func].N_part] == Info_pn[func].N_orb);  /* end of particle loop */

   if(num == 0)  {
      printf("\n Error in the number of neutron configurations -- num = %d", num);
      printf("\n Message from function number_of_neutron_SD()");
      printf("\nfor case = %d\n",func);
      exit(1);
   }
   free(matr);  /* Release temporary memory */
   return num;
} /* End: function number_of_neutron_SD() */

  /*
   * The function                                                 
   *           neutron_SD_configuration()                       
   * calculates all neutron SD configurations. Then the  SD are organized  into groups 
   * after m-value and parity. The first group has min. m_value parity = total parity, 
   * nex min. m_value and parity = - total parity, down to max. m_value and
   * parity = - total parity.
   * Some groups may not contain any SD.
   * The sequence of neutron SD are finally stored in SD_N[ ]. 
   */

void neutron_SD_configuration(int func)
{
   register int     loop,  num, group, m, m_max_limit, phase, start_pos, par;
   int              *matr;
   ULL               *st_ptr, high;
   struct sd        *sd_N,  *sd_ptr;
   struct group     *gr;

   gr  = Group[func];

    /*
     * Allocate temporary memory to store SD with parity 
     * and m-value plus particle bit configurations.       
     */

  sd_N = (struct sd *)vector(sizeof(struct sd),Info_pn[func].N_num + 1);

  matr = (int *) vector(sizeof(int), Info_pn[func].N_part + 1);

      /* the lowest particle configuration */

  for(loop = 0; loop < Info_pn[func].N_part; loop++) matr[loop] = loop; /* lowest config. */
  matr[Info_pn[func].N_part] = N_orbit.num;

  sd_ptr = sd_N;     /* initialization of pointer */
  num    = 0;        /* check the number of neutron SD */
  do  {               /* the particle loop */
    sd_ptr->m_value = 0;
    sd_ptr->parity  = +1;
    sd_ptr->config  = ULL_ZERO;
    for(loop = 0; loop < Info_pn[func].N_part; loop++)  {
      sd_ptr->m_value += N_orbit.mbas[matr[loop]].m;
      sd_ptr->parity  *= N_orbit.mbas[matr[loop]].par;
      sd_ptr->config  += ULL_ONE << matr[loop];
    }
       // check for particle limitations

    for(loop = 0; loop < N_orbit.numMask; loop++) {
      high =    sd_ptr->config & N_orbit.list[loop];
      phase = 0;
      while(high) {
	high &= high - 1;
	phase++;
      }
      if(  (phase < N_orbit.lim[2*loop])
	 ||(phase > N_orbit.lim[2*loop + 1])) break;
    }
    // end  Masking on current config

    if(loop ==  N_orbit.numMask) {

      /* Check if parity and m-value is acceptable */

      switch(Info_pn[func].N_par * sd_ptr->parity)  {
        case +1: if(   (sd_ptr->m_value <= Info_pn[func].Nm_max[0]) 
		    && (sd_ptr->m_value >= Info_pn[func].Nm_min[0]))  {
	           sd_ptr++;
                   num++;
                 }
	         break;
        case  0: if((sd_ptr->parity == +1) && (sd_ptr->m_value <= Info_pn[func].Nm_max[0])
		                           && (sd_ptr->m_value >= Info_pn[func].Nm_min[0]))  {
                   sd_ptr++;
                   num++;
                 }    
                 else if((sd_ptr->parity == -1) && (sd_ptr->m_value <= Info_pn[func].Nm_max[1])
                                                && (sd_ptr->m_value >= Info_pn[func].Nm_min[1])){
                   sd_ptr++;
                   num++;
                 }
        case -1: break;
      } /* end of switch loop */
    }
      /* new particle configuration */

    for(loop = 0; (loop < Info_pn[func].N_part) && ((matr[loop]+1) >= matr[loop+1]); loop++);
    matr[loop]++;
    for(loop--; loop >= 0; loop--)  {
      matr[loop] = loop;
    }
  } while(matr[Info_pn[func].N_part]== Info_pn[func].N_orb);  /* end of particle loop */

    /* Check the number of neutron SD */

  if(num != Info_pn[func].N_num)  {
    printf("\nError in calculating the neutron SD");
    printf("\n From function number_of_neutron_SD(): num = %d", Info_pn[func].N_num);
    printf("\n From function neutron_SD_configuration(): num = %d\n",num);
    exit(1);
  }

    /* Release temporary memory for the particle bit configurations.*/

  free(matr);

   /* All neutron SD are sorted into groups with given M_N value and with min M_N first.*/

  qsort(sd_N, (size_t) num,(size_t) sizeof(struct sd),
            (int(*)(const void  *,const void  *)) minus_m_value);

  group     = 0;    /* initialization of group counter, etc. */
  par       = Info_com[func].par;  /* parity of first group of neutron SD */
  start_pos = 0;
  sd_ptr    = sd_N;
  
       /* Stop test for the do-loop */

  m_max_limit = MAX(Info_pn[func].Nm_max[0], Info_pn[func].Nm_max[1]) + 2;
  sd_N[Info_pn[func].N_num].m_value = m_max_limit; /* stop test for later loop */
  
  do  {
     m = sd_ptr->m_value;
     num = 0;
     while( sd_ptr[++num].m_value == m);
       
         /* 
          * Sort the SD with M_N = m into two parity groups with 
          * parity = total parity first.
          */

     if(num > 1)  { 
        if(par == +1)  qsort(sd_ptr,(size_t) num,(size_t) sizeof(struct sd),
                          (int(*)(const void *,const void *)) plus_parity);
        else          qsort(sd_ptr,(size_t) num,(size_t) sizeof(struct sd),
                        (int(*)(const void *,const void *)) minus_parity);
     }
     if(sd_ptr[0].parity == -par)   {

          /* No SD with current m_value has parity = total parity */

        gr[group].N_m       = m;   /* store the necessary data */
        gr[group].N_par     = par;
        gr[group].N_start   = 0;
        gr[group++].N_numSD = 0;

        /* All SD with current m_value has parity = - total parity */

        gr[group].N_m       = m;   /* store the necessary data */
        gr[group].N_par     = -par;
        gr[group].N_start   = start_pos;
        gr[group++].N_numSD = num;
      
          /* Sort the SD with parity = - total parity after increasing config. number */

        if(num > 1) qsort(sd_ptr,(size_t) num,(size_t) sizeof(struct sd),
                         (int(*) (const void  *,const void  *)) config_comp);
   
        sd_ptr    += num;   /* point to beginning of next group */
        start_pos += num;
     
     }  /* end if-case with only one type of parity SD */
     
     else   {
        
           /* find number of SD with parity = par */

        for( loop = 0; ((sd_ptr[loop].parity) == par) && (loop < num); loop++);

        gr[group].N_m       = m;   /* store the necessary data */
        gr[group].N_par     = par;
        gr[group].N_start   = start_pos;
        gr[group++].N_numSD = loop;

           /* sort the SD with parity = par after increasing configuration number.*/

        if(loop > 1) qsort(sd_ptr,(size_t) loop,(size_t) sizeof(struct sd),
                         (int(*) (const void  *,const void  *)) config_comp);

           /* next group of SD with parity = - par */

        sd_ptr    += loop;  /* point to beginning of next group */
        start_pos += loop;

        gr[group].N_m       = m;   /* store the necessary data */
        gr[group].N_par     = - par;
        gr[group].N_start   = start_pos;
        gr[group++].N_numSD = num - loop;

         /* sort the SD with  parity = - par after increasing configuration number.*/

        if((num - loop) > 1) qsort(sd_ptr,(size_t)(num - loop),(size_t) sizeof(struct sd),
                         (int(*)(const void  *,const void  *)) config_comp);
   
        start_pos += num - loop;   /* point to beginning of next group */
        sd_ptr    += num - loop;
     } /* end if-case with both types of parity SD */    
  
  } while(sd_ptr->m_value < m_max_limit);     

  if(Info_pn[func].num_gr != group)  {
     printf("\n\nError from function neutron_SD_configuration(): Group number is wrong!");
     printf("\nFrom function number_of_groups() - group number = %d", Info_pn[func].num_gr);
     printf("\nFrom function neutron_SD_configuration()-group number = %d\n",group);
     exit(1); 
  }
    /*
     * Transform all SD configurations to permanent storage in  
     * SD_B[] and release temporary memory allocation.
     */

  sd_ptr = sd_N;      /* initialization */
  st_ptr = SD_N[func];
  for(loop = 0; loop < Info_pn[func].N_num; loop++) {
    *(st_ptr++) = (sd_ptr++)->config;
  }
  free(sd_N);

} /* End: function neutron_SD_configuration() */

    /*
     * The function                                          
     *       void  pn_particle_transitions()                
     * calculate for proton - neutron case all one-particle transition matrix
     * elements  between eigenstates generated through a Lanczos procedure.
     */

void pn_particle_transitions(int trans_type)
{
   int       trans, loop, init_num, init_save, final_num,
             j_trans, del_m, del_p;
   double    clebsch, value = D_ZERO;

         /*
         ** Calculate the group change between initial and final states
         ** performed by the current operator for both proton operator
         ** OP(p,p') and neutron operator OP(n,n')
         */

   memory_allocation(GROUP_CHANGE,0); 

   del_m = Model[1].MJ - Model[0].MJ;
   del_p = (Model[0].P == Model[1].P) ? +1 : -1;   

                /* calculate group changes for protons */

   group_change(Info_pn[0].num_gr, Info_pn[1].num_gr, del_m, del_p, 0, +1, Group_change_p);

                /* calculate group changes for neutrons */

   group_change(Info_pn[0].num_gr, Info_pn[1].num_gr, 0, +1, del_m, del_p, Group_change_n);
                
         /*
         ** If initial and final |SD> basis are equal  calculate
         ** and store all diagonal matrix elements in SD_diag[] 
         */

   if(Model[0].same_basis == YES) pn_diag_matrix_elements();

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
             fprintf(Files.out_data," B(M1) in units (e*hbar)/(2*M*c)");
             break;
       case  E3_TRANSITION: 
             j_trans = 6;  /* twice angular momentum transfer  */
             fprintf(Files.out_data,"\n\n            E3 transitions:");
             fprintf(Files.out_data,"\nInit(no 2J) Final(no 2J) (gr.st.no=0):");
             fprintf(Files.out_data,"  B(E3) in  (e**2) * (b**6)");
             break;
       default:
          printf("\n\nError in function pn_particle_transitions():");
          printf("\n   transition type = %d is wrong!!!!!\n",trans_type);
          exit(1);
   } /* end switch() */

    /* Start the calculation of the transition matrix elements.*/

   for(trans = 0,init_save = -1; trans < Model[0].num_trans; trans++) {

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

         read_pn_eigenvector(File_eigen_init, init_num, Info_com[0].num, Eigen_vec); 

                   /* initialize the vector Trans_contr to ZERO. */

         for(loop = 0; loop < Info_com[1].num; loop++)  Trans_contr[loop] = D_ZERO;

                   /* Calculate: |Trans_contr> = OP * |Eigen_vec> */

         if(Model[0].same_basis == YES) 
             diag_contr(Info_com[0].num, 1, SD_diag, Eigen_vec, Trans_contr);
        pn_nondiag_iteration(Eigen_vec, Trans_contr);

      }  /* end if(new initial eigenvector) */

      final_num = Model[0].trans[trans].num_final;

      read_pn_eigenvector(File_eigen_final, final_num, Info_com[1].num, Eigen_vec); 

      value = pn_overlap(Info_com[1].num, Trans_contr, Eigen_vec) / clebsch;
   
      value *= value * (Model[0].trans[trans].j_final +1)/(Model[0].trans[trans].j_init +1);
      fprintf(Files.out_data,"\n    %2d %2d        %2d %2d            ", /* Print result */
          init_num, Model[0].trans[trans].j_init,final_num, Model[0].trans[trans].j_final);
      fprintf(Files.out_data,"               %12.8E",value);
      init_save = init_num; /* saved for next transition */
   } /* all transitions calculated */
} /* End: function pn_particle_transitions() */

      /*
      ** The function
      **      group_change()
      ** calculates and stores in Group_change_p[] and Group_change_n[]
      ** the group changes performed by the one-particle operator
      ** OP(p,p') and OP(n,n') fro initial to final states
      */

void group_change(int num_gr_i, int num_gr_f, int del_mz, int del_pz,
                                 int del_mn, int del_pn, int *result)
{
   int    gr_i, gr_f, m_zf, p_zf, m_nf, p_nf;

   for(gr_i = 0; gr_i < num_gr_i; gr_i++)  {
      m_zf = Group[0][gr_i].Z_m + del_mz;
      p_zf = Group[0][gr_i].Z_par * del_pz;
      m_nf = Group[0][gr_i].N_m + del_mn;
      p_nf = Group[0][gr_i].N_par * del_pz;

      for(gr_f = 0; gr_f < num_gr_f; gr_f++)  {
         if(   (m_zf == Group[1][gr_f].Z_m)
            && (p_zf == Group[1][gr_f].Z_par)
            && (m_nf == Group[1][gr_f].N_m)
            && (p_nf == Group[1][gr_f].N_par))   break;
      } /* end loop through all final groups */

      result[gr_i] = (gr_f < num_gr_f) ? gr_f : -1;

   } /* end loop through all init groups */
   
   return;
 } /* End: function group_change() */

  /*
   * The function                                               
   *            pn_diag_matrix_elements()            
   * calculate diagonal matrix elements <SD(N),SD(Z)| OP | SD(Z),SD(N)> and
   * store the result in the diagonal matrix element vector SD_diag[ ]. It 
   * is assumed that Model.Z * Model.N is different from ZERO.
   */

void pn_diag_matrix_elements(void)
{
   int          orb, par, loop, group;
   ULL           pos, *sd_ptr;
   double       *p_matrix, *n_matrix, *sd_diag; 

       /* maximum proton and neutron configuration in a group */

   orb = Group[0][0].Z_numSD;
   par = Group[0][0].N_numSD; 

   for(group = 1; group < Info_pn[0].num_gr; group++)  {
      orb = MAX(orb,Group[0][group].Z_numSD);
      par = MAX(par,Group[0][group].N_numSD); 
   }  
        /* save the two max. values for later use */
  
   Info_pn[0].max_Z_SD = orb;
   Info_pn[0].max_N_SD = par;

   for(group = 1; group < Info_pn[1].num_gr; group++)  {
      orb = MAX(orb,Group[1][group].Z_numSD);
      par = MAX(par,Group[1][group].N_numSD); 
   }  
        /* save the two max. values for later use */
  
   Info_pn[1].max_Z_SD = orb;
   Info_pn[1].max_N_SD = par;

       /* local memory to store <SD|OP|SD> for protons and neutrons */

   p_matrix = (double *) vector(orb, sizeof(double));
   n_matrix = (double *) vector(par, sizeof(double));

       /* run through all proton-neutron groups */

   for(group = 0, sd_diag = SD_diag; group < Info_pn[0].num_gr; group++)  {

         /* any proton and neutron |SD> in the present group */

      if(!(Group[0][group].Z_numSD * Group[0][group].N_numSD)) continue;  /* no ! */

           /* run through all|SD(Z)> in the group and calculate <SD(Z)|OP(p)|SD(Z)> */

      for(loop = 0, sd_ptr = SD_Z[0] + Group[0][group].Z_start; loop < Group[0][group].Z_numSD; 
                                                                   loop++, sd_ptr++)  {
         p_matrix[loop]  = 0.0;             /* initialization */
         pos             = ULL_ONE;
         orb             = 0;
         par             = 0;
         do   {            /* particle loop */
            for( ; !((*sd_ptr) & pos); orb++, pos <<= 1);
            p_matrix[loop] += Table_pn_l[0][orb++].diag;
            pos <<= 1;
         } while(++par < Info_pn[0].Z_part);  /* end proton particle loop */

      } /* end loop through all proton SD in a group */        
    
         /* run through all |SD(N)> in the group and calculate <SD(N)|OP(n)|SD(N)> */

      for(loop = 0, sd_ptr = SD_N[0] + Group[0][group].N_start; loop < Group[0][group].N_numSD;
                                                                 loop++, sd_ptr++)  {

         n_matrix[loop]  = 0.0;             /* initialization */
         pos             = ULL_ONE;
         orb             = 0;
         par             = 0;
         do   {            /* particle loop */
            for( ; !((*sd_ptr) & pos); orb++, pos <<= 1);
            n_matrix[loop] += Table_pn_l[1][orb++].diag;
            pos <<= 1;
         } while(++par < Info_pn[0].N_part);  /* end neutron particle loop */

      } /* end loop through all neutron SD in a  group */

                /* calculate all <SD(N),SD(Z)| OP(p) + OP(n)|SD(Z), SD(N)> */

      for(orb = 0; orb < Group[0][group].Z_numSD; orb++)  {
         for(par = 0; par < Group[0][group].N_numSD; par++)  {
            *(sd_diag++) = p_matrix[orb] + n_matrix[par];
	 } /* end neutron loop */
      } /* end proton loop */
   
   } /* end loop through all groups */

  free(n_matrix);   /* release temporary memory */
  free(p_matrix);

} /* End: function pn_diag_matrix_elements() */

    /*
     * The function            
     *    void read_pn_eigenvector()
     * reads the eigenvector number eigen_num from file_name into vec_ptr[]
     * genrated from the Lanczos proton-neutron program.
     */

void read_pn_eigenvector(char *file_name,int eigen_num, int elements, double *vec_ptr)
{
   int             n; 
   FILE    *file_ptr;

   printf("\n\n\nStart reading eigenvector number = %d\n\n",eigen_num);

    /* Open the file  for the  eigenvectors. */

   if((file_ptr = fopen(file_name,"rb")) ==  NULL_PTR)  {
      printf("\n\nError in function read_pn_eigenvector:");
      printf("\nWrong eigenvector file-name = %s\n", file_name);
      exit(1);
   }

   rewind(file_ptr);

   fseek(file_ptr,(long)(eigen_num * (2*sizeof(int) +  elements * sizeof(double))),SEEK_SET);

   if(fread((void *)&n,(size_t)sizeof(int),1,file_ptr) != 1) {/* read  vector number*/
      printf("\n\nError in function read_pn_eigenvector():");
      printf("\nSomething is wrong in reading the vector number n = %d\n",n);
      exit(1);
   }
   if(eigen_num != n)  {
      printf("\n\nError in function read_pn_eigenvector():");
      printf("\nWrong eigenvector number in file %s",file_name );
      printf("\n The correct eigenvector has not been found");
      printf("\neigenvector number = %d   number read = %d\n",eigen_num,n);
      exit(1);
   } /* end of if-test */

  if(fread((void *)&n,(size_t)sizeof(int),1,file_ptr) != 1) {/* read  vector dimension */
      printf("\n\nError in function read_pn_eigenvector():");
      printf("\nSomething is wrong in reading the dimension of the vector number n = %d\n",n);
      exit(1);
   }
   if(elements != n)  {
      printf("\n\nError in function read_pn_eigenvector():");
      printf("\nWrong dimension of eigenvector number %d  in file %s",eigen_num, file_name );
      printf("\n The correct dimension = %d   number read = %d\n",elements,n);
      exit(1);
   } /* end of if-test */

   if(fread((void *) vec_ptr, (size_t) sizeof(double),        /* read the vector */
           (size_t) elements, file_ptr) != (size_t) elements)    {
      printf("\n\n Error in function read_pneigenvector():");
      printf("\nWrong number of elements in file %s",file_name );
      printf("\nvector_num = %d elements = %d read = %d\n",eigen_num,elements,n);
      exit(1);
   }   /* end of if-test  */

   fclose(file_ptr);
   return;

} /* End: function read_pn_eigenvector() */

  /*
   * The function                                 
   *        pn_nondiag_iteration()                  
   * calculate and store contribution from the non-diagonal matrix 
   * elements <SD'(N),SD'(Z)|OP(p,p') + OP(n.n')|SD(Z),SD(N)> to
   * vector OP * |init_vec> --> |final_vec>   
   */

void pn_nondiag_iteration(double *init_vec, double *final_vec)
{
   int        group, num_zf, num_nf, num_zi, num_ni;
   double    *init_start, *final_start;

        /* run through all proton-neutron groups.*/

   for(group = 0; group < Info_pn[0].num_gr; group++)  {

      num_zi = Group[0][group].Z_numSD;             /* initialization */ 
      num_ni = Group[0][group].N_numSD;
      init_start = init_vec + Group[0][group].ampl_start;
      if(Group_change_p[group] != -1)   {                     /* proton contribution */

         num_zf = Group[1][Group_change_p[group]].Z_numSD;      /* initialization */   
         num_nf = Group[1][Group_change_p[group]].N_numSD;
         if(num_nf != num_ni)   {      /* test check */
            printf("\n\nError in function pn_nondiag_iteration():");
            printf("\nNumber of |SD(N)> in groups %d and %d should be equal",
                                                  group, Group_change_p[group]);
            printf("/n num_ni(%d) = %d and num_nf(%d) = %d\n", group, num_ni,
                                                  Group_change_p[group], num_nf);
            exit(1);
         }
           /* contribution from init to final group by proton operator OP(p,p') */

         final_start = final_vec + Group[1][Group_change_p[group]].ampl_start;

         nondiag_proton_contribution
             (Info_pn[0].Z_part, num_zi, SD_Z[0] + Group[0][group].Z_start,
                 num_zf, SD_Z[1] + Group[1][Group_change_p[group]].Z_start,
                                          num_nf, init_start, final_start);
      }  /* end proton contribution */

      if(Group_change_n[group] != -1)   {

         num_zf = Group[1][Group_change_n[group]].Z_numSD;      /* initialization */   
         num_nf = Group[1][Group_change_n[group]].N_numSD;
         if(num_nf != num_ni)   {      /* test check */
            printf("\n\nError in function pn_nondiag_iteration():");
            printf("\nNumber of |SD(N)> in groups %d and %d should be equal",
                                                  group, Group_change_n[group]);
            printf("/n num_ni(%d) = %d and num_nf(%d) = %d\n", group, num_ni,
                                                  Group_change_n[group], num_nf);
            exit(1);
         }
           /* contribution from init to final group by proton operator OP(p,p') */

         final_start = final_vec + Group[1][Group_change_n[group]].ampl_start;

         nondiag_neutron_contribution
              (Info_pn[0].N_part, num_ni, SD_N[0] + Group[0][group].N_start,
                  num_nf, SD_N[1] + Group[1][Group_change_n[group]].N_start,
                                           num_zf, init_start, final_start);
      }  /* end neutron contribution */
   }  /* end loop through all initial groups */ 

} /* End: function pn_nondiag_iteration() */

    /*
    ** The function                                 
    **        nondiag_proton_contribution()                  
    ** calculates non-diagonal proton one-particle contributions 
    ** from all initial states amp * |SD[Z]_i: SD[N],i> and 
    ** stores the result as a contribution in the final state
    **  amp_f * |SD[Z]_f: SD[N],f>
    */

void nondiag_proton_contribution(int part_num, int num_zi,  ULL *sd_zi,
                              int num_zf, ULL *sd_zf, int num_nf,
                                      double *init, double *final)
{
   int                     nf, loop, part, orb, phase, search_num, low, high;
   ULL                      pos, l_sd, new_sd, test;   
   struct  trans_nondiag   *l_ptr;

         /* run through all initial |SD(Z)> in present group */

   for( loop = 0; loop < num_zi; loop++, sd_zi++) {

      for(part = 0, orb = 0, pos = ULL_ONE; part < part_num; part++, orb++, pos <<= 1) {
         for( ; !((*sd_zi) & pos); orb++, pos <<= 1);   /* l loop */
         if((l_ptr = Table_pn_l[0][orb].nondiag) == NULL_PTR) continue;
         l_sd    = (*sd_zi) ^ pos;                                 /* l part. found */
         search_num = 0;                                     /* search in final |SD> */
         do   {                     /* run through all one-particle matrix elements */

            for( ;(l_sd & l_ptr->one); l_ptr++);               /* Pauli principle */
            if(!(l_ptr->one)) break;                     /* no more matrix elements */
            new_sd = l_sd ^ l_ptr->one;              /* contribution found, new |SD> */

	    //   Mask test on |new_sd>
	    {
	      int  num, phase;
	      ULL  high;

	      if(Z_orbit.numMask > 0) {
		for(num = 0; num < Z_orbit.numMask; num++) {
		  high = new_sd & Z_orbit.list[num];
		  phase = 0;
		  while(high) {
		    high &= high - 1;
		    phase++;
		  }
		  if(  (phase < Z_orbit.lim[2*num])
		       ||(phase > Z_orbit.lim[2*num + 1])) break;
		}
		if(num < Z_orbit.numMask) continue;
	      }
	    }
	    // end Mask test on   |new_sd>

            test    = new_sd & l_ptr->two;        /* Calculate the l  permutation phase. */
            phase = +1;
            while(test) {
               test &= test - 1;
               phase = -phase;
            }
            low  = search_num;        /* binary search */
            high = num_zf;
            while(1)  {
               search_num = (low + high)>> 1;
               if(new_sd < sd_zf[search_num])      high = search_num - 1;
               else if(new_sd > sd_zf[search_num]) low  = search_num + 1;
               else                                   break;
            }
              /* contribution |sd_init_ptr[loop]> --> |sd_final[search_num]> */

           for(nf = 0; nf < num_nf; nf++)  {
              final[search_num * num_nf + nf]
                 += phase * l_ptr->val * init[loop * num_nf + nf];
           } 
         } while((++l_ptr)->one); /* end do-loop for contribution from all matr. elem */

      }  /* end l particle loop */
   
   } /* end loop through all |sd_zi[loop]>*/

} /* End: function  nondiag_proton_contribution() */

    /*
    ** The function                                 
    **        nondiag_neutron_contribution()                  
    ** calculates non-diagonal neutron one-particle contributions 
    ** from all initial states amp * |SD[Z]_i: SD[N],i> and 
    ** stores the result as a contribution in the final state
    **  amp_f * |SD[Z]_f: SD[N],f>
    */

void nondiag_neutron_contribution(int part_num, int num_ni,  ULL *sd_ni,
                              int num_nf, ULL *sd_nf, int num_zf,
                                      double *init, double *final)
{
   int                     zf, loop, part, orb, phase, search_num, low, high;
   ULL                      pos, l_sd, new_sd, test;   
   struct  trans_nondiag   *l_ptr;

         /* run through all initial |SD(N)> in present group */

   for( loop = 0; loop < num_ni; loop++, sd_ni++) {

      for(part = 0, orb = 0, pos = ULL_ONE; part < part_num; part++, orb++, pos <<= 1) {
         for( ; !((*sd_ni) & pos); orb++, pos <<= 1);   /* l loop */
         if((l_ptr = Table_pn_l[1][orb].nondiag) == NULL_PTR) continue;
         l_sd    = (*sd_ni) ^ pos;                                      /* l part. found */
         search_num = 0;                                               /* search in final |SD> */
         do   {                                /* run through all one-particle matrix elements */

            for( ;(l_sd & l_ptr->one); l_ptr++);                            /* Pauli principle */
            if(!(l_ptr->one)) break;                                /* no more matrix elements */
            new_sd = l_sd ^ l_ptr->one;                        /* contribution found, new |SD> */

	    //   Mask test on |new_sd>
	    {
	      int  num, phase;
	      ULL  high;

	      if(N_orbit.numMask > 0) {
		for(num = 0; num < N_orbit.numMask; num++) {
		  high = new_sd & N_orbit.list[num];
		  phase = 0;
		  while(high) {
		    high &= high - 1;
		    phase++;
		  }
		  if(  (phase < N_orbit.lim[2*num])
		       ||(phase > N_orbit.lim[2*num + 1])) break;
		}
		if(num < N_orbit.numMask) continue;
	      }
	    }
	    // end Mask test on   |new_sd>

            test    = new_sd & l_ptr->two;              /* Calculate the l  permutation phase. */
            phase = +1;
            while(test) {
               test &= test - 1;
               phase = -phase;
            }
            low  = search_num;        /* binary search */
            high = num_nf;
            while(1)  {
               search_num = (low + high)>> 1;
               if(new_sd < sd_nf[search_num])      high = search_num - 1;
               else if(new_sd > sd_nf[search_num]) low  = search_num + 1;
               else                                   break;
            }
              /* contribution |sd_init_ptr[loop]> --> |sd_final[search_num]> */

           for(zf = 0; zf < num_zf; zf++)  {
              final[zf * num_nf + search_num]
                 += phase * l_ptr->val * init[zf * num_nf + loop];
           } 
         } while((++l_ptr)->one); /* end do-loop for contribution from all matr. elem */

      }  /* end l particle loop */
   
   } /* end loop through all |sd_ni[loop]>*/

} /* End: function  nondiag_neutron_contribution() */

  /*
   * The function 
   *             pn_overlap()
   * calculates overlap between  |trans_contr> = OP * |init_vec>
   * and a given final eigenvector |final_eigen>
   */

double pn_overlap(int num_sd, double *trans_contr, double *final_eigen)
{
   int         number;
   double      *eigen_ptr, *op_ptr, result;

   result  = 0.0;    /* initialization */
   number  = num_sd;
   if(number)   {
      eigen_ptr = final_eigen;
      op_ptr    = trans_contr;
      do {
         result += (*(op_ptr++)) * (*(eigen_ptr++));
      } while(--number);

   } 
   return result;

} /* End: function pn_overlap()  */
