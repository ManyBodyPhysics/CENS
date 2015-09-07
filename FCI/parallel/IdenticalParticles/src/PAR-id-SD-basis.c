/*******************  The module PAR-lanc-id-SD-basis.c  ******************/

#include "PAR-shell-model.h"

  /*
   * The entrance function                                                      
   *  void id_slater_determinant(SHELL *model, SP_BAS *sp_bas, SD_BAS *tot_bas)
   * allocates necessary memory, calculates and store the complete
   * set of |SD> for a system of basis.part identical particles. 
   * If possible time-reversed symmetry is assumed.
   */
	    /* local function declarations */

static void number_of_trsym_SD(SP_BAS *sp_bas, SD_BAS *sd_bas);
   /*
   ** calculates and stores in basis->num_SD[] the number of slater determinants (SD) for 
   ** even number of particles n with MJ = 0. It calculates the number of asymmetric
   ** (num_SD[0]) and the number of symmetric (num_SD[1]) SD. If total angular
   ** momentum J is odd, num_SD[1] = 0:
   **             J odd -> j_value = FALSE. J even -> j_value = TRUE
   */

static void number_of_nosym_SD(SP_BAS *sp_bas, SD_BAS *sd_bas);
   /*
   ** calculates the number of slater determinants (SD) for n particles
   ** with given parity and MJ = m_value.
   */

static void tr_sym_SD_config(SP_BAS *sp_bas, SD_BAS *sd_bas);
  /*
  ** calculates and stores the complete set of slater determinants |SD()> for n
  ** particles with MJ = 0. It calculates the asymmetric |SD()> and the symmetric 
  ** |SD()> and stored in basis->SD[], asymSD[]> first, The number of asymmetric 
  ** |SD()> is checked against basis->numSD[0] and the number of symmetric |SD()> 
  ** against basis->numSD[1].
  ** If total angular momentum J is odd, no sym |SD> is stored:
  **             J odd -> basis->J_type = FALSE. J even -> basis->J_type = TRUE
  */

static void no_sym_SD_config(SP_BAS *sp_bas, SD_BAS *sd_bas);
   /*
   ** calculates and stores the complete set of slater determinants (SD) 
   ** for n identical particles with angular projection m_J and fixed parity (parity).
   ** The resulting SD are stored in conf_SD[]. The number of no-symmetric
   ** SD is checked against num_SD.
   */

            /**** End: function declarations ****/


               /**** The function definitions  ****/ 

  /*
  ** The entrance function                                                      
  **             id_slater_determinant()                          
  ** allocates necessary memory, calculates and store the complete
  ** set of |SD> for a system of basis.part identical particles. 
  ** If possible time-reversed symmetry is assumed.
  */

void id_slater_determinant(SHELL *model, SP_BAS *sp_bas, SD_BAS *sd_bas)
{
  char   *func = {"id_slater_determinant(): "};

  // initialyze data in structure sd_bas from structur model

  sd_bas->part     = model->part;
  sd_bas->P        = model->P;
  sd_bas->MJ       = model->MJ;
  sd_bas->numm_orb = sp_bas->numm_orb;
  sd_bas->J_type   = model->J_type;

  if(sd_bas->MJ == 0)  {   // time-reversal symmetry  
    number_of_trsym_SD(sp_bas, sd_bas);
    if((sd_bas->numSD[0] + sd_bas->numSD[1]) == 0) {
      printf("\n\nError in function id_slater_determinant():");
      printf("\nNo basis |asymSD()> or |symSD()> found!!\n");
      MPI_Abort(MPI_COMM_WORLD,Rank);
    }
    if(strcmp(model->type_calc,"dimension")) {
      sd_bas->SD 
	= MALLOC(sd_bas->numSD[0] + sd_bas->numSD[1], ULL, func, "sd_bas->SD[]");
      tr_sym_SD_config(sp_bas, sd_bas);
    }
  } /* end asym case */
  else  {               
    sd_bas->numSD[1] = 0;                       /* no time-reversal symmetry */
    number_of_nosym_SD(sp_bas, sd_bas);
    if(sd_bas->numSD[0] == 0) {
      printf("\n\nError in function id_slater_determinant():");
      printf("\nNo basis |nosymSD()> found!!\n");
      MPI_Abort(MPI_COMM_WORLD,Rank);
    }
    if(strcmp(model->type_calc,"dimension")) {
      sd_bas->SD = MALLOC(sd_bas->numSD[0], ULL, func, "_bas.SD[]");
      no_sym_SD_config(sp_bas, sd_bas);
    }
  }/* end nosym case */

  sd_bas->tot_dimSD = sd_bas->numSD[0] + sd_bas->numSD[1];     // total |SD> dimension

  if(Rank == MASTER) {
    FILE  *file_ptr;

    if((file_ptr = fopen(model->title,"a"))== NULL) {
      printf("\n\nError in function  id_slater_determinant():");
      printf("\nWrong file = %s for the output data\n",model->title);
      MPI_Abort(MPI_COMM_WORLD,Rank);
    }
    if(sd_bas->MJ == 0)  {
      fprintf(file_ptr,"\n\nTotal number of asymmetric slater determinants   =%8d",
                                                                sd_bas->numSD[0]);
      fprintf(file_ptr,"\nTotal number of symmetric slater determinants    =%8d",
                                                             sd_bas->numSD[1]);
      fprintf(file_ptr,"\nTotal number of slater determinants              =%8d",
	      2 * sd_bas->numSD[0] + sd_bas->numSD[1]);
    }
    else {
      fprintf(file_ptr,"\n\nTotal number of slater determinants = %d", sd_bas->numSD[0]);
    }
    fclose(file_ptr);

/**********************************************
    if(!strcmp(model->type_calc,"dimension")) exit(0);   dimension calculation only
***********************************************/
  } // end MASTER print-block

 } /* End: function slater_determinant() */

  /*
  ** The function
  **          number_of_trsym_SD()                   
  ** calculates and stores in sd_bas->num_SD[] the number of slater determinants |SD()>
  ** for even number of particles n with MJ = 0. It calculates the number of asymmetric
  ** (sd_bas->numSD[0]) and the number of symmetric (sd_bas->num_SD[1]) |SD>. 
  ** If total angular momentum J is odd, num_SD[1] = 0:
  **             J odd -> j_value = FALSE. J even -> j_value = TRUE
  */

void number_of_trsym_SD(SP_BAS *sp_bas, SD_BAS *sd_bas)
{
  char
         *func = {"number_of_trsym_SD(): "};
  int
         loop, m, par, last_orbit, *matr, num, phase;
  ULL
         config, time_reversed_config, high;

  matr = MALLOC(sd_bas->part + 1, int,func, "matr[]");

	      /* the lowest configuration */

  for(loop = 0; loop < sd_bas->part; loop++) matr[loop] = loop;
  matr[loop] = sp_bas->numm_orb;

       /* calculate the total number of SD. */

  sd_bas->numSD[0] = 0;  /* initialization */
  sd_bas->numSD[1] = 0;
  last_orbit = sp_bas->numm_orb - 1;

  do  {   /* run through all possible configurations */
    for(m = 0, par = +1, loop = 0; loop < sd_bas->part; loop++)  { /* m-value and parity */
      m   += sp_bas->mbas[matr[loop]].m;
      par *=  sp_bas->mbas[matr[loop]].par;
    }
    if((m == 0) && (par == sd_bas->P))  {
      
        /* Code the configuration  and its time-reversed one */

      config               = ULL_ZERO;
      time_reversed_config = ULL_ZERO;
      for(loop = 0; loop < sd_bas->part; loop++)  {
	config               += ULL_ONE << matr[loop];
	time_reversed_config += ULL_ONE <<(last_orbit - matr[loop]);
      }
      if(config <  time_reversed_config)  {
	if(sp_bas->mask.num == 0) {
	  sd_bas->numSD[0]++;
	}
	else {       // check for particle limitations
	  for(num = 0; num < sp_bas->mask.num; num++) {
	    high = config & sp_bas->mask.list[num];
	    phase = 0;
	    while(high) {
	      high &= high - 1;
	      phase++;
	    }
	    if(  (phase > sp_bas->mask.lim[2*num + 1])
	       ||(phase < sp_bas->mask.lim[2*num])) break;
	  }
          if(num == sp_bas->mask.num) sd_bas->numSD[0]++;

	}  // end check for particle limitations

      }
      else if(sd_bas->J_type && (config == time_reversed_config))  {
	if(sp_bas->mask.num ==  0) {
	  sd_bas->numSD[1]++;
	}
	else  {     // check for particle limitations
	  for(num = 0; num < sp_bas->mask.num; num++) {
	    high = config & sp_bas->mask.list[num];
	    phase = 0;
	    while(high) {
	      high &= high - 1;
	      phase++;
	    }
	    if(  (phase > sp_bas->mask.lim[2*num + 1])
	       ||(phase < sp_bas->mask.lim[2*num])) break;
	  }
          if(num ==  sp_bas->mask.num) sd_bas->numSD[1]++;

	} // end check for particle limitations
      }
    } /* end parity and angular momentum test */

	     /* new particle configuration */

    for(loop = 0; (loop < sd_bas->part) && ((matr[loop]+1) >= matr[loop+1]); loop++);
    matr[loop]++;
    for(loop--; loop >= 0; loop--) matr[loop] = loop;

  }while(matr[sd_bas->part] ==  sp_bas->numm_orb);  
       /* end running through all configurations */

  free(matr);  /* release temporary memory */

  return;

} /* End: function number_of_trsym_SD() */

  /*
  ** The function                                           
  **          number_of_nosym_SD()                   
  ** calculates the number of slater determinants |SD> for n particles
  ** with given parity and MJ = m_value.
  */

static void number_of_nosym_SD(SP_BAS *sp_bas, SD_BAS *sd_bas)
{
  char
          *func = {"number_of_nosym_SD():"};
  int
          loop, m, par, *matr, num, phase;
  ULL
          config, high;

  matr = MALLOC(sd_bas->part + 1, int, func, "matr[]");

	      /* the lowest configuration */

  for(loop = 0; loop < sd_bas->part; loop++) matr[loop] = loop;
  matr[loop] = sp_bas->numm_orb;

       /* calculate the total number of SD */

  sd_bas->numSD[0] = 0; // initialization 
  sd_bas->numSD[1] = 0;

  do  {   /* run through all possible configurations */
    for(m = 0, par = +1, loop = 0; loop < sd_bas->part; loop++)  { /* m-value and parity */
      m   += sp_bas->mbas[matr[loop]].m;
      par *=  sp_bas->mbas[matr[loop]].par;
    }
    if((m == sd_bas->MJ) && (par == sd_bas->P)) {

      if(sp_bas->mask.num == 0) { 
	sd_bas->numSD[0]++;
      }
      else  {   // check for particle limitations
  	        // code the configuration
	
	for(config = ULL_ZERO, loop = 0; loop < sd_bas->part; loop++)  {
	  config += ULL_ONE << matr[loop];
	}
	for(num = 0; num < sp_bas->mask.num; num++) {
	  high = config & sp_bas->mask.list[num];
	  phase = 0;
	  while(high) {
	    high &= high - 1;
	    phase++;
	  }
	  if(  (phase < sp_bas->mask.lim[2*num])
	     ||(phase > sp_bas->mask.lim[2*num + 1])) break;
	}
	if(num == sp_bas->mask.num) sd_bas->numSD[0]++;

      }  // end check for particle limitations

    }
	     /* new particle configuration */

    for(loop = 0; (loop < sd_bas->part) && ((matr[loop]+1) >= matr[loop+1]); loop++);
    matr[loop]++;
    for(loop--; loop >= 0; loop--) matr[loop] = loop;
  } while(matr[sd_bas->part] == sp_bas->numm_orb);  /* end running through all configurations */

  free(matr);  /* release temporary memory */
  
} /* End: function number_of_nosym_SD() */

  /*
  ** The function                                           
  **                 tr_sym_SD_config()                   
  ** calculates and stores the complete set of slater determinants |SD()> for n
  ** particles with MJ = 0. It calculates the asymmetric |SD()> and the symmetric 
  ** |SD()> and stored in sd_bas->SD[], asymSD[]> first, The number of asymmetric 
  ** |SD()> is checked against sd_bas->numSD[0] and the number of symmetric |SD()> 
  ** against sd_bas->numSD[1].
  ** If total angular momentum J is odd, no sym |SD> is stored:
  **             J odd -> sd_bas->J_type = FALSE. J even -> sd_bas->J_type = TRUE
  */

void tr_sym_SD_config(SP_BAS *sp_bas, SD_BAS *sd_bas)
{
   char
            *func = {"tr_sym_SD_config(): "};
   int
            loop, m, par, last_orbit, count_asym, 
            count_sym, *matr, num, phase;
   ULL
            *ptr_asym, *ptr_sym, config, tr_config, high; 

   last_orbit =  sp_bas->numm_orb - 1;
   count_asym = 0;
   count_sym  = 0;
   ptr_asym   = sd_bas->SD;
   ptr_sym    = ptr_asym + sd_bas->numSD[0];

   matr = MALLOC(sd_bas->part + 1, int, func, "matr[]");

      /* the lowest configuration */

   for(loop = 0; loop < sd_bas->part; loop++) matr[loop] = loop;
   matr[loop] = sp_bas->numm_orb;

   do  {                                          /* run through all configuration */
            /* find m-value and parity */

      for(m = 0, par = +1, loop = 0; loop < sd_bas->part; loop++)  {
         m   += sp_bas->mbas[matr[loop]].m;
         par *=  sp_bas->mbas[matr[loop]].par;
      }
      if((m == 0) && (par == sd_bas->P))  {

        /* Code the configuration  and its time-reversed */

         config    = 0;
         tr_config = 0;
         for(loop = 0; loop < sd_bas->part; loop++)  {
            config    += ULL_ONE << matr[loop];
            tr_config += ULL_ONE <<(last_orbit - matr[loop]);
         }
         if(config < tr_config)  {

	   if(sp_bas->mask.num == 0) {
	     *ptr_asym = config;
	     ptr_asym++;
	     count_asym++;
	   }
	   else  {       // check for particle limitations
	     for(num = 0; num < sp_bas->mask.num; num++) {
	       high = config & sp_bas->mask.list[num];
	       phase = 0;
	       while(high) {
		 high &= high - 1;
		 phase++;
	       }
	       if(  (phase > sp_bas->mask.lim[2*num + 1])
		    ||(phase < sp_bas->mask.lim[2*num])) break;
	     }
	     if(num == sp_bas->mask.num) {  // |config> included
	       *ptr_asym = config;
	       ptr_asym++;
	       count_asym++;
	     }
	   }  // end check for particle limitations

	 } // end config < tr_config

         else if(sd_bas->J_type && (config == tr_config))  {

	   if(sp_bas->mask.num == 0) {
	     *ptr_sym = config;
	     ptr_sym++;
	     count_sym++;
	   }
	   else  {    // check for particle limitations
	     for(num = 0; num < sp_bas->mask.num; num++) {
	       high = config & sp_bas->mask.list[num];
	       phase = 0;
	       while(high) {
		 high &= high - 1;
		 phase++;
	       }
	       if(  (phase > sp_bas->mask.lim[2*num + 1])
		  ||(phase < sp_bas->mask.lim[2*num])) break;
	     }
	     if(num == sp_bas->mask.num) { // |config> included
	       *ptr_sym = config;
	       ptr_sym++;
	       count_sym++;
	     }
	   }  // end check for particle limitations

	 } // end config < tr_config

      } /* end parity and angular momentum test */

                 /* new particle configuration */

      for(loop = 0; (loop < sd_bas->part) && ((matr[loop]+1) >= matr[loop+1]); loop++);
      matr[loop]++;
      for(loop--; loop >= 0; loop--) matr[loop] = loop;
   } while(matr[sd_bas->part] == sp_bas->numm_orb);  /* end of particle loop */

   free(matr);  /* release temporary memory */

          /* Check the number of identical SD */

   if((count_asym != sd_bas->numSD[0]) || (count_sym != sd_bas->numSD[1]))   {
      printf("\n\nError in function tr_sym_SD_config():");
      printf("\nNumber of symmetric SD = %d, should be %d", count_sym, sd_bas->numSD[1]);
      printf("\nNumber of asymmetric SD = %d, should be %d\n", count_asym, sd_bas->numSD[0]);
      MPI_Abort(MPI_COMM_WORLD,Rank);
   }  
} /* End: function tr_sym_SD_config() */

  /*
  ** The function                                           
  **               no_sym_SD_config()                   
  ** calculates and stores the complete set of slater determinants |SD()> 
  ** for n identical particles with angular projection MJ and fixed parity (parity).
  ** The resulting |SD()> are stored in sd_bas->SD[0]. The number of no-symmetric
  ** SD is checked against sd_bas->numSD[0].
  */

void no_sym_SD_config(SP_BAS *sp_bas, SD_BAS *sd_bas)
{
   char
           *func = {"no_sym_SD_config(): "};
   int
           loop, m, par, count, *matr, num, phase;
   ULL
           *ptr_SD, config, high;

   matr = MALLOC(sd_bas->part + 1, int, func, "matr[]");

      /* the lowest neutron configuration */

   for(loop = 0; loop < sd_bas->part; loop++) matr[loop] = loop;
   matr[loop] =  sp_bas->numm_orb;

        /* calculate the total number of SD.*/

   count  = 0;   /* initialization */
   ptr_SD = sd_bas->SD;

        /* run through all configurations */

   do  { 
        /* find m-value and parity */

      for(m = 0, par = +1, loop = 0; loop < sd_bas->part; loop++)  {
        m   += sp_bas->mbas[matr[loop]].m;
       par *=  sp_bas->mbas[matr[loop]].par;
      }
      if((m == sd_bas->MJ) && (par == sd_bas->P))  {

        /* Code the configuration */

         for(config = ULL_ZERO, loop = 0; loop < sd_bas->part; loop++)  {
            config += ULL_ONE << matr[loop];
         }

	 if(sp_bas->mask.num == 0) { 
	   *ptr_SD = config;
	   ptr_SD++;
	   count++;
	 }
	 else  {   // check for particle limitations
	   for(num = 0; num < sp_bas->mask.num; num++) {
	     high = config & sp_bas->mask.list[num];
	     phase = 0;
	     while(high) {
	       high &= high - 1;
	       phase++;
	     }
	     if(  (phase > sp_bas->mask.lim[2*num + 1])
		  ||(phase < sp_bas->mask.lim[2*num])) break;;
	   }
	   if(num == sp_bas->mask.num) { // |config> included
	     *ptr_SD = config;
	     ptr_SD++;
	     count++;
	   }
	 } // end check for particle limitations

      } /* end check m-value and parity */

            /* new particle configuration */

      for(loop = 0; (loop < sd_bas->part) && ((matr[loop]+1) >= matr[loop+1]); loop++);
      matr[loop]++;
      for(loop--; loop >= 0; loop--)  {
         matr[loop] = loop;
      }
   } while(matr[sd_bas->part] == sp_bas->numm_orb);  /* end running through all configurations */

       /* Check the number of identical SD */

   if(sd_bas->numSD[0] != count)  {
      printf("\n\nError in function no_sym_SD_config():");
      printf("\nIn calculating the particle SD with no symmetry");
      printf("\nNumber of SD = %d from function number_of_nosym_SD()",sd_bas->numSD[0]);
      printf("\nThe current calculation gives num = %d\n", count);
      MPI_Abort(MPI_COMM_WORLD,Rank);
   }
   free(matr);  /* release temporary memory */

}   /* End: function no_sym_SD_config() */
