/*******************  The module lanc-pn-SD-basis.c  ******************/

#include "shell.h"

   /*
   ** The entrance function
   **  int pn_slater_determinant(PN_SHELL *model, GR_BASIS **gr_basis, GROUP **group)
   ** allocates memory, calculates and stores the complete set
   ** of slater determinants for a system of Z protons
   ** |SD(Z)M_Z P_Z> and N neutrons |SD(N)M_N P_N>
   ** If no slater determinants are found, the function
   ** returns FALSE, otherwise TRUE.
   */

             /**** local data definitions ****/

        /* data structure max and min |SD(Z)MZ> */

   typedef  struct {
      int      max_MZ[2],    /* parity plus(minus) = [0]([1]) */
               min_MZ[2];
      } M_LIMITS;

	    /* local function declarations */

   static int pn_m_value_limits(M_LIMITS *m_limits, PN_SHELL *model, GR_BASIS *gr_basis);
      /*
      ** calculates and return in model->max_MZ[] and model->min_MZ[]
      ** maximum and minimum m-values for proton |SD(Z)M P> with fixed
      ** parity. The calculation is based on the condition that the
      ** total proton/neutron basis states has the form
      **                   |SD(Z)M_Z P_Z,SD(N) M_N P_N: M P>
      ** with condition M = m_Z + M_N and P = P_Z * P_N
      ** If these conditions are not fullfilled, the function
      ** returns FALSE, otherwise TRUE.
      */

   static int pn_max_m_value(int type, int parity, GR_BASIS const *gr_basis);
      /*
      ** calculates for a given parity maximum M-values of a set
      ** |SD(type)M P> for a given number of identical particles
      ** in a single-particle orbits specified in mbas[]. 
      ** The function returns maximum M-value.
      */

   static int pn_modify_m_values(M_LIMITS *m_limits, GR_BASIS *gr_basis, int type);
      /*
      ** takes maximum and minimum m-values for |SD(Z)M P> and |SD(N)M P> for 
      ** a given proton parity  and modify the limits according to the condition 
      **                 tot_MJ = MZ + MN 
      ** The modified limits for protons are stored in model->max_MZ[type]
      ** and model->min_MZ[type].
      ** Note: All m_values are twice physical values.
      ** if (max_mz < min_mz) or (max_n < min_n) the function returns
      ** FALSE, otherwise TRUE.
      */

   static void pn_groups(M_LIMITS *m_limits, GR_BASIS *gr_basis, GROUP **group);
      /*
      ** calculates the complete set of data for all groups of a 
      ** given shell model system. Data are stored in group[].
      */
   static int pn_number_of_SD(int type, GR_BASIS const *gr_basis, GROUP const *group);

      /*
      ** calculates the number of slater determinants for a given set of identical
      ** particles (num_part) with fixed total m-value (m_value) and parity (parity).
      ** Possible limitations on number of particles in j-orbits are checked through
      ** par_lim[]. No time-reversal symmetry.
      */

   static void pn_particle_SD_config(int type, GR_BASIS const *gr_basis, GROUP const *group);

      /*
      ** calculates all slater determinants |SD()MP> for a given number of identical
      ** particle (num_part) with fixed total m-value (m_value) and parity (parity).
      ** No time-reversal symmetry.
      */

             /**** End: function declarations ****/


               /**** The function definitions  ****/ 

   /*
   ** The entrance function                                                      
   **            int pn_slater_determinant() 
   ** allocates memory, calculates and stores the complete set
   ** of slater determinants for a system of Z protons
   ** |SD(Z)M_Z P_Z> and N neutrons |SD(N)M_N P_N>
   ** If no slater determinants are found, the function
   ** returns FALSE, otherwise TRUE.
   */

int pn_slater_determinant(PN_SHELL *model, GR_BASIS *gr_basis, GROUP **group)
{
  char         *func = {"pn_slater_determinant(): "};
  int          loop, start_amp;
  M_LIMITS     m_limits;
  GROUP        *gr_ptr;

            /* max and min MZ_values for both plus and minus parity */ 

  if(pn_m_value_limits(&m_limits, model, gr_basis) == FALSE)  return FALSE;

               /* number of groups*/

  if(gr_basis->parZ == 0) {             /* both parities: [0] = +  / [1] = - */

    gr_basis->num_gr =  MAX(m_limits.max_MZ[0], m_limits.max_MZ[1]) 
              - MIN(m_limits.min_MZ[0], m_limits.min_MZ[1]) + 2;
  }
  else {                              /* only one parity */
    gr_basis->num_gr = (m_limits.max_MZ[0] - m_limits.min_MZ[0])/2 +1;
  }
               /* memory for all group data */

  *group = MALLOC(gr_basis->num_gr, GROUP, func, " group[]");

  pn_groups(&m_limits, gr_basis, group);

  gr_basis->tot_numSD_Z = 0;        /* initialization basic output data */
  gr_basis->tot_numSD_N = 0;
  gr_basis->tot_numSD   = 0;
  gr_basis->maxSD[0]    = 0;
  gr_basis->maxSD[1]    = 0;

  gr_ptr    = *group;
  start_amp = 0;                     /* start position for first group */
  for(loop = 0; loop < gr_basis->num_gr; loop++, gr_ptr++) {
    gr_basis->tot_numSD_Z += gr_ptr->numSD[0];
    gr_basis->tot_numSD_N += gr_ptr->numSD[1];
    gr_ptr->start_amp      = start_amp;
    gr_basis->tot_numSD   += gr_ptr->numSD[0] * gr_ptr->numSD[1];
    gr_basis->maxSD[0]     = MAX(gr_basis->maxSD[0], gr_ptr->numSD[0]); 
    gr_basis->maxSD[1]     = MAX(gr_basis->maxSD[1], gr_ptr->numSD[1]); 

    start_amp =  gr_basis->tot_numSD;  /* start position for next group */
 }
  return TRUE;
} /* End: function pn_slater_determinant() */

     /*
     ** The function 
     **           pn_m_value_limits()
     ** calculates and return in m_limits->max_MZ[] and m_limits->min_MZ[]
     ** maximum and minimum m-values for proton |SD(Z)M P> with fixed
     ** parity. The calculation is based on the condition that the
     ** total proton/neutron basis states has the form
     **                   |SD(Z)M_Z P_Z,SD(N) M_N P_N: M P>
     ** with condition M = M_Z + M_N and P = P_Z * P_N
     ** If these conditions are not fullfilled, the function
     ** returns FALSE, otherwise TRUE.
     */

static int pn_m_value_limits(M_LIMITS *m_limits, PN_SHELL *model, GR_BASIS *gr_basis)
{
  int
          plus, minus;

  plus  = TRUE;                                           /* initialization */
  minus = TRUE;

       /* Identification: plus parity = [0], minus parity = [1] */

  if(gr_basis->parZ == 0) { 

            /* proton parity + */

    m_limits->max_MZ[PLUS] = pn_max_m_value(PROTON, +1, gr_basis);
    m_limits->min_MZ[PLUS] = - m_limits->max_MZ[PLUS];

          /* possible modification of proton parity + due to neutrons  */

    plus = pn_modify_m_values(m_limits, gr_basis, PLUS);

                   /* proton parity - */

    m_limits->max_MZ[MINUS] = pn_max_m_value(PROTON, -1, gr_basis);
    m_limits->min_MZ[MINUS] = - m_limits->max_MZ[MINUS];

          /* possible modification of proton parity - due to neutrons */

    minus = pn_modify_m_values(m_limits, gr_basis, MINUS);

    if((plus == TRUE) && (minus == TRUE)) {
      return TRUE;
    }
    else if((plus == TRUE) && (minus == FALSE)) {
      gr_basis->parZ = +1;
      m_limits->max_MZ[1] = m_limits->max_MZ[0];          /* same values */
      m_limits->min_MZ[1] = m_limits->min_MZ[0];
      return TRUE;
    }
    else if((plus == FALSE) && (minus == TRUE)) {
      gr_basis->parZ = -1;
      m_limits->max_MZ[0] = m_limits->max_MZ[1];          /* same values */
      m_limits->min_MZ[0] = m_limits->min_MZ[1];
      return TRUE;
    }
   else if((plus == FALSE) && (minus == FALSE)) {
     return FALSE;                                     /* no basis states */
    }
         /* should never reach this point */

    printf("\n\nError in function  pn_m_value_limits():");
    printf("\nWrong return values from function (pn_modify_m_values()");
    printf("\nplus = %d  minus = %d -- should be 0 and/or +1\n", plus, minus);
    exit(1);
  }   /* end gr_basis->parZ = 0 */

  else if(abs(gr_basis->parZ) == 1) { 

    m_limits->max_MZ[0] = pn_max_m_value(PROTON, gr_basis->parZ, gr_basis);
    m_limits->min_MZ[0] = - m_limits->max_MZ[0];

    if(pn_modify_m_values(m_limits, gr_basis,
                ((gr_basis->parZ == +1) ? PLUS : MINUS)) == TRUE) {
      m_limits->max_MZ[1] = m_limits->max_MZ[0];          /* same values */
      m_limits->min_MZ[1] = m_limits->min_MZ[0];
      return TRUE;
    }
    else {
      return FALSE;
    }
  } /* end gr_basis->parZ = +-1 */

  else  {                             /* should never reach this point */
    printf("\n\nError in function  pn_m_value_limits():");
    printf("\nWrong value for model.parZ = %d -- should be 0 and/or +-1\n",
                                                               gr_basis->parZ);
    exit(1);
  }
  return TRUE;

} /* End: function pn_m_value_limits() */ 

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
      **          pn_modify_m_values()
      ** takes maximum and minimum m-values for |SD(Z)M P> and |SD(N)M P> for 
      ** a given proton parity  and modify the limits according to the condition 
      **                 tot_MJ = MZ + MN 
      ** The modified limits for protons are stored in m_limits->max_MZ[type]
      ** and m_limits->min_MZ[type].
      ** Note: All m_values are twice physical values.
      ** if (max_mz < min_mz) or (max_n < min_n) the function returns
      ** FALSE, otherwise TRUE.
      */

static int pn_modify_m_values(M_LIMITS *m_limits, GR_BASIS *gr_basis, int parity)
{
   int    max_mz, min_mz, tot_mj, max_mn, min_mn, mz, mn, step, num;

   /* Parity identification: proton-parity plus (minus)--> parity = 0, (1) */

   max_mz = m_limits->max_MZ[parity];            /* local variable */
   min_mz = m_limits->min_MZ[parity]; 
   tot_mj = gr_basis->MJ;

       /* initial maximum and minimum neutron m-values */

   max_mn = pn_max_m_value(NEUTRON, ((parity == 0) ? gr_basis->P : -gr_basis->P),
                                                                        gr_basis);
   min_mn = -max_mn;

	  /* find new max_mz and min_mn */
 
   if((step = max_mz + min_mn - tot_mj) > 0)  {
      max_mz -= step;                          /* reduce max_MZ */
   }
   else if(step < 0) {
      min_mn -= step;                          /* increase min_MN */
   }
     /* number of steps from max_mz to min_mz */

   for(num = 1, mz = max_mz, mn = min_mn;
           (mz >min_mz) && (mn < max_mn); mz -= 2, mn += 2, num++); 

   min_mz = max_mz - 2 * (num - 1);
   m_limits->max_MZ[parity] = max_mz;        /* store and return the result */
   m_limits->min_MZ[parity] = min_mz;

   if((max_mz < min_mz) || (max_mn < min_mn)) return FALSE;

   return TRUE; 

} /* End: function pn_modify_m_values() */

      /*
      ** The function
      **           pn_groups()
      ** calculates the complete set of data for all groups of a 
      ** given shell model system. Data are stored in group[].
      */

static void pn_groups(M_LIMITS *m_limits,  GR_BASIS *gr_basis,GROUP **group)
{
   char     *func = {"pn_groups(): "};
   int      m, max_mz, min_mz, numSD, k;
   GROUP    *gr_ptr, *init_gr_ptr;

   numSD  = 0;                          /* number of basis vectors: initialization */

   max_mz = MAX(m_limits->max_MZ[PLUS], m_limits->max_MZ[MINUS]);
   min_mz = MIN(m_limits->min_MZ[PLUS], m_limits->min_MZ[MINUS]);

        /* Identification: proton = [0] and neutron = [1] */
 
   for(m = max_mz, gr_ptr = *group; m >= min_mz; m -= 2) {

        /* first group for given m-value */
      
      gr_ptr->m[0]              = m;              
      gr_ptr->m[1]              = gr_basis->MJ - m;
      gr_ptr->par[0]            = (gr_basis->parZ == 0) ? +1 : gr_basis->parZ;
      gr_ptr->par[1]            = gr_basis->P * gr_ptr->par[0];

         /* calculate and check number of |SD(Z)> and |SD(N)> */

      gr_ptr->numSD[0] = pn_number_of_SD(PROTON, gr_basis, gr_ptr);
      gr_ptr->numSD[1] = pn_number_of_SD(NEUTRON, gr_basis, gr_ptr);

      if((gr_ptr->numSD[0] * gr_ptr->numSD[1])== 0) {
         gr_ptr->numSD[0] = 0;
         gr_ptr->numSD[1] = 0;
         gr_ptr->SD[0]    = NULL_PTR;
         gr_ptr->SD[1]    = NULL_PTR;
      }
      else {
         numSD        += gr_ptr->numSD[0] * gr_ptr->numSD[1];
         gr_ptr->SD[0] = MALLOC(gr_ptr->numSD[0] + 1, UL, func,"gr_ptr.SD[0][]");
         pn_particle_SD_config(PROTON, gr_basis, gr_ptr);

         gr_ptr->SD[1] = MALLOC(gr_ptr->numSD[1] + 1, UL, func,"gr_ptr.SD[1][]");
         pn_particle_SD_config(NEUTRON, gr_basis, gr_ptr);
      }
      gr_ptr++;                                    /* point to next group[] */

           /* second  group for given m-value if both parities is allowed */

     if(gr_basis->parZ != 0) continue;
      
      gr_ptr->m[0]              = m;
      gr_ptr->m[1]              = gr_basis->MJ - m;
      gr_ptr->par[0]            = -1;
      gr_ptr->par[1]            = gr_basis->P * gr_ptr->par[0];

         /* calculate and check number of |SD(Z)> and 1SD/N)> */

      gr_ptr->numSD[0] = pn_number_of_SD(PROTON, gr_basis, gr_ptr);
      gr_ptr->numSD[1] = pn_number_of_SD(NEUTRON, gr_basis, gr_ptr);

      if((gr_ptr->numSD[0] * gr_ptr->numSD[1]) == 0) {
         gr_ptr->numSD[0] = 0;
         gr_ptr->numSD[1] = 0;
         gr_ptr->SD[0]    = NULL_PTR;
         gr_ptr->SD[1]    = NULL_PTR;
      }
      else {
         numSD             += gr_ptr->numSD[0] * gr_ptr->numSD[1];
         gr_ptr->SD[0] = MALLOC(gr_ptr->numSD[0] + 1, UL, func,"gr_ptr.SD[0][]");
         pn_particle_SD_config(PROTON, gr_basis, gr_ptr);
         gr_ptr->SD[1] = MALLOC(gr_ptr->numSD[1] + 1, UL, func,"gr_ptr.SD[1][]");
         pn_particle_SD_config(NEUTRON, gr_basis, gr_ptr);
      }
      gr_ptr++;        /* point to next group[] */

   } /* end loop through all groups */

        /*
        ** calculate final group number for nondiag
	** matrix elemenst within each group
        */

   init_gr_ptr = *group;
   max_mz      = gr_basis->num_gr - 1;
   for(m = 0; m < max_mz; m++, init_gr_ptr++) {
      gr_ptr = init_gr_ptr + 1;
      init_gr_ptr->max_gr = m;
      for(k = m + 1; k <= max_mz; k++, gr_ptr++) {
         if(gr_basis->mbas[0]->m < (init_gr_ptr->m[0] - gr_ptr->m[0])) break;
         init_gr_ptr->max_gr++;
      } /* end k-loop */
   } /* end m-loop */

} /* End: function pn_groups() */

     /*
     ** The function 
     **           pn_number_of_SD()
     ** calculates the number of slater determinants for a given set of identical
     ** particles (num_part) with fixed total m-value (m_value) and parity (parity).
     ** Possible limitations on number of particles in j-orbits are checked through
     ** par_lim[]. No time-reversal symmetry.
     */

static int pn_number_of_SD(int type, GR_BASIS const *gr_basis, GROUP const *group)
{
   char         *func = {"pn_number_of_SD(): "};
   register int loop, number, m, par, m_value, parity, num_part, num_orb, num_mask;
   int          *matr, *ptr;
   UL           sd;
   MBAS        *mbas;
   PART_LIM     part_lim;
   MASK         *table_ptr;
   
   m_value  = group->m[type];                         /* initialization */
   parity   = group->par[type];
   num_part = gr_basis->part[type];
   num_orb  = gr_basis->m_orb[type];
   mbas     = gr_basis->mbas[type];
   part_lim = gr_basis->part_lim[type];

   if(num_part == 0) return 1;                             /* no particles */

   matr     = MALLOC(num_part + 1, int, func, "matr[]");     /* temporary memory */
   num_mask = part_lim.num;
   for(loop = 0; loop < num_part; loop++) matr[loop] = loop; /* lowest config. */
   matr[num_part] = num_orb;

   number = 0;                                          /* initialization */
   do  {   /* particle loop */
      sd = UL_ZERO;
      for(m = 0, par = +1, ptr = matr, loop = 0; loop < num_part; loop++) {
         m   += mbas[*ptr].m; 
         par *= mbas[*ptr].par; 
         sd  += UL_ONE << *(ptr++); 
      }
      if((m == m_value) && (par == parity))   {
         if(num_mask == 0) {
            number++;
         }
         else {
            table_ptr = part_lim.table;
	    for(loop = 0; loop < num_mask; loop++, table_ptr++) {
               m = sd & table_ptr->mask;
               par = 0;
               while(m) {
                  m &= m-1;
                  par++;
               }
               if((par < table_ptr->min) || (par > table_ptr->max)) break;
            }
            if(loop == num_mask) number++;
         } 
      }
                /* new particle configuration */

      for(loop = 0; (loop < num_part) && ((matr[loop]+1) >= matr[loop+1]); loop++);
      matr[loop]++;
      for(loop--; loop >= 0; loop--)   matr[loop] = loop;
   } while(matr[num_part] == num_orb);  /* end of particle loop */

   free(matr);  /* release temporary memory */

   return number;

} /* End: function pn_number_of_SD() */

     /*
     ** The function 
     **         pn_particle_SD_config()
     ** calculates all slater determinants |SD()MP> for a given number of identical
     ** particle (num_part) with fixed total m-value (m_value) and parity (parity).
     ** No time-reversal symmetry.
     */

static void pn_particle_SD_config(int type, GR_BASIS const *gr_basis, GROUP const *group)
{
   char          *func = {"pn_particle_SD_config(): "};
   register int  loop, number, m, par, m_value, parity, num_part, num_orb, numSD;
   int           *matr, *ptr;
   UL            *sd;
   MBAS          *mbas;
   PART_LIM      part_lim;
   MASK          *table_ptr;

   m_value  = group->m[type];                       /* initialization */
   parity   = group->par[type];
   num_part = gr_basis->part[type];
   num_orb  = gr_basis->m_orb[type];
   mbas     = gr_basis->mbas[type];
   part_lim = gr_basis->part_lim[type];
   numSD    = group->numSD[type];
   sd       = group->SD[type];

   if(!num_part) {
     if(numSD != 1) {
        printf("\n\nError in function pn_particle_SD_config():");
        printf("\nIn calculating particle type = %d  |SD> with no symmetry", type);
        printf("\nParticle number = ZERO, number of |SD> = %d\n", numSD) ;
        exit(1);
     }
      *sd = UL_ONE;                /* |SD> = |vacuum> */
      return;
   }

   matr     = MALLOC(num_part + 1, int, func, "matr[]");     /* temporary memory */
   for(loop = 0; loop < num_part; loop++) matr[loop] = loop; /* lowest config. */
   matr[num_part] = num_orb;

   number = 0;                                            /* check number of SD */
   do  {                                                  /* particle loop */
      *sd = UL_ZERO;                                   /* initialization of |SD> */
      for(m = 0, par = +1, ptr = matr,loop = 0; loop < num_part; loop++) {
         m   += mbas[*ptr].m; 
         par *= mbas[*ptr].par; 
         *sd += UL_ONE << *(ptr++); 
      }
      if((m == m_value) && (par == parity))  {
         if(part_lim.num == 0) {
            sd++;
            number++;
         }
         else {
            table_ptr = part_lim.table;
	    for(loop = 0; loop < part_lim.num; loop++, table_ptr++) {
               m = *sd & table_ptr->mask;
               par = 0;
               while(m) {
                  m &= m-1;
                  par++;
               }
               if((par < table_ptr->min) || (par > table_ptr->max)) break;
            }
            if(loop == part_lim.num)  {
               sd++;
               number++;
            }
         }

      } /* end check m-value and parity */

              /* new particle configuration */

      for(loop = 0; (loop < num_part) && ((matr[loop]+1) >= matr[loop+1]); loop++);
      matr[loop]++;
      for(loop--; loop >= 0; loop--)  matr[loop] = loop;

   } while(matr[num_part] == num_orb);  /* end of particle loop */

   if(numSD != number)  {                   /* Check the number of identical SD */
      printf("\n\nError in function pn_particle_SD_config():");
      printf("\nIn calculating particle |SD> with no symmetry");
      printf("\nNumber of type = %d |SD> = %d",type, numSD);
      printf(" - current calculation gives num = %d\n", number);
      exit(1);
   }
   free(matr);  /* release temporary memory */

} /* End: function pn_particle_SD_config() */
