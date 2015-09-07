
/*******************  The module lanc-pn-SD-basis.c  ******************/

#include  "PAR-pnShellModel.h"


   /*
   ** The entrance function
   **  int pn_slater_determinant( GR_BAS **grBas, GROUP**group)
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

static void  oscillatorRefState(GR_BAS *grBas);
    /*
    ** calculates the refOscN  i.e summed N for 
    ** the lowest possible state with 
    ** required parity grBas->P
    */

static int sort_oscOrb(const JBAS *one, const JBAS *two);
    /*
    ** is a utility function for the library function qsort() in order to
    ** sort single-part elements of type SP_BAS tempJbas[]
    ** after increasing tempJbas[].N
    */

static int pn_m_value_limits(M_LIMITS *m_limits, GR_BAS *grBas);
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

   static int pn_max_m_value(int type, int parity, SP_BAS const *spBas);
      /*
      ** calculates for a given parity maximum M-values of a set
      ** |SD(type)M P> for a given number of identical particles
      ** in a single-particle orbits specified in mbas[]. 
      ** The function returns maximum M-value.
      */

   static int pn_modify_m_values(M_LIMITS *m_limits, GR_BAS *grBas, int type);
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

static int calcGroupData(GR_BAS *grBas,GROUP *grPtr);
        /*
	** calculates and store in group[]
	** all specific group data
	** ig current group does not contain 
	** any basis states |SD(Z),SD(N)>
	** function returns FALSE, otherwise TRUE
       	*/

static int pn_number_of_SD(int type, GR_BAS const *grBas, GROUP const *group);

      /*
      ** calculates the number of slater determinants for a given set of identical
      ** particles (num_part) with fixed total m-value (m_value) and parity (parity).
      ** Possible limitations on number of particles in j-orbits are checked through
      ** par_lim[]. No time-reversal symmetry.
      */

static void pn_particle_SD_config(int type, GR_BAS const *grBas, GROUP const *group);

      /*
      ** calculates all slater determinants |SD()MP> for a given number of identical
      ** particle (num_part) with fixed total m-value (m_value) and parity (parity).
      ** No time-reversal symmetry.
      */

static void pnGlobalData(GR_BAS *grBas, GROUP *group);
      /*
      ** calculates and store in grBas
      ** data common for all groups
      */

static void occ_orbits(int type, GR_BAS *grBas,GROUP *group, int *occ);
    /* 
    ** calculates and stores orbit numbers of all occupied
    ** orbits  in all |SD(type)> in present group 
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

int pn_slater_determinant(GR_BAS *grBas,GROUP **group)
{
  char         *func = {"pn_slater_determinant(): "};
  int          max_mz, min_mz, m, new_gr_num;
  M_LIMITS     mLim;
  GROUP        *grPtr, **old_group;

  // calculate the refernce oscillator energy

  oscillatorRefState(grBas);

  // energy excitation limits for basis |SD(Z);SD(N)>

  grBas->maxOscN = grBas->refOscN + grBas->delN;
 
  // max and min MZ_values for both plus and minus parity

  if(pn_m_value_limits(&mLim, grBas) == FALSE)  return FALSE;

  // number of groups

  if(grBas->parZ == 0) {  // both parities: [0] = +  / [1] = -

    grBas->num_gr =  MAX(mLim.max_MZ[0], mLim.max_MZ[1]) 
              - MIN(mLim.min_MZ[0], mLim.min_MZ[1]) + 2;
  }
  else { // only one parity
    grBas->num_gr = (mLim.max_MZ[0] - mLim.min_MZ[0])/2 +1;
  }

  // memory for all group data

   *group = MALLOC(grBas->num_gr, GROUP, func, " group[]");

  // group limitation in mz-values
 
   max_mz = MAX(mLim.max_MZ[PLUS], mLim.max_MZ[MINUS]);
   min_mz = MIN(mLim.min_MZ[PLUS], mLim.min_MZ[MINUS]);

   // Identification: proton = [0] and neutron = [1]

   new_gr_num = 0; // Only groups with basis states are included

   grPtr = *group;
   grBas->maxGroupNumSD_Z = 0;   // initialization
   grBas->maxGroupNumSD_N = 0;

   for(m = max_mz; m >= min_mz; m -= 2) {
     grPtr->m[0] = m;              
     grPtr->m[1] = grBas->MJ - m;

     // parity for first group for given m-value

     grPtr->par[0]  = (grBas->parZ == 0) ? +1 : grBas->parZ;
     grPtr->par[1]  = grBas->P * grPtr->par[0];

     if(calcGroupData(grBas,grPtr) == TRUE) {
       grBas->maxGroupNumSD_Z = MAX(grBas->maxGroupNumSD_Z, grPtr->numSD[0]);
       grBas->maxGroupNumSD_N = MAX(grBas->maxGroupNumSD_N, grPtr->numSD[1]);
       grPtr++;
       new_gr_num++; // group OK
     }
     
     //  Next group with change of parity

     if(grBas->parZ != 0) continue;

     // m-values  for second group

     grPtr->m[0] = m;              
     grPtr->m[1] = grBas->MJ - m;

     // parity for second group

     grPtr->par[0]  = -1;
     grPtr->par[1]  = grBas->P * grPtr->par[0];

     if(calcGroupData(grBas,grPtr) == TRUE) {
       grBas->maxGroupNumSD_Z = MAX(grBas->maxGroupNumSD_Z, grPtr->numSD[0]);
       grBas->maxGroupNumSD_N = MAX(grBas->maxGroupNumSD_N, grPtr->numSD[1]);
       grPtr++;
       new_gr_num++; // group OK
     }
   } // end m-loop through all groups


   if(new_gr_num < grBas->num_gr) {   // Reduce number of groups 

     grBas->num_gr = new_gr_num;
     old_group = group;
     *group = REALLOC(*group, grBas->num_gr, GROUP, func,
                         "pn_slater__determinant()");
     if(group != old_group) { 
       printf("\n\nRank%d. Error in pn_slater__determinant()",Rank);
       printf("\nREALLOC() does not work properly - move memory loc");
       printf("\nof group[] to new place\n");
       MPI_Abort(MPI_COMM_WORLD,Rank);
     }
   }

   pnGlobalData(grBas,*group);

   return TRUE;
} // End: function pn_slater_determinant()

    /*
    ** The function 
    **   oscillatorRefState()
    ** calculates the refOscN  i.e summed N for 
    ** the lowest possible state with 
    ** required parity grBas->P
    */

static void  oscillatorRefState(GR_BAS *grBas)
{
  char   *func = {"oscillatorRefState()"};
  int    type, k, numj_orb, refOscN_Z, refOscN_N, totNumPart, numPart;
  JBAS   *tempJbas;


  refOscN_Z = 0;  // sum N for all proton particles
  refOscN_N = 0;  // sum N for all neutron particles

  // type = 0(proton), = 1(neutron)

  for(type = 0; type <= 1; type++)  { 
    numj_orb = grBas->spBas[type]->numj_orb;
    tempJbas = MALLOC(numj_orb,JBAS,func,"tempJbas");

    for(k = 0; k < numj_orb; k++) {
      tempJbas[k].N        = grBas->spBas[type]->jbas[k].N;
      tempJbas[k].osc      = grBas->spBas[type]->jbas[k].osc;
      tempJbas[k].l        = grBas->spBas[type]->jbas[k].l;
      tempJbas[k].j        = grBas->spBas[type]->jbas[k].j;
      tempJbas[k].min_part = grBas->spBas[type]->jbas[k].min_part;
      tempJbas[k].max_part = grBas->spBas[type]->jbas[k].max_part;
      tempJbas[k].e        = grBas->spBas[type]->jbas[k].e;
    }

    // sorting tenpJbas[] after increasing oscNum N

    qsort(tempJbas, (size_t)numj_orb ,(size_t)sizeof(JBAS),
	    (int(*)(const void *, const void *))sort_oscOrb);

    totNumPart = grBas->spBas[type]->part;

    for(k = 0; k < numj_orb; k++) {
      numPart = MIN(totNumPart, tempJbas[k].j + 1);
      if(type == 0) {
	refOscN_Z += numPart * tempJbas[k].N;
      }
      else {
	refOscN_N += numPart * tempJbas[k].N;
      }
      totNumPart -= numPart;
      if(totNumPart <= 0) break;
    }
    free(tempJbas);

  }  // end proton/neutron ref state


  // total reference oscillator configuration

  grBas->refOscN = refOscN_Z + refOscN_N;

  // Max proton oscillator excitation


  grBas->maxOscN_Z = refOscN_Z + grBas->delN_Z;

  // max neutron oscillator excitation 

  grBas->maxOscN_N =  refOscN_N + grBas->delN_N;

  // Max Oscillator N for total P/N config

  grBas->maxOscN = grBas->refOscN + grBas->delN;
 
} // end function  oscillatorRefState()

    /*
    ** The function                         
    **        int sort_oscOrb()                  
    ** is a utility function for the library function qsort() in order to
    ** sort single-part elements of type SP_BAS tempJbas[]
    ** after increasing tempJbas[].N
    */

static int sort_oscOrb(const JBAS *one, const JBAS *two)
{
  if(one->N > two->N)       return +1;
  else  if(one->N < two->N) return -1;
  else                      return  0;

} /* End: function sort_oscOrb() */

     /*
     ** The function 
     **           pn_m_value_limits()
     ** calculates and return in mLim->max_MZ[] and mLim->min_MZ[]
     ** maximum and minimum m-values for proton |SD(Z)M P> with fixed
     ** parity. The calculation is based on the condition that the
     ** total proton/neutron basis states has the form
     **                   |SD(Z)M_Z P_Z,SD(N) M_N P_N: M P>
     ** with condition M = M_Z + M_N and P = P_Z * P_N
     ** If these conditions are not fullfilled, the function
     ** returns FALSE, otherwise TRUE.
     */

static int pn_m_value_limits(M_LIMITS *mLim, GR_BAS *grBas)
{
  int    plus, minus;

  plus  = TRUE;  // initialization
  minus = TRUE;

  // Identification: + par = [0], - par = [1]

  if(grBas->parZ == 0) { 

    // proton parity plus

    mLim->max_MZ[PLUS] = pn_max_m_value(PROTON,+1, grBas->spBas[0]);
    mLim->min_MZ[PLUS] = - mLim->max_MZ[PLUS];

    // possible modification of proton parity + due to neutrons

    plus = pn_modify_m_values(mLim, grBas, PLUS);

    // proton parity  minus

    mLim->max_MZ[MINUS] = pn_max_m_value(PROTON,-1, grBas->spBas[0]);
    mLim->min_MZ[MINUS] = - mLim->max_MZ[MINUS];

    // possible modification of proton parity - due to neutrons

    minus = pn_modify_m_values(mLim, grBas, MINUS);

    if((plus == TRUE) && (minus == TRUE)) {
      return TRUE;
    }
    else if((plus == TRUE) && (minus == FALSE)) {
      grBas->parZ = +1;
      mLim->max_MZ[1] = mLim->max_MZ[0]; // same value
      mLim->min_MZ[1] = mLim->min_MZ[0];
      return TRUE;
    }
    else if((plus == FALSE) && (minus == TRUE)) {
      grBas->parZ = -1;
      mLim->max_MZ[0] = mLim->max_MZ[1]; // same value
      mLim->min_MZ[0] = mLim->min_MZ[1];
      return TRUE;
    }
    else if((plus == FALSE) && (minus == FALSE)) {
     return FALSE;  // no basis states */
    }
    // should never reach this point

    printf("\n\nRank%d. Error in function  pn_m_value_limits():",Rank);
    printf("\nWrong return values from function (pn_modify_m_values()");
    printf("\nplus = %d  minus = %d -- should be 0 and/or +1\n", plus, minus);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  } // end grBas->parZ = 0

  else if(abs(grBas->parZ) == 1) { 
    mLim->max_MZ[0] = pn_max_m_value(PROTON,grBas->parZ, grBas->spBas[0]);
    mLim->min_MZ[0] = - mLim->max_MZ[0];

    if(pn_modify_m_values(mLim, grBas,
                ((grBas->parZ == +1) ? PLUS : MINUS)) == TRUE) {
      mLim->max_MZ[1] = mLim->max_MZ[0]; // same values
      mLim->min_MZ[1] = mLim->min_MZ[0];
      return TRUE;
    }
    else {
      return FALSE;
    }
  } // end grBas->parZ = +-1

  else  {  // should never reach this point
    printf("\n\nRank%d. Error in function  pn_m_value_limits():",Rank);
    printf("\nWrong value for model.parZ = %d -- should be 0 and/or +-1\n",
                                                               grBas->parZ);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }
  return TRUE;

} // End: function pn_m_value_limits()

     /*
     ** The function 
     **           pn_max_m_value()
     ** calculates for a given parity maximum M-values of a set
     ** |SD(type)M P> for a given number of identical particles
     ** in a single-particle orbits specified in mbas[]. 
     ** The function returns maximum M-value.
     */

static int pn_max_m_value(int type,int parity, SP_BAS  const *spBas)
{ 
   char          *func = {"pn_max_m_value(): "};
   int           loop, m, par, num_part, num_orb;
   int           *matr, *ptr;
   MBAS          *mbas;

   num_part = spBas->part;   // initialization
   num_orb  = spBas->numm_orb;
   mbas     = spBas->mbas;

   matr = MALLOC(num_part + 1, int, func, "matr[]"); // temporary memory

   for(loop = 0; loop < num_part; loop++) matr[loop] = loop; // lowest config
   matr[num_part] = num_orb;

   do  {   // particle loop
      for(m = 0, par = +1, ptr = matr, loop = 0; loop < num_part; loop++) {
         m   += mbas[*ptr].m; 
         par *= mbas[*(ptr++)].par; 
      }
      if(par == parity) break;

      // new particle configuration

      for(loop = 0; (loop < num_part && ((matr[loop]+1) >= matr[loop+1])); loop++);
      matr[loop]++;
      for(loop--; loop >= 0; loop--)  matr[loop] = loop;

   } while(matr[num_part] == num_orb);  // end of particle loop

   if(matr[num_part] != num_orb) {
     printf("\n\nRank%d. Error in function pn_max_m_value():",Rank);
      printf("\nNo %s %d particle configuration with parity %c is possible\n",
	     (type == 0) ? "PROTON": "NEUTRON",num_part,
                          ((parity == +1) ? '+' : '-'));
       MPI_Abort(MPI_COMM_WORLD,Rank);
   }
   free(matr); // release temporary memory

   return m;  // return max m-value

} // End: function pn_max_m_value()

      /*
      ** The function 
      **          pn_modify_m_values()
      ** takes maximum and minimum m-values for |SD(Z)M P> and |SD(N)M P> for 
      ** a given proton parity  and modify the limits according to the condition 
      **                 tot_MJ = MZ + MN 
      ** The modified limits for protons are stored in mLim->max_MZ[type]
      ** and mLim->min_MZ[type].
      ** Note: All m_values are twice physical values.
      ** if (max_mz < min_mz) or (max_n < min_n) the function returns
      ** FALSE, otherwise TRUE.
      */

static int pn_modify_m_values(M_LIMITS *mLim, GR_BAS *grBas, int parity)
{
   int    max_mz, min_mz, tot_MJ, tot_P,
          max_mn, min_mn, mz, mn, step, num;

   // Parity identification: proton-parity plus (minus)--> parity = 0,(1)

   max_mz = mLim->max_MZ[parity];  // local variable
   min_mz = mLim->min_MZ[parity]; 
   tot_MJ = grBas->MJ;
   tot_P  = grBas->P;

   // initial maximum and minimum neutron m-values

   max_mn = pn_max_m_value(NEUTRON,
     ((parity == 0)? tot_P : -tot_P),grBas->spBas[1]);
   min_mn = -max_mn;

   // find new max_mz and min_mn
 
   if((step = max_mz + min_mn - tot_MJ) > 0)  {
     max_mz -= step; // reduce max_MZ
   }
   else if(step < 0) {
     min_mn -= step; // increase min_MN
   }
   // number of steps from max_mz to min_mz

   for(num = 1, mz = max_mz, mn = min_mn;
      (mz >min_mz) && (mn < max_mn); mz -= 2, mn += 2, num++); 

   min_mz = max_mz - 2 * (num - 1);
   mLim->max_MZ[parity] = max_mz; // store and return the result
   mLim->min_MZ[parity] = min_mz;

   if((max_mz < min_mz) || (max_mn < min_mn)) return FALSE;

   return TRUE; 

} // End: function pn_modify_m_values()

        /*
	** The function 
	**    calcGroupData()
	** calculates and store in group[]
	** all specific group data for current
	** group pointed to by grPtr.
	** if current group does not contain 
	** any basis states |SD(Z),SD(N)>
	** function returns FALSE, otherwise TRUE
       	*/

static int  calcGroupData(GR_BAS *grBas,GROUP *grPtr)
{
   char         *func = {" calcGroupData(): "};
   int          k, num,dim,pSD, nSD,*occZ, *occN;

   //******** Part 1 - same for all groups 

   grPtr->mask[0] = grBas->spBas[0]->mask;
   grPtr->mask[1] = grBas->spBas[1]->mask;

   // calculate and check number of |SD(Z)> and |SD(N)>

   grPtr->numSD[0] = pn_number_of_SD(PROTON, grBas, grPtr);
   grPtr->numSD[1] = pn_number_of_SD(NEUTRON, grBas, grPtr);

   grPtr->numAmplStep = MIN(NumProcs,grPtr->numSD[1]);
 
   if((grPtr->numSD[0] * grPtr->numSD[1])== 0) return FALSE;

   grPtr->numSD_ZN = grPtr->numSD[0] * grPtr->numSD[1];
   grPtr->SD[0] = MALLOC(grPtr->numSD[0] + 1, ULL, func,"grPtr.SD[0][]");
   pn_particle_SD_config(PROTON, grBas, grPtr);
   grPtr->SD[1] = MALLOC(grPtr->numSD[1] + 1, ULL, func,"grPtr.SD[1][]");
   pn_particle_SD_config(NEUTRON, grBas, grPtr);
     
   // calculate and store sumOscN for each |SD(Z)> and |SD(N)>
     
   occZ = MALLOC(grPtr->numSD[0],int,func,"occ[PROTON]");
   occN = MALLOC(grPtr->numSD[1],int,func,"occ[NEUTRON]");
   occ_orbits(PROTON,grBas,grPtr,occZ);
   occ_orbits(NEUTRON,grBas,grPtr,occN);

       /* 
       ** calculate number of  |SD(Z),SD(N)> 
       ** with sum_ZN_Osc_N < grBas->maxOscN and
       ** a grPtr->bitArray[] to specify the excluded
       ** basis states
       */

   dim = grPtr->numSD_ZN;

   grPtr->bitArray =
       MALLOC(BARR_ARRAYSIZE(dim),BARR_ELTYPE,func,"bitArray[]");
   BARR_CLEARARRAY(grPtr->bitArray,dim);
   
        /*
	** All bitArray[] elements are initialized to 0  
	** For all |SD(Z),SD(N)> with SumOscN > maxOscN, 
	** the corresponding  bitArray[] sre set to 1 
	*/

   grPtr->numSD_ZN_oscLim = grPtr->numSD_ZN;  
   num = 0;
   for(pSD = 0; pSD <grPtr->numSD[0]; pSD++)   {
     for(nSD = 0; nSD < grPtr->numSD[1];nSD++,num++) {
       if((occZ[pSD] + occN[nSD]) > grBas->maxOscN)  {
	 BARR_SET(grPtr->bitArray,num);
	 grPtr->numSD_ZN_oscLim--;
       } // end reduce
     } // end nSD
   } // end pSD

   if( grPtr->numSD_ZN_oscLim == 0) {
     free(grPtr->bitArray); 
     free(occN);
     free(occZ);
     free(grPtr->SD[1]);
     free(grPtr->SD[0]);
     return FALSE;   //  NO basis states in current group
   }

   // ******* Part 2 - different calculation for each Rank
     
   grPtr->numSD[2] = 0;
   for(k = Rank; k < grPtr->numSD[1]; k += NumProcs) {
     grPtr->numSD[2]++;
   }
   grPtr->RankNumSD_ZN = grPtr->numSD[0] * grPtr->numSD[2];
   grPtr->RankNumSD_ZN_oscLim = grPtr->RankNumSD_ZN;  
   num = 0;
   for(pSD = 0; pSD < grPtr->numSD[0]; pSD++)   {
     for(nSD = Rank; nSD < grPtr->numSD[1];nSD += NumProcs,num += NumProcs) {
       if((occZ[pSD] + occN[nSD]) > grBas->maxOscN)  {
	 grPtr->RankNumSD_ZN_oscLim--;
       } // end reduce
     } // end nSD
   } // end pSD

   return TRUE;

}  // End: function calcGroupData()

     /*
     ** The function 
     **           pn_number_of_SD()
     ** calculates the number of slater determinants for a 
     ** given set of identical particles (num_part) with 
     ** fixed total m-value (m_value) and parity (parity).
     ** Possible limitations on number of particles in 
     ** j-orbits are checked through grBas->spBas[type].
     ** mask. No time-reversal symmetry.
     */

static int pn_number_of_SD(int type, GR_BAS const *grBas, GROUP const *group)
{
   char   *func = {"pn_number_of_SD(): "};
   int    loop, number, m, par, m_value, parity, num_part,
          num_orb,num_mask,maxOscN, oscN,*matr,*ptr;
   ULL    sd, sd_mask;
   MBAS   *mbas;
   MASK   mask;
   
   m_value  = group->m[type];  // initialization
   parity   = group->par[type];
   num_part = grBas->spBas[type]->part;
   num_orb  = grBas->spBas[type]->numm_orb;
   mbas     = grBas->spBas[type]->mbas;
   mask     = grBas->spBas[type]->mask;

   if(type == PROTON) {
     maxOscN = grBas->maxOscN_Z;
   }
   else {
     maxOscN = grBas->maxOscN_N;
   }

   if(num_part == 0) return FALSE; // no particles

   matr     = MALLOC(num_part + 1, int, func,"matr[]"); // temporary memory
   num_mask = mask.num;
   for(loop = 0; loop < num_part; loop++) matr[loop] = loop; // lowest config
   matr[num_part] = num_orb;

   number = 0;   // initialization
   do  { // particle loop
      sd   = ULL_ZERO;
      m    = 0;
      par  = +1;
      oscN = 0;
      ptr = matr;

      for(loop = 0; loop < num_part; loop++) {
	m    += mbas[*ptr].m; 
	par  *= mbas[*ptr].par; 
	oscN += mbas[*ptr].N;
	sd   += ULL_ONE << *(ptr++); 
      }
      if(  (m == m_value) && (par == parity)
	 &&(oscN <= maxOscN))   {
        if(num_mask == 0) {
	  number++;
	}
	else {  // check for particle limitations
	  for(loop = 0; loop < num_mask; loop++) {
	    sd_mask = sd & mask.list[loop];
	    par = 0;
	    while(sd_mask) {
	      sd_mask &= sd_mask - 1;
	      par++;
	    }
	    if(  (par < mask.lim[2*loop])
	       ||(par > mask.lim[2*loop + 1])) break;
	  }
	  if(loop == num_mask) number++;
	}
      }
      // new particle configuration

      for(loop = 0;
          (loop < num_part) && ((matr[loop]+1) >= matr[loop+1]); loop++);
      matr[loop]++;
      for(loop--; loop >= 0; loop--)   matr[loop] = loop;
   } while(matr[num_part] == num_orb);  // end of particle loop

   free(matr);  // release temporary memory

   return number;

} // End: function pn_number_of_SD()

     /*
     ** The function 
     **         pn_particle_SD_config()
     ** calculates all slater determinants |SD()MP> for a  
     ** given number of identical particle (num_part) with 
     ** fixed total m-value (m_value) and parity (parity).
     */

static void pn_particle_SD_config(int type, GR_BAS const *grBas, GROUP const *group)
{
   char          *func = {"pn_particle_SD_config(): "};
   int    loop, number, m, par, m_value, parity, num_part,
          num_orb,num_mask,numSD,oscN,*matr,*ptr, maxOscN;
   ULL    *sd, sd_mask;
   MBAS   *mbas;
   MASK   mask;

   m_value  = group->m[type];  // initialization
   parity   = group->par[type];
   num_part = grBas->spBas[type]->part;
   num_orb  = grBas->spBas[type]->numm_orb;
   mbas     = grBas->spBas[type]->mbas;
   mask     = grBas->spBas[type]->mask;

   if(type == PROTON) {
     maxOscN = grBas->maxOscN_Z;
   }
   else {
     maxOscN = grBas->maxOscN_N;
   }
   numSD    = group->numSD[type];
   sd       = group->SD[type];

   if(!num_part) {
     if(numSD != 1) {
       printf("\n\nRank%d. Error in function pn_particle_SD_config():",Rank);
        printf("\nIn calc part type = %d  |SD> with no symmetry", type);
        printf("\nPart numb = ZERO, number of |SD> = %d\n", numSD) ;
         MPI_Abort(MPI_COMM_WORLD,Rank);
     }
     *sd = ULL_ONE; // |SD> = |vacuum>
      return;
   }

   matr     = MALLOC(num_part + 1, int, func, "matr[]"); // temporary memory 
   for(loop = 0; loop < num_part; loop++) matr[loop] = loop; // lowest config
   matr[num_part] = num_orb;

   num_mask = mask.num;  // mask check
   number   = 0;         //|SD> check

   do  {        // particle loop
     *sd = ULL_ZERO; // initialization of |SD>
     m = 0;
     par = +1;
     ptr = matr;
     oscN = 0;

     for(loop = 0;loop < num_part; loop++) {
       m    += mbas[*ptr].m; 
       par  *= mbas[*ptr].par;
       oscN += mbas[*ptr].N;
       *sd += ULL_ONE << *(ptr++); 
     }
      if(  (m == m_value) && (par == parity)
         &&(oscN <= maxOscN))  {
         if(num_mask == 0) {
	   sd++;
	   number++;
	 }
         else {
	  for(loop = 0; loop < num_mask; loop++) {
	    sd_mask = (*sd) & mask.list[loop];
	    par = 0;
	    while(sd_mask) {
	      sd_mask &= sd_mask - 1;
	      par++;
	    }
	    if(  (par < mask.lim[2*loop])
	       ||(par > mask.lim[2*loop + 1])) break;
            }
            if(loop == num_mask) {
               sd++;
               number++;
            } // end new |SD>
         }  // end mask check

      } // end check m-value and parity

      // new particle configuration

      for(loop = 0; 
         (loop < num_part) && ((matr[loop]+1) >= matr[loop+1]); loop++);
      matr[loop]++;
      for(loop--; loop >= 0; loop--)  matr[loop] = loop;
   
   } while(matr[num_part] == num_orb);  // end of particle loop

   if(numSD != number)  { // Check the number of identical SD
     printf("\n\nRank%d. Error in function pn_particle_SD_config():",Rank);
      printf("\nIn calculating particle |SD> with no symmetry");
      printf("\nNumber of type = %d |SD> = %d",type, numSD);
      printf(" - current calculation gives num = %d\n", number);
      MPI_Abort(MPI_COMM_WORLD,Rank);
   }
   free(matr); // release temporary memory

} // End: function pn_particle_SD_config()

    /*
    ** The function 
    **      pnGlobalData()
    ** calculates and store in grBas
    ** data common for all groups
    */

static void pnGlobalData(GR_BAS *grBas, GROUP *group) 
{
  int      grNum;
  ULL      startAmp;
  GROUP    *grPtr;

  grBas->totNumSD_Z             = 0; // initialization 
  grBas->totNumSD_N             = 0;
  grBas->RankTotNumSD_N         = 0;
  grBas->totNumSD_ZN            = 0; 
  grBas->totNumSD_ZN_oscLim     = 0; 
  grBas->RankTotNumSD_ZN        = 0;
  grBas->RankTotNumSD_ZN_oscLim = 0;
  grBas->maxSD[0]               = 0;
  grBas->maxSD[1]               = 0;

  grPtr    = group;
  startAmp = 0;
  for(grNum = 0; grNum < grBas->num_gr; grNum++, grPtr++) {
    grBas->totNumSD_Z             += grPtr->numSD[0];
    grBas->totNumSD_N             += grPtr->numSD[1];
    grBas->RankTotNumSD_N         += grPtr->numSD[2];
    grBas->totNumSD_ZN            += grPtr->numSD_ZN;
    grBas->totNumSD_ZN_oscLim     += grPtr->numSD_ZN_oscLim;
    grBas->RankTotNumSD_ZN        += grPtr->RankNumSD_ZN;
    grBas->RankTotNumSD_ZN_oscLim += grPtr->RankNumSD_ZN_oscLim;
    grBas->maxSD[0]                = MAX(grBas->maxSD[0], grPtr->numSD[0]);
    grBas->maxSD[1]                = MAX(grBas->maxSD[1], grPtr->numSD[1]);

    grPtr->startAmp  = startAmp;
    startAmp        += (ULL)grPtr->numSD_ZN_oscLim;
  }

} // End: funtion  pnGlobalData()

    /* 
    ** The function
    **       occ_orbits()
    ** calculates and stores orbit numbers of all occupied
    ** orbits  in all |SD(type)> in present group 
    */

static void occ_orbits(int type, GR_BAS *grBas,GROUP *group, int *occ)
{
   int   num_part, loop,par, orb, sumN;
   ULL   *sd, pos, sd_value;
   MBAS  *mBasPtr;

   num_part = grBas->spBas[type]->part; // initialization
   sd       = group->SD[type];
   mBasPtr  = grBas->spBas[type]->mbas;

   for(loop = 0; loop < group->numSD[type]; loop++, sd++) {
     sd_value = *sd;
     sumN     = 0;
     for(par = 0, orb = 0, pos = ULL_ONE; par < num_part;
                                 par++, orb++, pos <<= 1)  {
       for(;!(sd_value & pos); orb++, pos <<= 1);
       sumN    += mBasPtr[orb].N; 
     }
     *(occ++) = sumN;
   } // end loop through all |SD> in sd[]

} // End: function occ_orbits()


