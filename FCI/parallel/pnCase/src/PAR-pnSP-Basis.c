
/*******************  The module pn-sp-basis.c  ******************/

#include "PAR-pnShellModel.h"

   /*
   ** The entrance function 
   **   void pn_single_particle_basis(GR_BAS *grBas)
   ** calculates proton/neutron single-particle basis:
   **    1,  proton single-particle m-orbits
   **    2.  particle limitations in proton orbits 
   **    3.  neutron single-particle m-orbits
   **    3.  particle limitations in neutron orbits
   */
             /**** local data definitions ****/


	     /**** local function declarations ****/

static void pn_single_particle_orbits(SP_BAS *spBas);
    /*
    ** calculate and store for particle type = 0 (proton),
    ** 1 (neutron) the complete single-particle orbits 
    */

static void pn_single_m_states(int num_j, JBAS *jbas, MBAS *mbas);
    /*
    ** takes for particle type = 0 (proton), 1 (neutron) the set of spherical
    ** orbits and creates an m-basis for the slater determinants. The result
    ** are stored in model.mbas[type][].
    */

static int m_comp(const MBAS *one, const MBAS *two);
    /*
    ** is a utility function for the library function qsort() in order to
    ** sort the single-particle m-orbits in mbas[] according to the rules:
    **   1.  orbits are ordered after decreasing single-particle m-value.
    **       Orbits k and l with k > l has m_k < m_l
    **   2.  For orbits with same m-values, highest j comes first
    **   3.  For orbits with the same (j,m)--value, 
    **       positive parity comes first
    */

static void pn_part_orb_limitations(SP_BAS *spBas);
    /*
    ** analyzes all spherical input j-orbits for possible restrictions on
    ** particle occupations. For orbits with restrictions a mask to be used
    ** for the m-orbits is generated and stored in PART_LIM part_lim[].
    */

static void print_single_particle_orbits(int type, GR_BAS const *grBas); 
    /*
    ** writes m-basis single-particle orbits for m > 0  to output file.
    */

static void pn_parity_check(GR_BAS *grBas);
    /*
    ** calculates and return the possible parities for proton |SD(Z)M P>.
    ** The information is stored as 
    **           parity  = -1    negativ parity only
    **                   =  0    both parities       
    **                   = +1    positiv parity only
    ** Possible parities for both |SD(Z)P_Z> and |SD(N)P_N> is determined
    ** and the restriction P(total_parity) = P_Z * P_N is taken into account.
    */

             /**** End: function declarations ****/

  /**** The function definitions  ****/ 

    /*
    ** The entrance function 
    **         pn_single_particle_basis()
    ** calculates proton/neutron single-particle basis: 
    **    1,  proton single-particle m-orbits
    **    2.  particle limitations in proton orbits 
    **    3.  neutron single-particle m-orbits
    **    3.  particle limitations in neutron orbits
    */

void pn_single_particle_basis(GR_BAS *grBas)
{

  pn_single_particle_orbits(grBas->spBas[PROTON]);

  pn_part_orb_limitations(grBas->spBas[PROTON]);

  // neutron m-scheme single-particle basis

  pn_single_particle_orbits(grBas->spBas[NEUTRON]);


  pn_part_orb_limitations(grBas->spBas[NEUTRON]);

       /* 
       ** possible parities of proton |SD(Z)>               
       ** coding grBas.parZ: -1: parZ  -1; 
       ** both parZ = 0; +1: parZ = +1
       */

  pn_parity_check(grBas);

  // print initial proton-neutron data

  if(Rank == MASTER) {
    print_single_particle_orbits(PROTON, grBas);
    print_single_particle_orbits(NEUTRON,grBas);
  } // end print MASTER

} // End_function pn_single_particle_basis() 

    /*
    ** The function
    **     pn_single_particle_orbits()
    ** calculates and stores the complete
    **  single-particle m-scheme basis 
    */

static void pn_single_particle_orbits(SP_BAS *spBas)
{
  char    *func = {"pn_single_particle_orbits(): "};
  int     numm_orb, loop;
  JBAS    *ptr_jbas;
  MBAS    *mbas;

        /*
        ** From spherical orbits create a m-basis for |SD>.
        ** Store the result in spBas->numm_orb and 
	** spBas->mbas[].
        */

  numm_orb = 0;
  ptr_jbas =  spBas->jbas;
  for(loop = 0; loop < spBas->numj_orb; loop++) {
    numm_orb += (ptr_jbas++)->j + 1;  // 2*j+1
  }

  // number of m-orbits must not exceed size of MAX_M_ORBITS

  if(numm_orb > (int)MAX_M_ORBITS) {
    printf("\nRank%d: Error in function pn_single_particle_orbits()",Rank);
    printf("\nNumber of single-partice m-orbits is too large");
    printf("\n number of m-orbits = %d  size of int = %d bits\n\n\n",
                                          numm_orb, (int)MAX_M_ORBITS);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }
  if(spBas->part > numm_orb) {  // check part num  against num m-orbits
    printf("\nRank%d: Error from function pn_ single_particle_orbits()",Rank);
    printf("  number of m_orbits = %d\n\n",numm_orb);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }
 
  mbas = MALLOC((numm_orb + 1), MBAS, func," mbas[]");

  pn_single_m_states(spBas->numj_orb, spBas->jbas, mbas);

  spBas->numm_orb = numm_orb; // return data for m-sheme orbits;
  spBas->mbas     = mbas;

} // End: function pn_single_particle_orbits()

    /*
    ** The function                                     
    **             pn_single_m_states()                    
    ** takes for particle type = 0 (proton), 1 (neutron) the set
    ** of spherical orbits and creates an m-basis for the slater 
    ** determinants. The result are stored in model.mbas[type][].
    */

static void pn_single_m_states(int num_j, JBAS *jbas, MBAS *mbas)
{
  int      k, m, num;
  JBAS     *ptr_j;
  MBAS     *ptr_m, *ptr;

  for(k = 0,num = 0,ptr_j = jbas,ptr_m = mbas;k < num_j;
                                          ptr_j++, k++) { //all m > 0
    for(m = ptr_j->j; m > 0; ptr_m++, m -= 2) {
      ptr_m->orb = k;
      ptr_m->N   = ptr_j->N;
      ptr_m->osc = ptr_j->osc;
      ptr_m->l   = ptr_j->l;
      ptr_m->par = (ptr_j->l % 2) ? -1 : +1;
      ptr_m->j   = ptr_j->j;
      ptr_m->m   = m;
      ptr_m->e   = ptr_j->e;

      // "Whitehead" phase factor"

      ptr_m->phase = (MOD(ptr_j->j+1,4)/2) ? -1: +1;

      num++;
    }
  }
       /* 
       ** sort all m > 0 m-orbits according to the rules:   
       **   1.  orbits are ordered after decreasing single-particle
       **       m-value. Orbits k and l with k > l has m_k < m_l
       **   2.  For orbits with same m-values, highest j comes first
       **   3.  For orbits with the same (j,m)--value, positive 
       **       parity comes first
       */

  qsort(mbas,(UL) num,sizeof(MBAS),(int(*)(const void *,const void *))m_comp);

  ptr = ptr_m - 1;
  do {                                                       /* all m < 0 */
    ptr_m->orb       = ptr->orb;
    ptr_m->N         = ptr->N;
    ptr_m->osc       = ptr->osc;
    ptr_m->l         = ptr->l;
    ptr_m->par       = ptr->par;
    ptr_m->j         = ptr->j;
    ptr_m->m         = -(ptr->m);
    ptr_m->e         = (ptr--)->e;
    (ptr_m++)->phase = +1;
  } while(--num);
  return;

} // End: function pn_single_m_states()

    /*
    ** The function                         
    **        int m_comp()                  
    ** is a utility function for the library function qsort()
    ** in order to sort the single-particle m-orbits in mbas[] 
    ** according to the rules:
    **   1.  orbits are ordered after decreasing single-particle
    **       m-value.Orbits k and l with k > l has m_k < m_l
    **   2.  For orbits with same m-values, highest j comes first
    **   3.  For orbits with the same (j,m)--value, positive parity
    **       comes first
    */

static int m_comp(const MBAS *one, const MBAS *two)
{
  if(one->m > two->m)            return -1;
  else if(one->m < two->m)       return +1;
  else if(one->m == two->m)   {        // one->m == two->m
    if(one->j > two->j)         return -1;
    else if(one->j < two->j)    return +1;
    else {         // one->j == two->j
      if(one->par == +1)       return -1;
      else                     return +1;
    }
  }
  else {
    printf("\n\nRank%d: Error in function m_comp():",Rank);
    printf("\nTwo single-particle orbits have same quantum numbers!");
    printf("\norbit one: par = %d  2j = %d  2m = %d",
                           one->par, one->j, one->m);
    printf("\norbit two: par = %d  2j = %d  2m = %d\n\n",
                              two->par, two->j, two->m);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }
  return 0;  // should never reach here

} // End: function m_comp()

   /*
   ** The function 
   **               part_orb_limitations()
   ** analyzes all spherical input j-orbits for possible
   ** restrictions on particle occupations. For orbits with 
   ** restrictions a mask to be used for the m-orbits is
   ** generated and stored in sp_bas->mask_list[] and
   ** upper and lower limits sp_bas->mask_lim[]
   */

static void pn_part_orb_limitations(SP_BAS *spBas)
{
  char    *func = {"pn_part_orb_limitations(): "};
  int     loop, k, m, osc, l, j, num_mask;
  ULL     mask;
  JBAS    *jbasPtr;
  MBAS    *mbasPtr;

  // mask orbits: check for any particle limitation in jBas

  spBas->mask.num = 0;
 jbasPtr          = spBas->jbas;
 for(loop = 0; loop < spBas->numj_orb; loop++,jbasPtr++) {
    if(  (jbasPtr->min_part > 0) 
       ||(jbasPtr->max_part < jbasPtr->j + 1)) spBas->mask.num++;
  }

  if(spBas->mask.num > 0) {
    spBas->mask.lim  = MALLOC(2*spBas->mask.num,int, func,"mask_lim[]");
    spBas->mask.list   = MALLOC(spBas->mask.num,ULL, func,"mask_list[]");
  }
        /*
        ** From spherical orbits create a m-basis for |SD>.
        ** Store the result in grBas->spbas[].num_m and 
	** grBas->spbas[].mbas[].
        */

 jbasPtr   = spBas->jbas;
  num_mask = 0; // number of mask orbits
  for(k= 0; k < spBas->numj_orb; k++, jbasPtr++) {
    mask = ULL_ZERO;
    if((jbasPtr->min_part > 0) || (jbasPtr->max_part < (jbasPtr->j + 1)))  {

      // particle limitation

      spBas->mask.lim[2*num_mask]     = jbasPtr->min_part;
      spBas->mask.lim[2*num_mask + 1] = jbasPtr->max_part;

      osc     = jbasPtr->osc;
      l       = jbasPtr->l;
      j       = jbasPtr->j;
      mbasPtr = spBas->mbas;

      for(m = 0; m < spBas->numm_orb; m++, mbasPtr++) {
	if((mbasPtr->osc == osc) && (mbasPtr->l == l) && (mbasPtr->j == j)) {
	  mask |= (ULL_ONE << m);
	}
      } // end loop through all m_orbits

      // save MASK value

      spBas->mask.list[num_mask] = mask;  
      num_mask++;
    } // if-test for j_orbit reduction
  
  } // loop through all j_orbits

  if(spBas->mask.num != num_mask) {
    printf("\n\nRank%d: Error in function part_orb_limitations()",Rank);
    printf("\n Wrong number of mask orbits found ");
    printf("\n sp_bas->num_mask = %d  -- current value = %d\n\n",
                       spBas->mask.num, num_mask);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }

} // End: function  pn_part_orb_limitations()

   /*
   ** The function 
   **     print_single_particle_orbits()
   ** writes m-basis single-particle orbits for m > 0  to output file.
   */

static void print_single_particle_orbits(int const type, GR_BAS const *grBas)
{
  char       filename[ONE_LINE];
  int        num_m, k, max;
  div_t      division;
  SP_BAS     *spBas;
  MBAS       *mbas, *ptr_mbas;
  FILE       *out_data; 

  spBas = grBas->spBas[type];

  num_m = spBas->numm_orb;   // initialization
  mbas  = spBas->mbas;
  sprintf(filename,"%s%s",grBas->title,RESULT_OUTPUT); 
  if( (out_data = fopen(filename,"a")) == NULL) {
    printf("\n\nRank%d: Error in function print_single_particle_orbits()",Rank);
    printf("\nWrong file = %s for the output data\n\n",filename);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }
  fprintf(out_data,   //  Out data text heading
      "\n\n  %s single-particle orbits for m >= 1/2 (symmetric around zero)\n",
	                                   (type == 0) ? "Proton" : "Neutron"); 
  fprintf(out_data,"\nm_orb j_orb N  n  l par 2*j 2*m energy *");
  fprintf(out_data,"  m_orb j_orb N  n  l par 2*j 2*m energy\n");

  division = div(num_m / 2, 2);
  max      = division.quot + division.rem;
  ptr_mbas = mbas;

  for(k = 0; k < max; k++)  {
    fprintf(out_data,
             "\n%2d    %2d   %2d %2d %2d  %c  %2d  %3d %6.2f *",
          k, ptr_mbas[k].orb, ptr_mbas[k].N, 
             ptr_mbas[k].osc,  ptr_mbas[k].l,
             ((ptr_mbas[k].par == +1) ? '+': '-'),ptr_mbas[k].j,
             ptr_mbas[k].m, ptr_mbas[k].e);
    if((k == max - 1) && (division.rem)) continue;
    fprintf(out_data,
          "  %2d    %2d   %2d %2d %2d  %c  %2d  %3d  %6.2f",
           k+max, ptr_mbas[k+max].orb,  ptr_mbas[k+max].N,
               ptr_mbas[k+max].osc,     ptr_mbas[k+max].l,
              ((ptr_mbas[k+max].par == +1) ? '+': '-'),
              ptr_mbas[k+max].j,   ptr_mbas[k+max].m,
              ptr_mbas[k+max].e);
  }
 
  if(fclose(out_data)) {  // close output data file
    printf("\n\nRank%d: Error in function print_single_particle_orbits():",
                                                                     Rank);
    printf("\nIn closing the output data file\n\n");
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }
} // End: print_single_particle_orbits()

     /*
     ** The function 
     **        pn_parity_check()
     ** calculates and return the possible parities for proton |SD(Z)M P>.
     ** The information is stored as 
     **           parity  = -1    negativ parity only
     **   |<                =  0    both parities       
     **                   = +1    positiv parity only
     ** Possible parities for both |SD(Z)P_Z> and |SD(N)P_N> is determined
     ** and the restriction tot_parity = P_Z * P_N is taken into account.
     */

static void pn_parity_check(GR_BAS *grBas)
{
  int     k, l, tot_parity, par[2], part[2], num_m[2];
  MBAS    *ptr, *mbas[2];

  tot_parity = grBas->P;    // initialization


  part[0]    = grBas->spBas[PROTON]->part;
  part[1]    = grBas->spBas[NEUTRON]->part;
  num_m[0]   = grBas->spBas[PROTON]->numm_orb;
  num_m[1]   = grBas->spBas[NEUTRON]->numm_orb;
  mbas[0]    = grBas->spBas[PROTON]->mbas;
  mbas[1]    = grBas->spBas[NEUTRON]->mbas;

  for(k = 0; k < 2; k++) {   // proton = [0], neutron = [1]
    ptr    = mbas[k];      // initialization
    par[k] = (ptr++)->par;
    for(l = 1; (l < num_m[k]) && (par[k] == ptr->par); l++, ptr++);
    if(l < num_m[k]) par[k] = 0;

    // possible modification for odd particle number

    if((par[k] == -1) && !(part[k] % 2)) par[k] = +1; 

    // possible modification for closed particle shell

    if(part[k] == num_m[k]) par[k] = +1;
  }
  // possible reduction: tot_parity = par[0] * par[1]

  if((par[0] == 0) && (par[1] != 0)) {
    par[0] = tot_parity * par[1];  // reduction
  }
  if((par[0] != 0) && (par[1] != 0) && (par[1] != tot_parity * par[0])) {
    printf("\n\nRank%dError in function check_parity():",Rank);
    printf("\n proton-parity = %d  neutron-parity = %d", par[0], par[1]);
    printf("\n can not produce total parity = %d\n\n", tot_parity);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }       
  grBas->parZ = par[0];

}  // End: function pn_parity_check()
