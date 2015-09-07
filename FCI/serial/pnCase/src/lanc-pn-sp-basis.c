/*******************  The module lanc-pn-sp-basis.c  ******************/

#include "shell.h"

   /*
   ** The entrance function 
   **   void pn_single_particle_basis(char *file_name, PN_SHELL *model)
   ** calculates proton/neutron single-particle basis:
   **    1,  proton single-particle m-orbits
   **    2.  particle limitations in proton orbits 
   **    3.  neutron single-particle m-orbits
   **    3.  particle limitations in neutron orbits
   */
             /**** local data definitions ****/


	     /**** local function declarations ****/

static void pn_single_particle_orbits(int type, PN_SHELL *model, GR_BASIS *gr_basis);
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
    **   3.  For orbits with the same (j,m)--value, positive parity comes first
    */

static PART_LIM pn_part_orb_limitations(int type, PN_SHELL *model, GR_BASIS *gr_basis);
    /*
    ** analyzes all spherical input j-orbits for possible restrictions on
    ** particle occupations. For orbits with restrictions a mask to be used
    ** for the m-orbits is generated and stored in PART_LIM part_lim[].
    */

static void print_single_particle_orbits(char const *out_file, int type, GR_BASIS const *gr_basis); 
    /*
    ** writes m-basis single-particle orbits for m > 0  to output file.
    */

static int pn_parity_check(PN_SHELL const *model, GR_BASIS *gr_basis);
    /*
    ** calculates and return the possible parities for the proton |SD(Z)M P>.
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

void pn_single_particle_basis(PN_SHELL *model, GR_BASIS *gr_basis)
{
   	 /* proton m-scheme single-particle basis */

  pn_single_particle_orbits(PROTON, model, gr_basis);
         
  gr_basis->part_lim[PROTON] = pn_part_orb_limitations(PROTON, model, gr_basis);

   	 /* neutron m-scheme single-particle basis */

  if(model->same_pn_orb == NO) {
    gr_basis->m_orb[1]    = gr_basis->m_orb[0];
    gr_basis->mbas[1]     = gr_basis->mbas[0];
    gr_basis->part_lim[1] = gr_basis->part_lim[0];
  }
  else {
    pn_single_particle_orbits(NEUTRON, model, gr_basis);
    gr_basis->part_lim[NEUTRON] = pn_part_orb_limitations(NEUTRON, model, gr_basis);
  }
           /* 
	   ** possible parities of proton |SD(Z)>               
           ** coding: -1: parZ  -1; both parZ = 0; +1: parZ = +1
           */

  gr_basis->parZ = pn_parity_check(model, gr_basis);               /* proton parities */

      /* print initial proton-neutron data */

  print_single_particle_orbits(model->out_data, PROTON, gr_basis);
  print_single_particle_orbits(model->out_data, NEUTRON, gr_basis);

} /* End_function pn_single_particle_basis() */

    /*
    ** The function
    **     pn_single_particle_orbits()
    ** calculate and store for particle type = 0 (proton),
    ** 1 (neutron) the complete single-particle basis 
    */

static void pn_single_particle_orbits(int type, PN_SHELL *model, GR_BASIS *gr_basis)
{
  char
             *func = {"pn_single_particle_orbits(): "};
  int
             part, num_j, num_m, loop;
  JBAS
             *jbas, *ptr_jbas;
  MBAS
             *mbas;

  part  = model->part[type];                          /* initialization */
  num_j = model->num_j[type];
  jbas  = model->jbas[type];

        /*
        ** From spherical orbits create a m-basis for |SD>. Store
        ** the result in model.num_m and model.mbas[].
        */

  ptr_jbas = jbas;           
  for(loop = 0, num_m = 0, ptr_jbas = jbas; loop < num_j; loop++) {
    num_m += (ptr_jbas++)->j + 1;        /* 2*j+1 */
  }
       /* Number of m-orbits must not exceed size of int */ 

  if(num_m > MAX_M_ORBITS) {
    printf("\nError in function pn_single_particle_orbits():");
    printf("\nNumber of single-partice m-orbits for %s is too large.",
                                           (type ? "neutron" : "proton"));
    printf("\nnumber of m-orbits = %d  size of int = %d bits\n",
                                              num_m, MAX_M_ORBITS);
    exit(1);
  }
  if(part > num_m) {        /* Check number of particles against number of m-orbits */
    printf("\nError from function pn_ single_particle_orbits():");
    printf("\nType of particle = %d number of particles = %d", type, part);
    printf("  number of m_orbits = %d\n",num_m);
    exit(1);
  }
 
  mbas = MALLOC((num_m + 1), MBAS, func, " mbas[]"); /* memory for m-scheme orbits */ 
  pn_single_m_states(num_j, jbas, mbas);                /* calculate m-scheme orbits */

  gr_basis->m_orb[type] = num_m;             /* return data for m-sheme orbits */
  gr_basis->mbas[type]  = mbas;
 
}  /* End: function pn_single_particle_orbits() */

    /*
    ** The function                                     
    **             pn_single_m_states()                    
    ** takes for particle type = 0 (proton), 1 (neutron) the set of spherical
    ** orbits and creates an m-basis for the slater determinants. The result
    ** are stored in model.mbas[type][].
    */

static void pn_single_m_states(int num_j, JBAS *jbas, MBAS *mbas)
{
  int
              k, m, num;
  JBAS
              *ptr_j;
  MBAS
              *ptr_m, *ptr;

  for(k = 0,num = 0,ptr_j = jbas,ptr_m = mbas;k < num_j; ptr_j++, k++) { /* all m > 0 */
    for(m = ptr_j->j; m > 0; ptr_m++, m -= 2) {
      ptr_m->orb   = k;
      ptr_m->osc   = ptr_j->osc;
      ptr_m->l     = ptr_j->l;
      ptr_m->par   = (ptr_j->l % 2) ? -1 : +1;
      ptr_m->j     = ptr_j->j;
      ptr_m->m     = m;
      ptr_m->e     = ptr_j->e;
      ptr_m->phase = (MOD(ptr_j->j+1,4)/2) ? -1: +1; /* "Whitehead" phase factor" */   
      num++;
    }
  }
       /* 
       ** sort all m > 0 m-orbits according to the rules:   
       **   1.  orbits are ordered after decreasing single-particle m-value.
       **       Orbits k and l with k > l has m_k < m_l
       **   2.  For orbits with same m-values, highest j comes first
       **   3.  For orbits with the same (j,m)--value, positive parity comes first
       */

  qsort(mbas,(UL) num,sizeof(MBAS),(int(*)(const void *,const void *))m_comp);

  ptr = ptr_m - 1;
  do {                                                       /* all m < 0 */
    ptr_m->orb       = ptr->orb;
    ptr_m->osc       = ptr->osc;
    ptr_m->l         = ptr->l;
    ptr_m->par       = ptr->par;
    ptr_m->j         = ptr->j;
    ptr_m->m         = -(ptr->m);
    ptr_m->e         = (ptr--)->e;
    (ptr_m++)->phase = +1;
  } while(--num);
  return;

} /* End: function pn_single_m_states()  */

    /*
    ** The function                         
    **        int m_comp()                  
    ** is a utility function for the library function qsort() in order to
    ** sort the single-particle m-orbits in mbas[] according to the rules:
    **   1.  orbits are ordered after decreasing single-particle m-value.
    **       Orbits k and l with k > l has m_k < m_l
    **   2.  For orbits with same m-values, highest j comes first
    **   3.  For orbits with the same (j,m)--value, positive parity comes first
    */

static int m_comp(const MBAS *one, const MBAS *two)
{
  if(one->m > two->m)            return -1;
  else if(one->m < two->m)       return +1;
  else if(one->m == two->m)   {                             /* one->m == two->m */ 
    if(one->j > two->j)         return -1;
    else if(one->j < two->j)    return +1;
    else {                                                /* one->j == two->j */
      if(one->par == +1)       return -1;
      else                     return +1;
    }
  }
  else {
    printf("\n\nError in function m_comp():");
    printf("\nTwo single-particle orbits have same quantum numbers!");
    printf("\norbit one: par = %d  2j = %d  2m = %d", one->par, one->j, one->m);
    printf("\norbit two: par = %d  2j = %d  2m = %d", two->par, two->j, two->m);
    exit(1);
  }
} /* End: function m_comp() */

   /*
   ** The function 
   **               pn_part_orb_limitations()
   ** analyzes all spherical input j-orbits for possible restrictions on
   ** particle occupations. For orbits with restrictions a mask to be used
   ** for the m-orbits is generated and stored in PART_LIM part_lim[].
   */

static PART_LIM pn_part_orb_limitations(int type, PN_SHELL *model, GR_BASIS *gr_basis)
{
  char
             *func = {"pn_part_orb_limitations(): "};
  int 
             num_j, num_m, k, m, osc, l, j, num_mask, max;
  UL
             mask;
  JBAS
             *jbas, *ptr_jbas;
  MBAS
             *mbas, *ptr_mbas;
  PART_LIM
             part_lim;
  MASK
             *table_ptr;

  num_j    = model->num_j[type];                        /* initialization */
  jbas     = model->jbas[type];
  num_m    = gr_basis->m_orb[type];
  mbas     = gr_basis->mbas[type];

  num_mask = 0;                                     /* number of mask orbits */
  for(k= 0, ptr_jbas = jbas; k < num_j; k++, ptr_jbas++) {
    max = jbas[k].j + 1;
    if(   (ptr_jbas->min_part < 0) || (ptr_jbas->min_part > max)  /* error test */
       || (ptr_jbas->max_part < 0) || (ptr_jbas->max_part > max)
       || (ptr_jbas->max_part < ptr_jbas->min_part))      {

      printf("\n\nError in function pn_part_orb_limitations():");
      printf("\nlimitations on %s particle number in j-orbits are wrong!!",
                                                   (type ? "neutron" : "proton"));
      printf("\norb_num = %d  j_value = %d",k, ptr_jbas->j);
      printf(" min_particle_number = %d  max_particle_numbler = %d",
                                     ptr_jbas->min_part, ptr_jbas->max_part);  
      exit(1);
    }
    if((ptr_jbas->min_part > 0) || (ptr_jbas->max_part < max)) num_mask++;
  }
  part_lim.num = num_mask;
  if(num_mask >0 ) {
    part_lim.table = MALLOC(num_mask, MASK, func,"part_lim.table[];");
    table_ptr = part_lim.table;
    for(k = 0, ptr_jbas = jbas; k < num_j; k++, ptr_jbas++) {
      max = ptr_jbas->j + 1;
      if((ptr_jbas->min_part > 0) || (ptr_jbas->max_part < max)) {
	table_ptr->orb = k;
	table_ptr->min = ptr_jbas->min_part;
	table_ptr->max = ptr_jbas->max_part;
	osc  = ptr_jbas->osc;
	l    = ptr_jbas->l;
	j    = ptr_jbas->j;
	mask = 0;
	for(m = 0, ptr_mbas = mbas; m < num_m; m++, ptr_mbas++) {
	  if((ptr_mbas->osc == osc) && (ptr_mbas->l == l) && (ptr_mbas->j == j)) {
	    mask |= (UL_ONE << m);
	  }
	}
	(table_ptr++)->mask = mask;
      }
    }
  }
  else   part_lim.table = NULL_PTR;                 /* no mask's */

  return part_lim;

} /* End: function  pn_part_orb_limitations() */

   /*
   ** The function 
   **     print_single_particle_orbits()
   ** writes m-basis single-particle orbits for m > 0  to output file.
   */

static void print_single_particle_orbits(char const *out_file, int const type, GR_BASIS const *gr_basis)
{
  int
             num_m, k, max;
  div_t
             division;
  MBAS
             *mbas, *ptr_mbas;
  FILE
             *out_data; 

  num_m = gr_basis->m_orb[type];                                /* initialization */
  mbas  = gr_basis->mbas[type];

  if( (out_data = fopen(out_file,"a")) == NULL) {           /* open out_data file */
    printf("\n\nError in function print_single_particle_orbits():");
    printf("\nWrong file = %s for the output data\n", out_file);
    exit(1);
  }
  fprintf(out_data,                       /*  Out data text heading */
           "\n\n  %s single-particle orbits for m >= 1/2 (symmetric around zero)\n",
	                                         (type == 0) ? "Proton" : "Neutron"); 
  fprintf(out_data,"\nm_orb j_orb osc l par 2*j 2*m energy *");
  fprintf(out_data,"  m_orb j_orb osc l par 2*j 2*m energy\n");

  division = div(num_m / 2, 2);
  max      = division.quot + division.rem;
  ptr_mbas = mbas;
  for(k = 0; k < max; k++)  {
    fprintf(out_data,
             "\n %2d    %2d   %2d %2d  %c  %2d  %3d  %5.2f *",
          k, ptr_mbas[k].orb,   ptr_mbas[k].osc,
             ptr_mbas[k].l,   ((ptr_mbas[k].par == +1) ? '+': '-'),
             ptr_mbas[k].j,     ptr_mbas[k].m,
             ptr_mbas[k].e);
    if((k == max - 1) && (division.rem)) continue;
    fprintf(out_data,
          "   %2d    %2d   %2d %2d  %c  %2d  %3d  %5.2f",
           k+max, ptr_mbas[k+max].orb,  ptr_mbas[k+max].osc,
              ptr_mbas[k+max].l, ((ptr_mbas[k+max].par == +1) ? '+': '-'),
              ptr_mbas[k+max].j,   ptr_mbas[k+max].m,
              ptr_mbas[k+max].e);
  }
  if(fclose(out_data)) {                                  /* close output data file */
    printf("\n\nError in function print_single_particle_orbits():");
    printf("\nIn closing the output data file\n");
    exit(1);
  }
} /* End: print_single_particle_orbits() */

     /*
     ** The function 
     **        pn_parity_check()
     ** calculates and return the possible parities for the proton |SD(Z)M P>.
     ** The information is stored as 
     **           parity  = -1    negativ parity only
     **                   =  0    both parities       
     **                   = +1    positiv parity only
     ** Possible parities for both |SD(Z)P_Z> and |SD(N)P_N> is determined
     ** and the restriction tot_parity = P_Z * P_N is taken into account.
     */

static int pn_parity_check(PN_SHELL const *model, GR_BASIS *gr_basis)
{
  int
           k, l, tot_parity, par[2], part[2], num_m[2];
  MBAS
           *ptr, *mbas[2];

  tot_parity = model->P;                                  /* initialization */
  part[0]    = model->part[0];
  part[1]    = model->part[1];
  num_m[0]   = gr_basis->m_orb[0];
  num_m[1]   = gr_basis->m_orb[1];
  mbas[0]    = gr_basis->mbas[0];
  mbas[1]    = gr_basis->mbas[1];

  for(k = 0; k < 2; k++) {                 /* proton = [0], neutron = [1] */
    ptr    = mbas[k];                           /* initialization */
    par[k] = (ptr++)->par;
    for(l = 1; (l < num_m[k]) && (par[k] == ptr->par); l++, ptr++);
    if(l < num_m[k]) par[k] = 0;

            /* possible modification for odd particle number */

    if((par[k] == -1) && !(part[k] % 2)) par[k] = +1; 

            /* possible modification for closed particle shell */

    if(part[k] == num_m[k]) par[k] = +1;
  }
      /* possible reduction: tot_parity = par[0] * par[1] */

  if((par[0] == 0) && (par[1] != 0)) {
    par[0] = tot_parity * par[1];           /* reduction */
  }
  if((par[0] != 0) && (par[1] != 0) && (par[1] != tot_parity * par[0])) {  /* error */
    printf("\n\nError in function check_parity():");
    printf("\n proton-parity = %d  neutron-parity = %d", par[0], par[1]);
    printf("\n can not produce total parity = %d\n", tot_parity);
    exit(1);
  }       
  return par[0];

}  /* End: function pn_parity_check() */
