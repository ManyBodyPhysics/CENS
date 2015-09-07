//      *******  The module lanc-sp-basis.c  ******

#include "shell-model.h"

   /*
   ** The entrance functions 
   **  void single_particle_basis(SHELL *model, SP_BAS *sp_bas)
   ** calculates and stores in basis.sp_basis 
   ** identical particle m-scheme single-particle basis:
   */

         //  **** local function declarations ****

static void single_particle_orbits(int part, int numj_orb, JBAS *jbas, SP_BAS *sp_bas);
     /*
     ** calculates and stores the complete 
     ** set of m-scheme single-particle states 
     */

static void single_m_states(int num_j, JBAS *jbas, MBAS *mbas);
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

static void part_orb_limitations(SP_BAS *sp_bas);
    /*
    ** analyzes all spherical input j-orbits for possible restrictions on
    ** particle occupations. For orbits with restrictions a mask to be used
    ** for the m-orbits is generated and stored in PART_LIM part_lim[].
    */

static void print_single_particle_orbits(char *file_name, SP_BAS *sp_bas);
   /*
   ** writes m-basis single-particle orbits for m > 0  to output file.
   */

            //  **** End: function declarations ****

            //  **** The function definitions  ****

    /*
    ** The entrance function 
    **         single_particle_basis()
    ** calculates and stores in basis.sp_basis 
    ** identical particle m-scheme single-particle basis:
    */

void single_particle_basis(SHELL *model, SP_BAS *sp_bas)
{
  single_particle_orbits(model->part, model->numj_orb, model->jbas, sp_bas);
  
  part_orb_limitations(sp_bas);              // possible j-orbit mask structure

  print_single_particle_orbits(model->title, sp_bas);

} // End: function single_particle_basis()

    /*
    ** The function
    **     single_particle_orbits()
    ** calculates and stores the complete 
    ** set of m-scheme single-particle states 
    */

static void single_particle_orbits(int part, int numj_orb, JBAS *jbas, SP_BAS *sp_bas)
{
  char
             *func = {"pn_single_particle_orbits(): "};
  int
             loop;

  sp_bas->numj_orb = numj_orb;         // number of j-orbits

        /*
        ** From spherical orbits create a m-basis for |SD>. Store
        ** the result in sp_bas.numm_orm and sp_bas.mbas[].
        */

  sp_bas->numm_orb = 0;                     // number of m-orbits
  for(loop = 0; loop < sp_bas->numj_orb; loop++) {
    sp_bas->numm_orb += jbas[loop].j + 1;                               // 2*j+1
  }
             // Number of m-orbits must not exceed size of int

  if(sp_bas->numm_orb > MAX_M_ORBITS) {
    printf("\nError in function single_particle_orbits():");
    printf("\nNumber of single-partice m-orbits is too large.");
    printf("\nnumber of m-orbits = %d  size of int = %d bits\n",
                                              sp_bas->numm_orb, MAX_M_ORBITS);
    exit(1);
  }
  if(part > sp_bas->numm_orb) {        /* Check number of particles against number of m-orbits */
    printf("\nError from function single_particle_orbits():");
    printf("\nNumber of particles = %d", part);
    printf("  number of m_orbits = %d\n", sp_bas->numm_orb);
    exit(1);
  }
  sp_bas->mask.num = 0;                                     // number of mask orbits
  for(loop = 0; loop < sp_bas->numj_orb; loop++) {
    if(  (jbas[loop].min_part > 0) 
       ||(jbas[loop].max_part < jbas[loop].j + 1)) sp_bas->mask.num++;
  }

  sp_bas->jbas = MALLOC(sp_bas->numj_orb, JBAS, func,"jbas[]");
  sp_bas->mbas = MALLOC(sp_bas->numm_orb + 1, MBAS, func, " mbas[]"); 

  if(sp_bas->mask.num > 0) {
    sp_bas->mask.lim  = MALLOC(2*sp_bas->mask.num,int, func,"mask_lim[]");
    sp_bas->mask.list = MALLOC(sp_bas->mask.num,ULL, func,"mask_list[]");
  }

  for(loop = 0; loop < sp_bas->numj_orb; loop++) {      // transfer all jbas[] data
    sp_bas->jbas[loop].osc      = jbas[loop].osc;
    sp_bas->jbas[loop].l        = jbas[loop].l;
    sp_bas->jbas[loop].j        = jbas[loop].j;
    sp_bas->jbas[loop].min_part = jbas[loop].min_part;
    sp_bas->jbas[loop].max_part = jbas[loop].max_part;
    sp_bas->jbas[loop].e        = jbas[loop].e;
  }
      // calculate sp_bas->mbas[] m-scheme orbits

  single_m_states(sp_bas->numj_orb, sp_bas->jbas, sp_bas->mbas); 

}  // End: function single_particle_orbits()

    /*
    ** The function                                     
    **             single_m_states()                    
    ** takes the set of spherical orbits and creates an m-basis 
    **for the slater determinants. The result are stored in 
    ** sp_bas->mbas[].
    */

static void single_m_states(int num_j, JBAS *jbas, MBAS *mbas)
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
      ptr_m->phase = (MOD(ptr_j->j+1,4)/2) ? -1: +1; // "Whitehead" phase factor"
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

  qsort(mbas,(size_t) num, (size_t)sizeof(MBAS),
            (int(*)(const void *,const void *))m_comp);

  ptr = ptr_m - 1;
  do {                                                       // all m < 0 
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

} // End: function single_m_states()

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
  else if(one->m == two->m)   {                             // one->m == two->m 
    if(one->j > two->j)         return -1;
    else if(one->j < two->j)    return +1;
    else {                                                // one->j == two->j
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
  return +1;
} // End: function m_comp()

   /*
   ** The function 
   **               part_orb_limitations()
   ** analyzes all spherical input j-orbits for possible restrictions on
   ** particle occupations. For orbits with restrictions a mask to be used
   ** for the m-orbits is generated and stored in sp_bas->mask_list[]
   ** and upper and lower limits sp_bas->mask_lim[]
   */

static void part_orb_limitations(SP_BAS *sp_bas)
{
  int 
             k, m, osc, l, j, num_mask;
  ULL
             mask;
  JBAS
             *ptr_jbas;
  MBAS
             *ptr_mbas;

  ptr_jbas = sp_bas->jbas;                      // initialization

  num_mask = 0;                                     // number of mask orbits
  for(k= 0; k < sp_bas->numj_orb; k++, ptr_jbas++) {
    mask = ULL_ZERO;
    if((ptr_jbas->min_part > 0) || (ptr_jbas->max_part < (ptr_jbas->j + 1)))  {
      sp_bas->mask.lim[2*num_mask]     = ptr_jbas->min_part;     // particle limitations
      sp_bas->mask.lim[2*num_mask + 1] = ptr_jbas->max_part;

      osc  = ptr_jbas->osc;
      l    = ptr_jbas->l;
      j    = ptr_jbas->j;
      for(m = 0, ptr_mbas = sp_bas->mbas; m < sp_bas->numm_orb; m++, ptr_mbas++) {
	if((ptr_mbas->osc == osc) && (ptr_mbas->l == l) && (ptr_mbas->j == j)) {
	  mask |= (ULL_ONE << m);
	}
      } // loop through all m_orbits
      sp_bas->mask.list[num_mask] = mask;      // save MASK value
      num_mask++;
    } // if-test for j_orbit reduction

  } // loop through all j_orbits

  if(sp_bas->mask.num != num_mask)  {
    printf("\n\nError in function part_orb_limitations():");
    printf("\n Wrong number of mask orbits found ");
    printf("\n sp_bas->num_mask = %d  -- current value = %d\n\n",
            	                       sp_bas->mask.num, num_mask);
    exit(1);
  }

} // End: function  part_orb_limitations()

   /*
   ** The function 
   **     print_single_particle_orbits()
   ** writes m-basis single-particle orbits for m > 0  to output file.
   */

static void print_single_particle_orbits(char *title, SP_BAS *sp_bas)
{
  char
             filename[ONE_LINE];
  int
             num_m, k, max;
  div_t
             division;
  MBAS
             *mbas, *ptr_mbas;
  FILE
             *out_data; 

  num_m = sp_bas->numm_orb;
  mbas  = sp_bas->mbas;

  strcpy(filename, title);
  strcat(filename,RESULT_OUTPUT);
  if( (out_data = fopen(filename,"a")) == NULL) {           // open out_data file 
    printf("\n\nError in function print_single_particle_orbits():");
    printf("\nWrong file = %s for the output data\n", filename);
    exit(1);
  }
  fprintf(out_data,                       //  Out data text heading
	  "\n\nsingle-particle orbits for m >= 1/2 (symmetric around zero)\n");
  fprintf(out_data,"\nm_orb j_orb osc l par 2*j 2*m energy *");
  fprintf(out_data,"  m_orb j_orb osc l par 2*j 2*m energy\n");

  division = div(num_m / 2, 2);
  max      = division.quot + division.rem;
  ptr_mbas = mbas;
  for(k = 0; k < max; k++)  {
    fprintf(out_data,"\n %2d    %2d   %2d %2d  %c  %2d  %3d  %5.2f *",
          k, ptr_mbas[k].orb,   ptr_mbas[k].osc,
             ptr_mbas[k].l,   ((ptr_mbas[k].par == +1) ? '+': '-'),
             ptr_mbas[k].j,     ptr_mbas[k].m,
             ptr_mbas[k].e);
    if((k == max - 1) && (division.rem)) continue;
    fprintf(out_data,"   %2d    %2d   %2d %2d  %c  %2d  %3d  %5.2f",
           k+max, ptr_mbas[k+max].orb,  ptr_mbas[k+max].osc,
              ptr_mbas[k+max].l, ((ptr_mbas[k+max].par == +1) ? '+': '-'),
              ptr_mbas[k+max].j,   ptr_mbas[k+max].m,
              ptr_mbas[k+max].e);
  }
  if(fclose(out_data)) {                                  // close output data file
    printf("\n\nError in function print_single_particle_orbits():");
    printf("\nIn closing the output data file\n");
    exit(1);
  }
} // End: print_single_particle_orbits()
