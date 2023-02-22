  /*
   * The program module                                 
   *           transition-lib.c                               
   * contains the following functions:                  
   *                                                    
   *      // Memory allocation functions //       
   *                                                    
   *    1. a. void   *vector( )                            
   *       b. void free_vector( )                          
   *       c. void   *vector_realloc( )                    
   *       d. void   **matrix( )                           
   *       e. void free_matrix( )                          
   *
   *      // Functions to read input data //
   *
   *    2. a. void read_orbit()
   *       b. int j_comp()                                
   *       c. int read_line( )                           
   *       d. void read_numbers_into_vector()
   *       e. char *removeBlanks() 
   *
   *      // Calculate m-states //
   *
   *    3. a. void single_m_states()
   *       b. part_orb_limitations()
   *       c. int mbas_parity()
   *
   *     // Create the single-particle basis //  
   *
   *     void read_transitions(void)
   *    4. a. int trans_comp()
   *       b. int plus_parity()
   *       c. int minus_parity()
   *       d. int plus_m_value()
   *       e. int minus_m_value()
   *       f. int config_comp()
   *    5. a. int binomial( )    
   *       b. void run_time( )  
   *       c. double clebsch_gordan()
   *       d. double fac_ratio()                  
   */

#include "shell.h"
#include "sm.h"

  /*
   * The function                         
   *         void  *vector()              
   * reserves dynamic memory in  heap for a vector using the 
   * function malloc(). No initialization of the elements.   
   * Input data:                          
   *   int num_bytes - number of bytes for each elements       
   *   int num_elem  - number of elements 
   * Returns a void pointer to the memory location.
   */

void  *vector(int num_bytes, int num_elem)
  {
    void  *pointer;

    if((pointer = (void *) malloc((size_t)(num_bytes * num_elem))) == NULL)  {
      printf("\nNo enough memory to store a vector with");
      printf(" \n%d number of %d bytes \n",num_elem,num_bytes);
      exit(1);
    }
    return pointer;
  } /* End: function  void *vector()  */

    /*
     * The function                        
     *      void free_vector()             
     * releases the memory reserved by the function vector() 
     * for the vector pointed to by void *vec             
     */

void free_vector(void *vec)
  {
    free(vec);
  } /* End: function free_vector()  */

  /*
   * The function                                       
   *               void  *vector_realloc()              
   * increases the dynamic memory in  heap for a vector 
   * using malloc(). No initialization of the elements. 
   * Input data:                                        
   *    void *old_ptr - pointer to the old vector       
   *          num     - total number of elements        
   */

void  *vector_realloc(void *old_ptr,int num)
{
  void  *new_ptr;

  new_ptr = (void  *) realloc( old_ptr,num);
  if( !new_ptr )  {
      printf("\nProblem with memory reallocation");
      printf(" \n%ul number of bytes \n",num);
      exit(1);
  }
  return new_ptr;
} /* End: function vector_realloc() */

  /*
   * The function                             
   *      void  **matrix()                    
   * reserves dynamic memory in far heap for a two-dimensional matrix 
   * using the function MALLOC(). No initialization of the elements. 
   * Input data:                      
   *  int row      - number of  rows          
   *  int col      - number of columns        
   *  int num_bytes- number of bytes for each 
   *                 element                  
   * Returns a void   **pointer to the location.                                
   */

void  **matrix(int row, int col, int num_bytes)
  {
    int          loop;
    void    **pointer;

    pointer = (void  **) malloc((ULL)( row * sizeof(void  *)));
    if(!pointer )  {
      printf("\nNo enough memory to store a matrix with");
      printf(" \nrow = %d  col = %d num_bytes = %d\n",row, col, num_bytes);
      exit(1);
    }
    for( loop = 0; loop < row; loop++)  {
      pointer[loop] = (void  *) malloc((ULL)(col * num_bytes));
      if( !pointer[loop] )  {
        printf("\nNo enough memory to store a matrix with");
        printf(" \ncol = %d\n",col);
        exit(1);
      }
    }
    return pointer;
  } /* End: function  void  **matrix()  */

    /*
     * The function                         
     *      void free_matrix()              
     * releases the memory reserved by the function matrix() for the 
     * two-dimensional matrix with number of rows given by the input data.
     * Input data:                          
     *  void far **pointer - pointer to the matrix
     *  int row          - number of rows
     */

void free_matrix(void **pointer, int row)
  {
    int loop;

    for( loop = 0; loop < row; loop++)  {
      free(pointer[loop]);
    }
    free(pointer);
  }  /* End:  function free_matrix() */

    /*
     * The function                       
     *         read_orbit()               
     * reads the data from file and store them in the structure Model.jbas[] 
     * after decreasing j-values.         
     */

void read_orbit(int num_j, struct jbas *jbas)
{
  char       ch;
  int        loop;
  double     e;

  for(loop = 0; loop < num_j; loop++) {
    fscanf(Files.in_data,"%d%d%d%d%d%lf",
        &(jbas[loop].osc),&(jbas[loop].l),&(jbas[loop].j),
	&(jbas[loop].min_part), &(jbas[loop].max_part), &e);
    do   {
       ch = fgetc(Files.in_data);
    } while((ch != 'a') && (ch != 'n') && (ch != '\n'));
    if(ch == '\n') {
      printf("\n\nError in in-data file for the single-particle orbits.\n");
      exit(1);
    }
    jbas[loop].e      = e;
    jbas[loop].status = ch;
  }     /* end of loop  */

    /* Re-order the j_obits in Model.j_bas[] after decreasing j-values.*/

  qsort(jbas,(ULL) num_j,sizeof(struct jbas),(int(*)(const void *,const void *)) j_comp);

} /* End: function read_orbit()  */

    /*
     * The function                         
     *        int j_comp()                  
     * is a utility function for the library function qsort() in order to
     * sort the single-particle orbits in jbas[] after decreasing j-values.
     */

int j_comp(const struct jbas *one, const struct jbas *two)
{
  if(one->j > two->j)       return -1;
  else  if(one->j < two->j) return +1;
  else if(one->j == two->j)  {
       if(PARITY(one->l)) return +1;
       else               return -1;
  }
  else                       return  0;
} /* End: function j_comp() */

    /*
     * The function                            
     *           read_line()                   
     * reads a line of characters from file and splits it in two - the first is a   
     * text string and the second is a data string. They are separated by : .
     *  Input data:                            
     *    Files.in_data     - pointer to file      
     *    char *text_string - contain the text by return   
     *    char *data_string - contain the data by return, not converted.   
     */

int read_line(char *text_string, char *data_string)
{
  char total_string[100], *ptr;

  if(fgets(total_string,100, Files.in_data) == NULL)  /* Read a line from file */
    return EOF;   /*Return EOF if end-of-file is reached */

  ptr =total_string;   /* initialization */

  while((*ptr) == ' ') ptr++;         /*Remove leading blanks */

  /* Extract the text string */ 

  while(((*ptr) != ':') && ((*ptr) != '\n') && ((*ptr) != '\r') && ((*ptr)!= 0x0))
     (*text_string++) = *ptr++;
   *text_string = 0x0;

  /* The data string */

  if(*ptr != ':') return TRUE;
  ptr++;                         /* step over :  */

  while((*ptr) == ' ') ptr++;  /* Remove leading blanks */

  /* Extract the data string */ 

  while(((*ptr) != '\n') && ((*ptr) != 0x0)) *data_string++ = *ptr++;
  *data_string = 0x0;

  return TRUE;
}  /* End: function read_line()  */

  /*
    * The function
    *         read_numbers_into_vector()
    * converts a character string containing 0 and 1 into numbers
    * and stored them in an integer vector list[].
    */

void read_numbers_into_vector(char *string, int number, int *list)
{
   char    *ptr;
   int     loop;

   ptr = string;      /* initialization */
   loop = 0;

   while((*ptr) != '\0')  {
      switch( *(ptr++))  {
         case  '0': list[loop++] = 0;
                    break;
         case  '1': list[loop++] = 1;
                    break;
         default  : break;
      } /* end switch()  */
   };  /* end while - loop */
   if(loop != number)  {
      printf("\n Error: Wrong input data in Start vector mixture!!");
      printf("\nNumber of start vectors mixture = %d   found = %d\n",number,loop);
      exit(1);
   }
   return;
} /* End: function read_number_into_vectors() */


   /*
    ** The function
    **           removeBlanks()
    ** removes all blanks in inString[] and return
    ** the result in outString[]
    */

void removeBlanks(char *string)
{
  char    *inPtr, *outPtr;

  inPtr  = string;
  outPtr = string;
  for( ; *inPtr != '\0'; ) {
    if(*inPtr != ' ') 
      *(outPtr++) = *(inPtr++);
    else inPtr++;
  }
  *outPtr = '\0';
} // End: function removeBlanks()  

    /*
     * The function                                     
     *             single_m_states()                    
     * takes the set of spherical orbits and creates an m-basis for the 
     * slater determinants. The result are store in orbit.mbas[].
     * A list of the m-states is written to the output file.
     */

void single_m_states(int num_j, struct jbas *jbas, struct m_orbit *orbit)
{
  int            k, m, count, max;
  struct mbas    *ptr_mbas;

  max = orbit->num - 1;
  ptr_mbas = orbit->mbas;

  for(count = 0, m = jbas[0].j; m > 0; m -= 2)  {
    for(k = 0;(jbas[k].j >= m)&&(k < num_j); count++,k++) {

              /* contribution to + m */

	ptr_mbas[count].orb    = k;
	ptr_mbas[count].osc    = jbas[k].osc;
	ptr_mbas[count].l      = jbas[k].l;
	ptr_mbas[count].par    = (jbas[k].l % 2) ? -1: +1;
	ptr_mbas[count].j      = jbas[k].j;
	ptr_mbas[count].m      = m;
	ptr_mbas[count].e      = jbas[k].e;
	ptr_mbas[count].status = (jbas[k].status == 'a') ? 1 : 0;
	ptr_mbas[count].phase  = (MOD(jbas[k].j+1,4)/2) ? -1:+1;

              /* contribution to - m  */

	ptr_mbas[max - count].orb    = k;
	ptr_mbas[max - count].osc    = ptr_mbas[count].osc;
	ptr_mbas[max - count].l      = ptr_mbas[count].l;
	ptr_mbas[max - count].par    = ptr_mbas[count].par;
	ptr_mbas[max - count].j      = ptr_mbas[count].j;
	ptr_mbas[max - count].m      = - m;
	ptr_mbas[max - count].e      = ptr_mbas[count].e;
	ptr_mbas[max - count].status = ptr_mbas[count].status;
	ptr_mbas[max - count].phase  = +1;

    } /* end of loop k  */
  } /* end of loop m  */

      /* print a list of m-states to the output file */

  fprintf(Files.out_data,"\nm_orb j_orb osc l par 2*j 2*m energy *");
  fprintf(Files.out_data,"  m_orb j_orb osc l par 2*j 2*m energy\n");

  max = (max = orbit->num / 2, (max/ 2) + (max % 2));

  for(k = 0; k < max; k++)  {
     fprintf(Files.out_data,
          "\n %2d    %2d   %2d %2d %2d   %2d  %3d  %5.2f *",
          k, ptr_mbas[k].orb,  ptr_mbas[k].osc,
             ptr_mbas[k].l,    ptr_mbas[k].par,
             ptr_mbas[k].j,    ptr_mbas[k].m,
             ptr_mbas[k].e);

     if( !((max % 2) && ( k == max - 1)) )   {
        fprintf(Files.out_data,
             "   %2d    %2d   %2d %2d %2d    %2d  %3d  %5.2f",
             k+max, ptr_mbas[k+max].orb,  ptr_mbas[k+max].osc,
                ptr_mbas[k+max].l,    ptr_mbas[k+max].par,
                ptr_mbas[k+max].j,    ptr_mbas[k+max].m,
                ptr_mbas[k+max].e);
     }
  }

    /*
     * Set a ONE in the proper position in Passive_orbit
     * for all orbits which should be excluded from the
     * present calculation.
     */

  for(Passive_orbit = ULL_ZERO, k = 0; k < orbit->num; k++)   {
     if(!ptr_mbas[k].status)  Passive_orbit |= ULL_ONE << k;
  }

  part_orb_limitations(num_j, jbas, orbit);

} /* End: function single_m_states()  */

   /*
   ** The function 
   **               part_orb_limitations()
   ** analyzes all spherical input j-orbits for possible restrictions on
   ** particle occupations. For orbits with restrictions a mask to be used
   ** for the m-orbits is generated and stored in m_orbit->mask_list[]
   ** and upper and lower limits m_orbit->mask_lim[]
   */

void part_orb_limitations(int numj_orb, struct jbas *jbas,
                                          struct m_orbit *m_orbit)
{
  int 
               k, m, osc, l, j, num_mask;
  ULL
               mask;
  struct jbas
               *ptr_jbas;
  struct mbas
               *ptr_mbas;

  ptr_jbas = jbas;                      // initialization

  num_mask = 0;                                     // number of mask orbits
  for(k= 0; k < numj_orb; k++, ptr_jbas++) {
    mask = ULL_ZERO;
    if((ptr_jbas->min_part > 0) || (ptr_jbas->max_part < (ptr_jbas->j + 1)))  {
      m_orbit->lim[2*num_mask]     = ptr_jbas->min_part;     // particle limitations
      m_orbit->lim[2*num_mask + 1] = ptr_jbas->max_part;

      osc  = ptr_jbas->osc;
      l    = ptr_jbas->l;
      j    = ptr_jbas->j;
      for(m = 0, ptr_mbas = m_orbit->mbas; m < m_orbit->num; m++, ptr_mbas++) {
	if((ptr_mbas->osc == osc) && (ptr_mbas->l == l) && (ptr_mbas->j == j)) {
	  mask |= (ULL_ONE << m);
	}
      } // loop through all m_orbits
      m_orbit->list[num_mask] = mask;      // save MASK value
      num_mask++;
    } // if-test for j_orbit reduction

  } // loop through all j_orbits

  if(m_orbit->numMask != num_mask)  {
    printf("\n\nError in function part_orb_limitations():");
    printf("\n Wrong number of mask orbits found ");
    printf("\n sp_bas->num_mask = %d  -- current value = %d\n\n",
            	                       m_orbit->num, num_mask);
    exit(1);
  }

} // End: function  part_orb_limitations()

    /*
     * The function                         
     *        int mbas_parity()                  
     * is a utility function for the library function qsort() in order to
     * sort the single-particle orbits in mbas[] such that orbits with the
     * same m-value and + parity comes before orbits with - parity.
     */

int mbas_parity(const struct mbas *one, const struct mbas *two)
{
  if(one->par> two->par)        return -1;
  else if(one->par < two->par)  return +1;
  else if(one->par == two->par)  {
     if(one->j > two->j)       return -1;
     else                       return +1;
  }    
  else                          return  0;
} /* End: function mbas_parity() */

     /*
      * The function 
      *         read_transitions()
      * allocates memory for a structure trans to
      * store all transition to be calculated.
      */

void read_transitions(void)
{
  char           ch;
  register int    i;

    /*
     * Allocate memory to store information 
     * about the transitions to be calculated 
     */
  
  Model[0].trans = (struct trans *) vector(sizeof(struct trans), Model[0].num_trans);

  for(i = 0; i < Model[0].num_trans; i++)  {
    fscanf(Files.in_data,"%d%d%d%d%c",&(Model[0].trans[i].num_init),
                            &(Model[0].trans[i].j_init),
                            &(Model[0].trans[i].num_final),
                            &(Model[0].trans[i].j_final),
                            &ch);
  }     /* end of loop  */

    /*
     * Re-order the sequence of transitions in 
     * Model.trans[] after increasing initial state.
     */

  qsort(Model[0].trans, (ULL) Model[0].num_trans, sizeof(struct trans),
      (int(*)(const void *, const void *)) trans_comp);

} /* End: read_transitions()  */

    /*
     * The function 
     *                 int trans_comp() 
     * is a utility function for the library function qsort()
     * in order to sort the transitions after increasing
     * number of the initial states. For equal initial states
     * it sorts the final states after increasing order.
     */

int trans_comp(const struct trans *one, const struct trans *two)
{
  if(one->num_init < two->num_init)       return -1;
  else  if(one->num_init > two->num_init) return +1;
  else  {
    if(one->num_final < two->num_final)
      return -1;
    else  if(one->num_final > two->num_final)
      return +1;
    else
      return 0;
  }
} /* End: function trans_comp() */


    /*
     * The function                                           
     *               int plus_parity()                     
     * is a utility function for the library function qsort() 
     * in order to sort the SD after parity such that SD with 
     * positive parity comes before SD with negative parity.  
     */

int plus_parity(const struct sd *one, const struct sd *two)
{
  if(one->parity > two->parity)       return -1;
  else  if(one->parity < two->parity) return +1;
  else                                return  0;
} /* End: function plus_parity() */

    /*
     * The function                                           
     *               int minus_parity()                     
     * is a utility function for the library function qsort() 
     * in order to sort the SD after parity such that SD with 
     * negative parity comes before SD with positive parity.  
     */

int minus_parity(const struct sd *one, const struct sd *two)
{
  if(one->parity < two->parity)       return -1;
  else  if(one->parity > two->parity) return +1;
  else                                return  0;
} /* End: function minus_parity() */

    /*
     * The function                                              
     *               int plus_m_value()                        
     * is a utility function for the library function qsort() in 
     * order to sort the SD after m_value decreasing m-value for 
     * each parity group starting with maximum m-value.          
     */

int plus_m_value(const struct sd *one, const struct sd *two)
{
  if(one->m_value > two->m_value)        return -1;
  else  if(one->m_value < two->m_value)  return +1;
  else                                   return  0;
} /* End: function plus_m_value() */

    /*
     * The function                                              
     *               int minus_m_value()                       
     * is a utility function for the library function qsort() in 
     * order to sort the SD after m_value increasing m-value for 
     * each parity group starting with minimum m-value.          
     */

int minus_m_value(const struct sd *one, const struct sd *two)
{
  if(one->m_value < two->m_value)       return -1;
  else  if(one->m_value > two->m_value) return +1;
  else                                   return  0;
} /* End: function minus_m_value() */

    /*
     * The function                                     
     *               int config_comp()                
     * is a utility function for the library function qsort() in order to 
     * sort the SD after increasing configuration number for each group.
     */

int config_comp(const struct sd *one, const struct sd *two)
{
  if(one->config < two->config)       return -1;
  else  if(one->config > two->config) return +1;
  else                                   return  0;
} /* End: function config_comp() */

  /*
   * The function                       
   *         binomial()              
   * calculates and return the binomial coefficient  " n over k ".         
   */

int binomial(int n, int k)
{
  register int up, down, loop, n_minus_k;

  for(up = 1, down = 1, n_minus_k = n - k, loop = 1; loop <= k; loop++) {
    up   *= n_minus_k + loop;
    down *= loop;
  }
  return (up / down);

} /* End: function binomial() */

    /*
     * The function                           
     *      void run_time(..)                 
     * calculates start and stop time  and  returns the difference.                
     * Input data:                            
     *    int number      - = 1 for start     
     *                        2 for stop      
     *    struct tid run  -  returns the time difference. 
     */

void run_time(int number, struct tid *run)
{
  int    num, stop = 0;
  static int start;

  if(number == 1) {
    start = clock();
    return;
  }
  else if(number == 2)  {
    stop = clock();
  }
  num = ((stop - start) /CLOCK_UNIT);
  run->tick = stop - start;
  run->sec = num % 60;
  run->min = num / 60;
  run->hour = run->min/ 60;
  run->min = run->min % 60;

} /* End: of function run_time(..)   */

    /*
     * The function                                 
     *           clebsch_gordan()                   
     * returns the value of the Clebsch-Gordan      
     * coefficient                                  
     * <j1/2, m1/2, j2/2, m2/2 | j3/2, (m1 + m2)/2> 
     */

double clebsch_gordan(int j1, int j2, int j3, int m1, int m2)
{
  static double fac[150] = {
    .1000000000000000E+01,.3162277660168379E-01,.2000000000000000E-02,
    .1897366596101028E-03,.2400000000000000E-04,.3794733192202056E-05,
    .7200000000000002E-06,.1593787940724864E-06,.4032000000000001E-07,
    .1147527317321902E-07,.3628800000000002E-08,.1262280049054092E-08,
    .4790016000000003E-09,.1969156876524384E-09,.8717829120000006E-10,
    .4135229440701207E-10,.2092278988800002E-10,.1124782407870728E-10,
    .6402373705728006E-11,.3846755834917892E-11,.2432902008176642E-11,
    .1615637450665515E-11,.1124000727777609E-11,.8175125500367506E-12,
    .6204484017332402E-12,.4905075300220504E-12,.4032914611266062E-12,
    .3443362860754794E-12,.3048883446117143E-12,.2796010642932893E-12,
    .2652528598121915E-12,.2600289897927591E-12,.2631308369336940E-12,
    .2745906132211536E-12,.2952327990396046E-12,.3267628297331728E-12,
    .3719933267899019E-12,.4352480892045862E-12,.5230226174666021E-12,
    .6450376682011969E-12,.8159152832478993E-12,.1057861775849963E-11,
    .1405006117752883E-11,.1910498367185033E-11,.2658271574788455E-11,
    .3782786767026367E-11,.5502622159812102E-11,.8178384990311005E-11,
    .1241391559253610E-10,.1923556149721149E-10,.3041409320171346E-10,
    .4905068181788930E-10,.8065817517094409E-10,.1351836790901029E-09,
    .2308436973392420E-09,.4014955268976057E-09,.7109985878048654E-09,
    .1281573721857157E-08,.2350561331282885E-08,.4385545276195193E-08,
    .8320987112741415E-08,.1605109571087441E-07,.3146997326038803E-07,
    .6269557984667544E-07,.1268869321858846E-06,.2608136121621699E-06,
    .5443449390774448E-06,.1153317792981115E-05,.2480035542436839E-05,
    .5411367084667394E-05,.1197857166996993E-04,.2689449441079695E-04,
    .6123445837688631E-04,.1413574626231488E-03,.3307885441519399E-03,
    .7845339175584759E-03,.1885494701666058E-02,.4591092485552201E-02,
    .1132428117820634E-01,.2829031189597267E-01,.7156945704626409E-01,
    .1833212210859029E+00,.4753643337012861E+00,.1247684230710655E+01,
    .3314240134565367E+01,.8908465407274079E+01,.2422709538367284E+02,
    .6665313817722467E+02,.1854826422573992E+03,.5220273782040236E+03,
    .1485715964481768E+04,.4275404227490954E+04,.1243841405464136E+05,
    .3658035857041260E+05,.1087366156656748E+06,.3266626020337846E+06,
    .9916779348709544E+06,.3041882150138602E+07,.9426890448883294E+07,
    .2951234062064472E+08,.9332621544394462E+08,.2980746402685117E+09,
    .9614466715035176E+09,.3131572170660985E+10,.1029901674514568E+11,
    .3419676810361796E+11,.1146280563734714E+12,.3878597438312349E+12,
    .1324641819451836E+13,.4565884904381298E+13,.1588245541522752E+14,
    .5574945468249565E+14,.1974506857221085E+15,.7055650984616650E+15,
    .2543559733472202E+16,.9249958440832429E+16,.3393108684451918E+17,
    .1255404359589777E+18,.4684525849754318E+18,.1762838801735966E+19,
    .6689502913449167E+19,.2559641940120622E+20,.9875044200833661E+20,
    .3840998695345006E+21,.1506141741511150E+22,.5953547977784760E+22,
    .2372173242880062E+23,.9526867474051174E+23,.3856204823625829E+24,
    .1573076357315330E+25,.6466855489220516E+25,.2678949036508007E+26,
    .1118248651196012E+27,.4703162928493458E+27,.1992942746161532E+28,
    .8508021737644667E+28,.3659042881952573E+29,.1585214610157954E+30,
    .6917786472619537E+30,.3040758665204989E+31,.1346201247571762E+32,
    .6002457605114648E+32,.2695364137888182E+33,.1218859041294581E+34,
    .5550293832739345E+34,.2544977678223085E+35,.1174997204390920E+36,
    .5462031093002385E+36,.2556323917872885E+37,.1204487096628886E+38};

  double value;       /* The return value */
  int temp, change;   /* interchange parameters */
  int jz1, jz2, jz3,
     jt1, jt2, jt3,
     jt4, jt5, j4, j5; /* ang. mom. relations */

    /* Local variables */

  int     loop_min, loop_max, phase, loop;
  double                           factor;

  value = 0.0;              /* initialization */

    /*
     * Interchange the angular momenta such that the smallest j-value 
     * is placed in position 2. A parameter is set to remember the interchange 
     *    change  = -1 - j1 is placed in position 2.  
     *            =  0 - no change                       
     *            =  1 - j3 is placed in position 2.  
     */

  if(   (j1 < j2 )  && (j1 < j3)) {

      /* interchange 1 and 2 */

    temp = j1; j1 = j2; j2 = temp;
    temp = m1; m1 = m2; m2 = temp;
    change = -1;
  }
  else if(   (j3 < j1) && (j3 < j2))  {

      /* interchange 2 and 3 */

    temp = j2; j2 = j3; j3 = temp;
    m2 = -(m1 + m2);
    change = 1;
  }
  else change = 0;

       /* Test of angular momentum relations  */

  jz1 = (j1 + j2 - j3)/2;
  jz2 = (j1 + j3 - j2)/2;
  jz3 = (j2 + j3 - j1)/2;

  if((jz1 < 0) || (jz2 < 0) || (jz3 < 0))  return  value ;

       /* Test of angular projection relations */

  if((j1 < abs(m1)) || (j2 < abs(m2)) || (j3 < abs(m1 + m2)))  return  value;

       /* Definition of loop parameters */

  jt1 = (j1 - j3 + m2)/2;
  jt2 = (j2 - j3 - m1)/2;
  jt3 = (j1 - m1)/2;
  jt4 = (j2 + m2)/2;
  jt5 = (j2 - m2)/2;
  j4  = j1/2;
  j5  = j3/2;

        /* Loop limits */

  loop_min = MAX( MAX(jt1, jt2) , 0);
  loop_max = MIN( MIN(jt3, jt4) , jz1);

          /* Loop test */

  if( loop_min > loop_max) return  value;

           /*  test for maximum factorials   */

  if ((jz1 > FAC_LIM) || (jz3 > FAC_LIM) || (jt4 > FAC_LIM) || (jt5 > FAC_LIM)) {
    printf("\n\nTermination from VCC. - too large factorials\n");
    exit(1);
  }
  for (loop = loop_min,phase = PHASE(loop_min); loop <= loop_max;
                              phase *= -1, loop++) {
    value +=  (fac_ratio(jt3-loop,j4) / fac [jt4-loop] )
          * (fac_ratio(loop-jt2,j5) / fac[loop-jt1] )
          * ((phase / fac[jz1-loop]) / fac[loop] );
  }
  factor =  fac_ratio(j4,(j1 + m1)/2) * fac_ratio(j4,jt3)
           * fac_ratio((j1 + j2 + j3)/2+1,jz2)
           * fac_ratio(j5,(j3 + m1 + m2)/2)
           * fac_ratio(j5,(j3 - m1 - m2)/2)
           * fac[jt4] * fac[jt5] * fac[jz1] * fac[jz3]
           * (j3 + 1);
  value *= sqrt(factor);

         /* Add the angular momentum interchange factor */

  if(change == -1) value *= PHASE(jz1);
  else if (change == 1)
    value *= PHASE(jt3) * sqrt(  (j2+1.0)
                      / (j3+1.0));
  return (value);
} /* End: function clebsch_gordan() */

    /*
     * The function               
     *      double fac_ratio()    
     * calculates and returns the ratio (n! / m!).           
     */

double fac_ratio(int m, int n)
{
  int i, high;
  double value;

  value = 1.0;      /* initialization */
  if ( n == m)    return (value);          /* Return for n equal m */
  high = MAX(n,m); 
  for ( i = MIN(n,m)+1; i <= high; i++) value*= i;
  if (n < m) return (1.0/value);       /* Return for n less  m  */
  return (value);                  /* Return for n greater than m */
} /* End: function fac_ratio() */
