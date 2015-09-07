   /******************  The module id-two-part-veffM.c  ******************/

#include "shell-model.h"

      /*
      ** The module entrance function 
      **    id_two_part_intM()
      ** reads from file uncoupled m-scheme two-particle
      ** identical particle matrix elements and store the
      ** result in files pointed to by MAT_OP *op_veff
      */
             /**** data definitions ****/

     /* Storage of the spherical j-j coupled effective
     ** matrix elements for identical particles
     **         <la, lb| V | ra, rb>
     */
 
   typedef struct  {
     int   
  	       la, // left orbit  - particle one
	       lb, // left orbit  - particle two
	       ra, // right orbit - particle one
	       rb; // right orbit - particle two
     ULL
             left, // left two-particle config
	    right; // right part-particle config
                   // config: b<<length + a 
                   // length = 8*(sizeof(ULL) << 1)
     double
            value; //  value of the matrix elemen
   } TWO_PART;

                  /*
                  ** local function declarations in 
		  **      file lanc-input.c
                  */


static int max_nondiag_id_elem(SP_BAS *sp_bas);
     /*
     ** calculates and returns the maximum  possible
     ** number non-diagonal m-scheme two-body 
     ** matrix elements for identical particles.
     */

static int read_id_two_part_matr_elem(char *filename, SP_BAS *sp_bas,
				              TWO_PART **two_part_ref);
     /*
     ** reads uncoupled m-scheme two-particle matrix
     ** elements for identical particles.
     ** Each matrix element is given an lfet and right 
     ** identity number <left_Indexid|V|right_id>
     ** Then they are sorted into monotonic increasing 
     ** group of right_id and for each left_id into
     ** increasing left_id an incr and finally 
     ** stored in TWO_PART two_part()
     ** The function return the total number of matrix
     ** elements stored.
     */
static void id_matr_elem_identity(SP_BAS *sp_bas, int *int_data,
				  TWO_PART *two_part);
     /*
     ** convert (n,l.j) input data to orbit numbers
     ** and store the result in stroct two_part
     */

static void id_orbit_interchange(TWO_PART *two_part);
     /*
     ** reorder if necessary two_part.la, two_part.lb, 
     ** two_part.ra and two_part.b such that    
     **    two_part.la <= two_part.lb and  
     **    two_part.ra <= two_part.rb
     ** If any change multiply matrix element with 
     ** a Pauli phase factor 
     ** Calculate and store the config two_part.left
     ** and two_part.right such that 
     **    two_part.left >=  two_part.right
     */

static void sort_matrix_elements(int num_elem, TWO_PART *two_part);
     /*
     ** sorts all matrix elements  in struct two_part()
     ** after increasing two_part()->right. 
     ** Next, for each group of two_part()->right sort the matrix 
     ** elements after increasing two_part()->left. 
      */

static int right_config(const TWO_PART *one, const TWO_PART *two);
     /*
     ** is a utility function for the library function qsort() in order to sort the
     ** two-body matrix elements into groups of increasing value of two_part->right.
     */

static int left_config(const TWO_PART *one, const TWO_PART *two);
     /*
     ** is a utility function for the library function qsort() in order to sort the
     ** two-body matrix elements for the same two_part->right 
     ** into sunbgroups of increasing value of two_part->left.
     */
static int id_compress_matr_elem(int num, TWO_PART *two_part);
     /*
     ** takes num list of sorted two-particle matrix 
     ** element and removes possibly equal elements 
     ** Note: The function must have num > 1;
     ** The final number of elements is returned.
     */

static void check_id_matr_elem(SP_BAS *sp_bas, int num_elem, TWO_PART *two_part);
     /*
     ** checks current two-particle matrix elements to see if 
     ** there are any elements missing compared to the available
     ** single-particle orbits in the model.
     */

static void id_m_pot_diag(SP_BAS *sp_bas, MATR_OP *op_veff, 
                          int num_elem, TWO_PART *two_part);
     /*
     ** transfers all m-scheme diagonal two-particle matrix
     ** elements to the struct MATR_OP op_veff->diag
     */

static int id_m_pot_nondiag(SP_BAS *sp_bas, MATR_OP *op_veff, TWO_PART *two_part);
     /*
     ** transfer all nondiagonal m-scheme two-particle matrix elements
     ** <i,j| VEFF |k.l> with (i,j) > (k,l) to op_veff->id_nondiag_bas[],
     ** one group of elements for each (k,l).
     ** op_veff->id_nondiag_table[((2*num_m - k - 3) * k)/2 + (l - 1)] 
     ** element points to the first element of the (k,l) group.
     */

                /**** End: function declarations ****/


               /**** The function definitions  ****/ 

     /*
     ** The entrance function 
     **         id_two_part_intM()
     ** reads from file uncoupled m-scheme identical 
     ** two-particle matrix elements and store the
     ** result in tables: ID_TABLE id_table[]
     */
     

void id_two_part_intM(char *filename, SP_BAS *sp_bas, MATR_OP *op_veff)
{
  char     *func = {"id_two_part_veffM(): "};
  int      num_diag_elem, num_nondiag_elem, tot_num_elem, num;
  TWO_PART *two_part;

             // memory for  kl-table

   num_diag_elem = (sp_bas->numm_orb * (sp_bas->numm_orb - 1)) >> 1;
   op_veff->id_diag
       = MALLOC(num_diag_elem, double, func, "id_kl_diag[]");  
   op_veff->id_nondiag_table
       = MALLOC(num_diag_elem, ID_INT *, func,"id_kl_nondiag_table[]");

   num_nondiag_elem // max number of m-scheme matrix elements 
       = max_nondiag_id_elem(sp_bas);
    op_veff->id_nondiag_bas 
       = MALLOC(num_nondiag_elem, ID_INT, func,"id_nondiag_int_elements[]");
 
   // read and store two-particle matrix elements in M-scheme

   tot_num_elem = read_id_two_part_matr_elem(filename, sp_bas, &two_part);

     /*
     ** store diag m-scheme two-particle matrix
     ** elements in MATR_OP op_veff.id_diag
     */

  id_m_pot_diag(sp_bas, op_veff, tot_num_elem, two_part);

  // calculate and store nondiag m-scheme two-particle matrix elements

   num = id_m_pot_nondiag(sp_bas,op_veff, two_part); 

   if(num > num_nondiag_elem) {
     printf("\n\nError in function id_effective_interaction(): ");
     printf(" = %d - calculated = %d",num_nondiag_elem, num); 
     printf("\n This may not occur!!!!!\n");
     exit(1);
   }

   if(num < num_nondiag_elem) { // reduce size of nondiag[] to new number
      op_veff->id_nondiag_bas = REALLOC( op_veff->id_nondiag_bas, num, ID_INT, func,
                                                     "new id_nondiag_bas[]");
   }

   free(two_part);  // release temporary memory

} // End: function  id_two_part_veffM()

  /*
  ** The function                                           
  **           max_nondiag_id_elem()
  ** calculates and returns the maximum  possible
  ** number non-diagonal m-scheme two-body 
  ** matrix elements for identical particles.
  ** Formula:  numTwoPartConfig*(numTwoPartConfig-1)/2
  */
static int max_nondiag_id_elem(SP_BAS *sp_bas)
{
  int    left_m, left_par, leftMatr[3],
         number, right_m, right_par, rightMatr[3];

  number = 0;            // initialization

  leftMatr[0] = 0; // the lowest configuration
  leftMatr[1] = 1;
  leftMatr[2] = sp_bas->numm_orb;

  //  loop through  all |left> config

  for( ; ; ) {
    left_m   = sp_bas->mbas[leftMatr[0]].m   + sp_bas->mbas[leftMatr[1]].m;
    left_par = sp_bas->mbas[leftMatr[0]].par * sp_bas->mbas[leftMatr[1]].par;

    // loop through all right config above current left config
    
    rightMatr[0] = leftMatr[0];     // initialization |right> = |left>
    rightMatr[1] = leftMatr[1];
    rightMatr[2] = leftMatr[2];

    //  loop through  all |right> config

    for( ; ; ) {
      right_m   = sp_bas->mbas[rightMatr[0]].m   + sp_bas->mbas[rightMatr[1]].m;
      right_par = sp_bas->mbas[rightMatr[0]].par * sp_bas->mbas[rightMatr[1]].par;

     if((left_m == right_m) && (left_par == right_par)) number++;

      // new right

      if((rightMatr[0]+1) < rightMatr[1]) {
	rightMatr[0]++;
      }
      else if((rightMatr[1]+1) < rightMatr[2])  {
	rightMatr[1]++;
	rightMatr[0] = 0;
      }
      else break;
      
    } // end through |right>


    // next |left>

    if((leftMatr[0]+1) < leftMatr[1]) {
	leftMatr[0]++;
      }
      else if((leftMatr[1]+1) < leftMatr[2])  {
	leftMatr[1]++;
	leftMatr[0] = 0;
      }
    else break;
    }  // end through |left>

  return number;
  
} // End: function max_nondiag_id_elem()

     /*
     ** The function
     **         read_id_two_part_matr_elem()
     ** reads uncoupled m-scheme two-particle matrix
     ** elements for identical particles.
     ** Each matrix element is given an lfet and right 
     ** identity number <left_id|V|right_id>
     ** Then they are sorted into monotonic increasing 
     ** group of right_id and for each left_id into
     ** increasing left_id an incr and finally 
     ** stored in TWO_PART two_part()
     ** The function return the total number of matrix
     ** elements stored.
     */

static int read_id_two_part_matr_elem(char *filename, SP_BAS *sp_bas,
                                               TWO_PART **two_part_ref)
{
  char      ch, *func = {"read_id_two_part_matr_elem(): "};
  int       loop, int_data[16], num_elem, new_num_elem;
  TWO_PART  *two_part, *twoPartPtr;
  FILE      *file_ptr;

  if( (file_ptr = fopen(filename,"r")) == NULL) {
    printf("\n\nError in function read_id_two_part_matr_elem();");
    printf("\nWrong file = %s for input of eff.interaction.\n",filename);
    exit(1);
  }

  // read total number of matrix elements

  if((fscanf(file_ptr,"%d%c",&num_elem,&ch)!=2)|(ch!='\n')) {
    printf("\n\n Error in function read_two_part_matr_elem(): ");
    printf("\nFirst number in file %s", filename);
    printf(" - number of matrix element - is not correct\n");
    exit(1);
  } 
  // memory to store two-particle matrix elements

  two_part = MALLOC(num_elem, TWO_PART, func,"two_part[]");

  new_num_elem = 0;       // initialization
  twoPartPtr   = two_part; 

  for(loop = 0; loop < num_elem; loop++) {

     // read (n1,l1,j1,m1), (n2,l2,j2.m2), (n3,l3,j3,m3), (n4,l4,j4,m4) 

    if(read_int_number(file_ptr, 16, int_data) == FALSE) {
      printf("\n\nError in function read_two_part_matr_elem();");
      printf("\nError in (n,l,j) - values");
      printf(" for matrix element number = %d\n",loop);
      exit(1);
    }
    // read effective matrix element

    if(read_float_number(file_ptr, 1, &twoPartPtr->value) == FALSE) {
      printf("\n\nError in function read_two_part_matr_elem();");
      printf("\nWrong value for two-particle matrix element no %d",loop);
      printf("\n read from file %s\n", filename);
      exit(1);
    }
    // test conservation of total M-value and parity

    if(  ((int_data[3] + int_data[7]) != (int_data[11] + int_data[15]))
       ||(  ((int_data[1] + int_data[5]) % 2) 
          !=((int_data[9] + int_data[13])% 2))) {
      printf("\n\nError in function read_id_two_part_matr_elem():");
      printf("\n  M-value and parity is wrong for matrix elemen no %d\n", loop);
      exit(1);
    }
       /*
       ** convert (n,l.j, m) input data to orbit numbers
       ** and store the result in struct two_part
       */
      
    id_matr_elem_identity(sp_bas, int_data, twoPartPtr);

       /*
       **  Only matrix elements relevant for the present
       **  set of single-particle orbits are saved
       */

    if(  (twoPartPtr->la == -1) || (twoPartPtr->lb == -1)
       ||(twoPartPtr->ra == -1) || (twoPartPtr->rb == -1)) continue;
                  
     /*
     ** reorder if necessary twoPartPtr.la, twoPartPtr.lb, 
     ** twoPartPtr.ra and twoPartPtr.b such that    
     **       twoPartPtr.la <= twoPartPtr.lb and  
     **       twoPartPtr.ra <= twoPartPtr.rb
     ** Add a Pauli phase factor to matrix element.
     ** Calculate and store the config twoPartPtr.left
     ** and twoPartPtr.right such that 
     **    twoPartPtr.left >=  twoPartPtr.right
     */
 
    id_orbit_interchange(twoPartPtr);

    twoPartPtr++;            // current matrix elements accepted
    new_num_elem++;
    
   } // end loop through all input two-particle matrix elements 

   fclose(file_ptr);

      /* 
      ** All matrix elements are sorted into increasing groups of 
      ** same two_part->right and for each group into subgroups
      ** after increasing two_part->left.
      */

   sort_matrix_elements(new_num_elem, two_part);

   // remove possibly equal matrix elements

   new_num_elem = id_compress_matr_elem(new_num_elem, two_part);

   if(new_num_elem < num_elem) {  // reduce size of two_part[]
     two_part = REALLOC(two_part, new_num_elem, TWO_PART, func,
                                      "realloc j_pot.two_part[]");
   }

   // check that all necessary two-particle matrix elements are stored
 
   check_id_matr_elem(sp_bas, new_num_elem, two_part);

   *two_part_ref = two_part; // return storage address to matrix elements
   
   return new_num_elem;

} // End: function  read_id_two_part_matr_elem()

     /*
     ** The function 
     **         id_matr_elem_identity()
     ** convert (n,l.j) input data to orbit numbers
     ** and store the result in stroct two_part
     */

static void id_matr_elem_identity(SP_BAS *sp_bas, int *int_data,
                                              TWO_PART *two_part)
{
  int        k;
  MBAS       *basPtr;

  two_part->la = -1;     // initialization
  two_part->lb = -1;
  two_part->ra = -1;
  two_part->rb = -1;
  basPtr       = sp_bas->mbas;

  for(k = 0; k < sp_bas->numm_orb; k++, basPtr++)     {
    if(   (basPtr->osc == int_data[0]) && (basPtr->l == int_data[1]) 
       && (basPtr->j   == int_data[2]) && (basPtr->m == int_data[3]))
      two_part->la = k;
    if(   (basPtr->osc == int_data[4]) && (basPtr->l == int_data[5]) 
       && (basPtr->j   == int_data[6]) && (basPtr->m == int_data[7]))
      two_part->lb = k;
    if(   (basPtr->osc == int_data[8]) && (basPtr->l == int_data[9]) 
       && (basPtr->j   == int_data[10])&& (basPtr->m == int_data[11]))
      two_part->ra = k;
    if(   (basPtr->osc == int_data[12]) && (basPtr->l == int_data[13]) 
       && (basPtr->j   == int_data[14]) && (basPtr->m == int_data[15]))
      two_part->rb = k;
  }

            /* "Whitehead" phase is included */
 
  two_part->value *=  sp_bas->mbas[two_part->la].phase
                    * sp_bas->mbas[two_part->lb].phase
                    * sp_bas->mbas[two_part->ra].phase
                    * sp_bas->mbas[two_part->rb].phase;
                
} // End function  id_matr_elem_identity()

     /*
     ** The function 
     **         id_orbit_interchange()
     ** reorder if necessary two_part.la, two_part.lb, 
     ** two_part.ra and two_part.b such that    
     **    two_part.la <= two_part.lb and  
     **    two_part.ra <= two_part.rb
     ** If any change multiply matrix element with 
     ** a Pauli phase factor 
     ** Calculate and store the config two_part.left
     ** and two_part.right such that 
     **    two_part.left >=  two_part.right
     */

static void id_orbit_interchange(TWO_PART *two_part)
{
  int        temp, phase;
  ULL        length, left, right;

  phase  = +1;               // initialization
  length = 8*(sizeof(ULL) >> 1);
 
  if(two_part->la > two_part->lb) {
    temp         = two_part->la;
    two_part->la =  two_part->lb;
    two_part->lb = temp;
    phase *= -1;
  }
  if(two_part->ra > two_part->rb) {
    temp         = two_part->ra;
    two_part->ra =  two_part->rb;
    two_part->rb = temp;
    phase *= -1;
  }
  two_part->value *= (double)phase;

       /* if necessary, interchange right and
       ** left orbit identity such that
       ** two_part->left < two_part->right
       */

  left  = (((ULL)two_part->lb) << length) + ((ULL)two_part->la);
  right = (((ULL)two_part->rb) << length) + ((ULL)two_part->ra);

  if(left < right) {
    temp         = two_part->la;
    two_part->la =  two_part->ra;
    two_part->ra = temp;
    temp         = two_part->lb;
    two_part->lb =  two_part->rb;
    two_part->rb = temp;
  }

  two_part->left  = (((ULL)two_part->lb) << length) + ((ULL)two_part->la);
  two_part->right = (((ULL)two_part->ra) << length) + ((ULL)two_part->rb);

} // End: function id_matr_elem_identity()

     /* 
     ** The function
     **        sort_matrix_elements()             
     ** sorts all matrix elements  in struct two_part()
     ** after increasing two_part()->right. 
     ** Next, for each group of two_part()->right sort the matrix 
     ** elements after increasing two_part()->left. 
     */

static void sort_matrix_elements(int num_elem, TWO_PART *two_part)
{
   int       loop, num;
   ULL       test;
   TWO_PART  *ptr_start, *ptr;

  // sort matrix elements into groups after increasing two_part->right

   qsort(two_part, (size_t)num_elem, (size_t)sizeof(TWO_PART),
                      (int(*)(const void *, const void *)) right_config);
   /*
   ** sort matrix elements with same two_part->right
   ** into increasing two_part->left index
   */

   ptr_start = two_part;    // initialization;
   ptr       = ptr_start;
   test      = ptr_start->right;
   num       = 0;

   for(loop = 0; loop < num_elem; loop++) {
     if(ptr->right == test) {
       num++;
       ptr++;
       continue;
     }
     if(num > 0)  {
       qsort(ptr_start, (size_t)num, (size_t)sizeof(TWO_PART),
            (int(*)(const void *, const void *)) left_config);
       ptr_start += num;
       ptr        = ptr_start;
       test       = ptr_start->right;
       ptr++;
       num        = 1;
     }
     else {
       ptr_start += num;
       ptr        = ptr_start;
       test       = ptr_start->right;
       num        = 0;
     }
   } // end loop;

   if(num > 0)  {        // last block
     qsort(ptr_start, (size_t)num, (size_t)sizeof(TWO_PART),
	    (int(*)(const void *, const void *)) left_config);
   }

} /* End: function sort_matrix_elements() */

     /*
     ** The function                         
     **        int right_config()                  
     ** is a utility function for the library function qsort() in order to sort the
     ** two-body matrix elements into groups of increasing value of two_part->right.
     */

static int right_config(const TWO_PART *one, const TWO_PART *two)
{
  if(one->right > two->right)       return +1;
  else  if(one->right < two->right) return -1;
  else                              return  0;
} // End: function right_config()

     /*
     ** The function                         
     **        int left_config()                  
     ** is a utility function for the library function qsort() in order to sort the
     ** two-body matrix elements for the same two_part->right 
     ** into sunbgroups of increasing value of two_part->left.
     */

static int left_config(const TWO_PART *one, const TWO_PART *two)
{
  if(one->left > two->left)       return +1;
  else  if(one->left < two->left) return -1;
  else                            return  0;     
} // End: function left_config()


      /*
      ** The function 
      **      id_compress_matr_elem()
      ** takes num list of sorted two-particle matrix 
      ** element and removes possibly equal elements 
      ** Note: The function must have num > 1;
      ** The final number of elements is returned.
      */

static int id_compress_matr_elem(int num, TWO_PART *two_part)
{
  int   k, m, limit;

  limit = num - 1;  // initialization
  k     = 0;

  while(k < limit) {
    if(  (two_part[k].left  == two_part[k+1].left)
       &&(two_part[k].right == two_part[k+1].right))  {
      for(m = k+1; m < limit; m++) {
	two_part[m].la    = two_part[m+1].la;
	two_part[m].lb    = two_part[m+1].lb;
	two_part[m].ra    = two_part[m+1].ra;
	two_part[m].rb    = two_part[m+1].rb;
	two_part[m].left  = two_part[m+1].left;
	two_part[m].right = two_part[m+1].right;
	two_part[m].value = two_part[m+1].value;
      }
      limit--;
    }
    else {
      k++;
    }
  } // end while() loop

  return (limit + 1);
   
} /* End: function id_compress_matr_elem() */

     /*
     ** The function 
     **        check_id_matr_elem()
     ** checks current two-particle matrix elements to see if 
     ** there are any elements missing compared to the available
     ** single-particle orbits in the model.
     */

static void check_id_matr_elem(SP_BAS *sp_bas, int num_elem, TWO_PART *two_part)
{
   int      length, k, l, i, j, m_kl, parity_kl, number, limit;
   ULL      left, right;

   number = 0;                        // initialization
   limit  =  sp_bas->numm_orb - 1;
   length = 8*(sizeof(ULL) >> 1);
 
   for(k = 0; k < limit; k++) {  // first right 
     for(l = k+1; l < sp_bas->numm_orb; l++) {  // second right 
       m_kl      = sp_bas->mbas[k].m + sp_bas->mbas[l].m;
       parity_kl = sp_bas->mbas[k].par * sp_bas->mbas[l].par;
       right = (((ULL) k) << length) + l;
       left  =  (((ULL) l) << length) + k;
       //  diagonal  matrix  element

       if(  (left  == two_part[number].left)
	    &&(right == two_part[number].right))  {
       }
       else {
	 printf("\n\nError in function check_id_matr_elem()");
	 printf("\n matrix element <%d, %d|V|%d, %d> is not found",k,l,k,l);
	 exit(1);
       }
       number++;

       for(j = l; j < sp_bas->numm_orb; j++)  { // second left
	 for(i = (j == l) ? k + 1: 0; i < j; i++){ // first left

	   // test mj and parity

	   if(   (m_kl == sp_bas->mbas[i].m + sp_bas->mbas[j].m)
	      && (parity_kl== sp_bas->mbas[i].par * sp_bas->mbas[j].par))  {
	     left = (((ULL) j) << length) + i;
	     if(  (left  == two_part[number].left)
		&&(right == two_part[number].right))  {
               number++;
	       continue;
	     }
	     else {
	       printf("\n\nError in function check_id_matr_elem()");
	       printf("\n matrix element <%d, %d|V|%d, %d> is not found",i,j,k,l);
	       exit(1);
	     }
	   }
	 } // end i orb
       } // end j orb
     } // end l orb
   } // end k orb

   if(number != num_elem) {
     printf("\n\nError in function check_id_matr_elem(): ");
     printf("Not enough two-particle matrix elements input file");
     printf("\nRead number = %d - should have been = %d\n", number, num_elem);
     exit(1);
   }
 
} // End: function check_id_matr_elem()

     /*
     ** The function                                          
     **              id_m_pot_diag(...)                    
     ** transfers all m-scheme diagonal two-particle matrix
     ** elements to the struct MATR_OP op_veff->diag
     */

static void id_m_pot_diag(SP_BAS *sp_bas, MATR_OP *op_veff, 
                      int tot_num_elem, TWO_PART *two_part)
{
  int         k, l, limit, step;
  ULL         length, left, right;
  double      *id_diag;
  TWO_PART    *ptr;

  ptr     = two_part;   // initialization
  id_diag = op_veff->id_diag;
  length  = 8*(sizeof(ULL) >> 1);
  step    = 0;
  limit   = sp_bas->numm_orb - 1;
   
  for(k = 0; k < limit; k++) {
    for(l = k+1; l <= limit; l++) {
      right = (((ULL)k) << length) + (ULL)l;
      left  = (((ULL)l) << length) + (ULL)k;

      while( !((right == ptr->right) && (left == ptr->left)))  {
	ptr++;
	if(++step >= tot_num_elem) {
	  printf("\n\n Error in function id_m_pot_diag():");
	  printf("The m-scheme diagonal two-particle matrix ");
	  printf("element <%d,%d|V|%d,%d> is missing",k,l,k,l);
	  exit(1);
	}
      }
      *(id_diag++) = ptr->value;
    } // end l index
  } // end k index

} /* End: function id_m_pot_diag() */

     /*
     ** The function                                          
     **             id_m_pot_nondiag(,...)                    
     ** transfer all nondiagonal m-scheme two-particle matrix elements
     ** <i,j| VEFF |k.l> with (i,j) > (k,l) to op_veff->id_nondiag_bas[],
     ** one group of elements for each (k,l).
     ** op_veff->id_nondiag_table[((2*num_m - k - 3) * k)/2 + (l - 1)] 
     ** element points to the first element of the (k,l) group.
     */

static int id_m_pot_nondiag(SP_BAS *sp_bas, MATR_OP *op_veff, TWO_PART *two_part)
{
  int       i, j, k, l, m_kl, parity_kl, limit, count, number;
  ULL       length, leftIndex, rightIndex;
  ID_INT    *nondiag, **nondiag_table = NULL;
  TWO_PART  *ptr;

  nondiag_table = op_veff->id_nondiag_table;
  op_veff-> maxNondiagElem = 0;

  limit   =  sp_bas->numm_orb - 1;  // initialization
  nondiag = op_veff->id_nondiag_bas;
  number  = 0;
  length  = 8*(sizeof(ULL) >> 1);
 
  for(k = 0; k < limit; k++) {   // first right
    for(l = k+1; l <= limit; l++, nondiag_table++)  { // second right
      m_kl = sp_bas->mbas[k].m + sp_bas->mbas[l].m;
      parity_kl = sp_bas->mbas[k].par * sp_bas->mbas[l].par;

      *nondiag_table = nondiag; // store pointer

      rightIndex = (((ULL)k) << length) + (ULL)l;
      ptr = two_part;
      while(ptr->right != rightIndex) ptr++;

      for(count = 0,j = l; j <= limit; j++)  {    // second left
	for(i = (j == l) ? k + 1: 0; i < j; i++){ // first left

	  // test mj and parity

	  if(   (m_kl == sp_bas->mbas[i].m + sp_bas->mbas[j].m)
	     && (parity_kl == sp_bas->mbas[i].par * sp_bas->mbas[j].par))  {
	    ptr++;
	    leftIndex = (((ULL)j) << length) + (ULL)i;
	    if(ptr->left != leftIndex) {
	      printf("\n\nError in function id_m_pot_nondiag():");
	      printf("\n no nondiag matrix element matching");
	      printf("\n two_part->left = %d - ij_left = %d\n",
                                     (int)ptr->left, (int)leftIndex);
	      exit(1);
	    }
	    nondiag->val = ptr->value;

	    if(fabs(nondiag->val) > MATRIX_LIMIT) {
	      nondiag->one    = ULL_ONE<<j ^ ULL_ONE<<i;
	      nondiag->two    = (ULL_ONE<<j) - (ULL_ONE<<(i+1));
	      count++;
	      number++;  
	      nondiag++;
	    }
	  } // end mj and parity test
		     	 	    
	}  // end index i
      } // end index j
      if(!count)  {// no elements for current (k,l)
	*nondiag_table = NULL_PTR;
      }
      else   { // count elements for current (k,l)
	op_veff-> maxNondiagElem = MAX(op_veff-> maxNondiagElem, count);
	number++; 
	(nondiag++)->one = ULL_ZERO;// terminate each (ij)-set with a ZERO
      }
    } // end of index l
  } // end of  index k

  return number; // return total number of elements

} // End: function id_m_pot_nondiag()
