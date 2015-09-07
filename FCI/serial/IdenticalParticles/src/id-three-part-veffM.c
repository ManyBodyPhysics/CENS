   /******************  The module id-two-part-veffM.c  ******************/

#include "shell-model.h"

      /*
      ** The module entrance function 
      **    id_three_part_veffM()
      ** reads from file uncoupled m-scheme three-particle
      ** identical particle matrix elements and store the
      ** result in files pointed to by MAT_OP *op_veff
      */
             /**** local data definitions ****/

     /* Storage of the spherical m-scheme three-particle 
     ** effective matrix elements for identical particles
     **         <la, lb, lc| V | ra, rb, rc>
     */
 
   typedef struct  {
     int   
  	       la, // left orbit  - particle one
	       lb, // left orbit  - particle two
	       lc, // left orbit  - particle three
	       ra, // right orbit - particle one
	       rb, // right orbit - particle two
	       rc; // right orbit - particle three
     ULL
             left, // left three-particle config
	    right; // right three-particle config
     double
            value; // value of threePart interaction
   } THREE_PART;

                  /*
                  ** local function declarations in 
		  **      file lanc-input.c
                  */


static int max_nondiag_id_elem(SP_BAS *sp_bas);
     /*
     ** calculates and returns the maximum  possible
     ** number non-diagonal m-scheme three-body 
     ** matrix elements for identical particles.
     */

static int read_id_three_part_matr_elem(SHELL *model, SP_BAS *sp_bas,
                                         THREE_PART **three_part_ref);
     /*
     ** reads uncoupled m-scheme three-particle matrix
     ** elements for identical particles.
     ** Each matrix element is given an lfet and right 
     ** identity number <left_Indexid|V|right_id>
     ** Then they are sorted into monotonic increasing 
     ** group of right_id and for each left_id into
     ** increasing left_id an incr and finally 
     ** stored in THREE_PART three_part()
     ** The function return the total number of matrix
     ** elements stored.
     */
static void id_matr_elem_identity(SP_BAS *sp_bas, int *int_data,
				  THREE_PART *three_part);
     /*
     ** convert (n,l.j) input data to orbit numbers
     ** and store the result in stroct three_part
     */

static void id_orbit_interchange(THREE_PART *three_part);
     /*
     ** reorder if necessary three_part.la, three_part.lb, 
     ** three_part.ra and three_part.b such that    
     **    three_part.la <= three_part.lb and  
     **    three_part.ra <= three_part.rb
     ** If any change multiply matrix element with 
     ** a Pauli phase factor 
     ** Calculate and store the config three_part.left
     ** and three_part.right such that 
     **    three_part.left >=  three_part.right
     */

static void sort_matrix_elements(int num_elem, THREE_PART *three_part);
     /*
     ** sorts all matrix elements  in struct three_part()
     ** after increasing three_part()->right. 
     ** Next, for each group of three_part()->right sort the matrix 
     ** elements after increasing three_part()->left. 
      */

static int right_config(const THREE_PART *one, const THREE_PART *two);
     /*
     ** is a utility function for the library function qsort() in order to sort the
     ** two-body matrix elements into groups of increasing value of three_part->right.
     */

static int left_config(const THREE_PART *one, const THREE_PART *two);
     /*
     ** is a utility function for the library function qsort() in order to sort the
     ** two-body matrix elements for the same three_part->right 
     ** into sunbgroups of increasing value of three_part->left.
     */
static int id_compress_matr_elem(int num, THREE_PART *three_part);
     /*
     ** takes num list of sorted three-particle matrix 
     ** element and removes possibly equal elements 
     ** Note: The function must have num > 1;
     ** The final number of elements is returned.
     */

static void id_m_pot_diag(SP_BAS *sp_bas, MATR_OP *op_veff, 
                                    int num_elem, THREE_PART *three_part);
     /*
     ** transfers all m-scheme diagonal three-particle matrix
     ** elements to the struct MATR_OP op_veff->diag
     */

static int id_m_pot_nondiag(SP_BAS *sp_bas, MATR_OP *op_veff,
                                                   THREE_PART *three_part);
     /*
     ** transfer all nondiagonal m-scheme three-particle matrix elements
     ** <h,i,j| VEFF |k.l,m> with (h,i,j) > (k,l,m) to op_veff->id_nondiag_bas[],
     ** one group of elements for each (k,l.m).
     ** op_veff->id_nondiag_table[] 
     ** element points to the first element of the (k,l.m) group.
     */

                /**** End: function declarations ****/


               /**** The function definitions  ****/ 

     /*
     ** The entrance function 
     **         id_three_part_veffM()
     ** reads from file uncoupled m-scheme identical 
     ** three-particle matrix elements and store the
     ** result in tables: ID_TABLE id_table[]
     */
     
void id_three_part_veffM(SHELL *model, SP_BAS *sp_bas, MATR_OP *op_veff)
{
  char       *func = {"id_three_part_veffM(): "};
  int        num_diag_elem, num_nondiag_elem, tot_num_elem, num;
  THREE_PART *three_part;

  // memory for klm-table

  num_diag_elem = (sp_bas->numm_orb*(sp_bas->numm_orb - 1) 
                  *(sp_bas->numm_orb - 2) ) / 6;
  op_veff->id_diag
       = MALLOC(num_diag_elem, double, func, "id_klm_diag[]");  
  op_veff->id_nondiag_table
       = MALLOC(num_diag_elem, ID_INT *, func,"id_klm_nondiag_table[]");

  num_nondiag_elem  // max number of m-scheme matrix elements 
       = max_nondiag_id_elem(sp_bas);
  op_veff->id_nondiag_bas 
       = MALLOC(num_nondiag_elem, ID_INT, func,"id_nondiag_int_elements[]");
 
  // read and store three-particle matrix elements in M-scheme

  tot_num_elem = read_id_three_part_matr_elem(model, sp_bas, &three_part);

     /*
     ** store diag m-scheme three-particle matrix
     ** elements in MATR_OP op_veff.id_diag
     */

  id_m_pot_diag(sp_bas, op_veff, tot_num_elem, three_part);

  num = id_m_pot_nondiag(sp_bas,op_veff, three_part); 

  if(num == 0) {
    free(op_veff->id_nondiag_bas);
    op_veff->id_nondiag_bas = NULL;
  }
  else {
    if(num > num_nondiag_elem) {
      printf("\n\nError in function id_three_part_veffM(): ");
      printf(" = %d - calculated = %d",num_nondiag_elem, num); 
      printf("\n This may not occur!!!!!\n");
      exit(1);
    }

    if(num < num_nondiag_elem) { // reduce size of nondiag[] to new number
      op_veff->id_nondiag_bas = REALLOC( op_veff->id_nondiag_bas, num, ID_INT, func,
                                                     "new id_nondiag_bas[]");
    }
  } 
  free(three_part);  // release temporary memory

} // End: function  id_three_part_veffM()

  /*
  ** The function                                           
  **           max_nondiag_id_elem()
  ** calculates and returns the maximum  possible
  ** number non-diagonal m-scheme three-particle 
  ** matrix elements for identical particles.
  **
  */
static int max_nondiag_id_elem(SP_BAS *sp_bas)
{
  int     left_m, right_m,  left_par, right_par, 
          leftMatr[3], rightMatr[3], countRight, number;

  number= 0;  // initialization

  leftMatr[0] = 0;   // lowest |left>
  leftMatr[1] = 1;
  leftMatr[2] = 2;

  do  {   // run through all |left>
    left_m   =  sp_bas->mbas[leftMatr[0]].m 
               +sp_bas->mbas[leftMatr[1]].m
               +sp_bas->mbas[leftMatr[2]].m;
    left_par =  sp_bas->mbas[leftMatr[0]].par 
               *sp_bas->mbas[leftMatr[1]].par
               *sp_bas->mbas[leftMatr[2]].par;

    rightMatr[0] = leftMatr[0];  // |right> = |left>  
    rightMatr[1] = leftMatr[1];
    rightMatr[2] = leftMatr[2];

   // next |right>
      
    if((rightMatr[0] + 1) < rightMatr[1]) {
      rightMatr[0]++;
    }
    else if((rightMatr[1] + 1) < rightMatr[2]) {
      rightMatr[1]++;
      rightMatr[0] = 0;
    }
    else {
      rightMatr[2]++;
      rightMatr[1] = 1;
      rightMatr[0] = 0;
    }
    if(rightMatr[2] <= sp_bas->numm_orb)  {

      // run through all |right> 

      countRight = 0;   // initialization

      do  {   // run through all |right>  for current |left>
	right_m   =  sp_bas->mbas[rightMatr[0]].m 
                    +sp_bas->mbas[rightMatr[1]].m
	            +sp_bas->mbas[rightMatr[2]].m;
	right_par =  sp_bas->mbas[rightMatr[0]].par 
                    *sp_bas->mbas[rightMatr[1]].par
                    *sp_bas->mbas[rightMatr[2]].par;
      
	//  test |left> and |right> for m- and parity-values
 
	if((right_m == left_m) && (right_par == left_par))  {
	  countRight++;
	  number++;
	}
	if(countRight > 0) number++; // extra NULL elements

	// next |right>
      
	if((rightMatr[0] + 1) < rightMatr[1]) {
	  rightMatr[0]++;
	}
	else if((rightMatr[1] + 1) < rightMatr[2]) {
	  rightMatr[1]++;
	  rightMatr[0] = 0;
	}
	else {
	  rightMatr[2]++;
	  rightMatr[1] = 1;
	  rightMatr[0] = 0;
	}
      } while(rightMatr[2] <= sp_bas->numm_orb); // loop |right>

    } // no more |right>

    // new |left>
      
    if((leftMatr[0] + 1) < leftMatr[1]) {
      leftMatr[0]++;
    }
    else if((leftMatr[1] + 1) < leftMatr[2]) {
      leftMatr[1]++;
      leftMatr[0] = 0;
    }
    else {
      leftMatr[2]++;
      leftMatr[1] = 1;
      leftMatr[0] = 0;
    }
  } while(leftMatr[2] <= sp_bas->numm_orb); // |left>

  return number;

} // End: function   max_nondiag_id_elem()

     /*
     ** The function
     **         read_id_three_part_matr_elem()
     ** reads uncoupled m-scheme three-particle matrix
     ** elements for identical particles.
     ** Each matrix element is given an left and right 
     ** identity number <left_id|V|right_id>
     ** Then they are sorted into monotonic increasing 
     ** group of right_id. For each right_id sorted
     ** into increasing left_id and finally 
     ** stored in THREE_PART three_part()
     ** The function return the total number of matrix
     ** elements stored.
     */

static int read_id_three_part_matr_elem(SHELL *model,SP_BAS *sp_bas,
                                          THREE_PART **three_part_ref)
{
  char        ch, *func = {"read_id_three_part_matr_elem(): "};
  int         loop, int_data[24], num_elem, new_num_elem;
  THREE_PART  *three_part, *threePartPtr;
  FILE        *file_ptr;

  if( (file_ptr = fopen(model->veffMatrElem,"r")) == NULL) {
    printf("\n\nError in function read_id_three_part_matr_elem():");
    printf("\nWrong file = %s for input of eff.interaction.\n",
                                              model->veffMatrElem);
    exit(1);
  }

  // read total number of matrix elements

  if((fscanf(file_ptr,"%d%c",&num_elem,&ch)!=2)|(ch!='\n')) {
    printf("\n\n Error in function read_three_part_matr_elem():");
    printf("\nFirst number in file %s",model->veffMatrElem);
    printf(" - number of matrix element - is not correct\n");
    exit(1);
  } 
  // memory to store three-particle matrix elements

  three_part = MALLOC(num_elem, THREE_PART, func,"three_part[]");

  new_num_elem = 0;       // initialization
  threePartPtr = three_part; 

  for(loop = 0; loop < num_elem; loop++) {

     // read left (na,la,ja,ma), (nb,lb,jb.mb), (nc,lc,jc,mc)
     // and right (na,la,ja,ma), (nb,lb,jb.mb), (nc,lc,jc,mc)

    if(read_int_number(file_ptr, 24, int_data) == FALSE) {
      printf("\n\nError in function read_three_part_matr_elem();");
      printf("\nError in (n,l,j) - values");
      printf(" for matrix element number = %d\n",loop);
      exit(1);
    }

    // Read effective matrix elements

    if(read_float_number(file_ptr, 1, &threePartPtr->value) == FALSE) {
      printf("\n\nError in function read_three_part_matr_elem():");
      printf("\nWrong value for three-particle matrix element no %d",loop);
      printf("\n read from file %s\n",model->veffMatrElem);
      exit(1);
    }
   /*
   **           TEST case
   ** The effective three-particle interaction is is obtained 
   ** through a transformation of a  two-particle veff.
   ** This requires a transformation factor
   */

    if(model->typeVeff == 4) {
      threePartPtr->value *= (1.0/((double)(model->part - 2)));
    }
    // test conservation of total M-value and parity

    if(      ((int_data[3] + int_data[7] + int_data[11])
          != (int_data[15] + int_data[19] + int_data[23]))
       ||(   ((int_data[1] + int_data[5] + int_data[9]) % 2) 
          != ((int_data[13] + int_data[17] + int_data[21])% 2))) {
      printf("\n\nError in function read_id_three_part_matr_elem():");
      printf("\nM-value and parity is wrong for element no %d\n", loop);
      exit(1);
    }
       /*
       ** convert (n,l.j, m) input data to orbit numbers
       ** and store the result in struct three_part
       */
      
    id_matr_elem_identity(sp_bas, int_data, threePartPtr);

       /*
       **  Only matrix elements relevant for the present
       **  set of single-particle orbits are saved
       */

    if(  (threePartPtr->la==-1) || (threePartPtr->lb==-1)
       ||(threePartPtr->lc==-1) || (threePartPtr->ra==-1)
       ||(threePartPtr->rb==-1) || (threePartPtr->rc==-1))
      continue;
                  
     /*
     ** reorder if necessary the orbits such that 
     **       threePartPtr.la <= threePartPtr.lb and  
     **       threePartPtr.ra <= threePartPtr.rb
     ** Add a Pauli phase factor to matrix element.
     ** Calculate and store the config threePartPtr.left
     ** and threePartPtr.right and possibly interchange 
     ** such that 
     **    threePartPtr.left >=  threePartPtr.right
     */
 
    id_orbit_interchange(threePartPtr);

    threePartPtr++;   // current matrix elements accepted
    new_num_elem++;
    
   } // end loop through all input three-particle matrix elements 

   fclose(file_ptr);

      /* 
      ** All matrix elements are sorted into increasing groups of 
      ** same three_part->right and for each group into subgroups
      ** after increasing three_part->left.
      */

   sort_matrix_elements(new_num_elem, three_part);

   // remove possibly equal matrix elements

   new_num_elem = id_compress_matr_elem(new_num_elem, three_part);

  if(new_num_elem < num_elem) {  // reduce size of three_part[]
     three_part = REALLOC(three_part, new_num_elem, THREE_PART, func,
                                      "realloc j_pot.three_part[]");
   }

   *three_part_ref = three_part; // return storage address to matrix elements
   
   return new_num_elem;

} // End: function  read_id_three_part_matr_elem()

     /*
     ** The function 
     **         id_matr_elem_identity()
     ** convert (n,l.j) input data to orbit numbers
     ** and store the result in struct three_part
     */

static void id_matr_elem_identity(SP_BAS *sp_bas, int *int_data,
                                              THREE_PART *three_part)
{
  int        k;
  MBAS       *basPtr;

  three_part->la = -1;     // initialization
  three_part->lb = -1;
  three_part->lc = -1;
  three_part->ra = -1;
  three_part->rb = -1;
  three_part->rc = -1;

  basPtr       = sp_bas->mbas;

  for(k = 0; k < sp_bas->numm_orb; k++, basPtr++)     {
    if(   (basPtr->osc == int_data[0]) && (basPtr->l == int_data[1]) 
       && (basPtr->j   == int_data[2]) && (basPtr->m == int_data[3]))
      three_part->la = k;
    if(   (basPtr->osc == int_data[4]) && (basPtr->l == int_data[5]) 
       && (basPtr->j   == int_data[6]) && (basPtr->m == int_data[7]))
      three_part->lb = k;
    if(   (basPtr->osc == int_data[8]) && (basPtr->l == int_data[9]) 
       && (basPtr->j   == int_data[10]) && (basPtr->m == int_data[11]))
      three_part->lc = k;
    if(   (basPtr->osc == int_data[12]) && (basPtr->l == int_data[13]) 
       && (basPtr->j   == int_data[14])&& (basPtr->m == int_data[15]))
      three_part->ra = k;
    if(   (basPtr->osc == int_data[16]) && (basPtr->l == int_data[17]) 
       && (basPtr->j   == int_data[18]) && (basPtr->m == int_data[19]))
      three_part->rb = k;
    if(   (basPtr->osc == int_data[20]) && (basPtr->l == int_data[21]) 
       && (basPtr->j   == int_data[22]) && (basPtr->m == int_data[23]))
      three_part->rc = k;
  }

            /* "Whitehead" phase is included */
 
  three_part->value *= sp_bas->mbas[three_part->la].phase
                      *sp_bas->mbas[three_part->lb].phase
                      *sp_bas->mbas[three_part->lc].phase
                      *sp_bas->mbas[three_part->ra].phase
                      *sp_bas->mbas[three_part->rb].phase
                      *sp_bas->mbas[three_part->rc].phase;
                
} // End function  id_matr_elem_identity()

     /*
     ** The function 
     **         id_orbit_interchange()
     ** reorder if necessary the orbits such that    
     **  three_part.la < three_part.lb < three_part.lc 
     ** and  
     **  three_part.ra < three_part.rb < three_part.rc 
     ** If any change multiply matrix element with 
     ** a Pauli phase factor 
     ** Calculate and store the config three_part.left
     ** and three_part.right and possibly interchange 
     ** such that 
     **    three_part.left >=  three_part.right
     */

static void id_orbit_interchange(THREE_PART *three_part)
{
  int        temp, phase;
  ULL        number;

  phase = +1; // initialization

  // interchange left config

  number = 0;
  if(three_part->la > three_part->lb) number += ULL_ONE << 2;
  if(three_part->la > three_part->lc) number += ULL_ONE << 1;
  if(three_part->lb > three_part->lc) number += ULL_ONE;

  switch(number) {
     case 0 :  break;
     case 1 :  temp = three_part->lb;
               three_part->lb =  three_part->lc;
               three_part->lc = temp;
               phase *= -1;
               break;  
     case 3:   temp = three_part->lc;
               three_part->lc =  three_part->lb;
               three_part->lb =  three_part->la;
               three_part->la = temp;
               break;
     case 4:   temp = three_part->la;
               three_part->la =  three_part->lb;
               three_part->lb = temp;
               phase *= -1;
               break;
     case 6:   temp = three_part->la;
               three_part->la =  three_part->lb;
               three_part->lb =  three_part->lc;
               three_part->lc = temp;
               break;
     case 7:   temp = three_part->la;
               three_part->la =  three_part->lc;
               three_part->lc = temp;
               phase *= -1;
               break;
    default:   printf("\n\n Error in function id_orbit_interchange():");
               printf("Not correct interchange number = %d\n",(int)number);
	       exit(1);
  }  // end switch() left config


  // interchange right config

  number = 0;               // initialization

  if(three_part->ra > three_part->rb) number += ULL_ONE << 2;
  if(three_part->ra > three_part->rc) number += ULL_ONE << 1;
  if(three_part->rb > three_part->rc) number += ULL_ONE;

  switch(number) {
     case 0 :  break;
     case 1 :  temp = three_part->rb;
               three_part->rb =  three_part->rc;
               three_part->rc = temp;
               phase *= -1;
               break;  
     case 3:   temp = three_part->rc;
               three_part->rc =  three_part->rb;
               three_part->rb =  three_part->ra;
               three_part->ra = temp;
               break;
     case 4:   temp = three_part->ra;
               three_part->ra =  three_part->rb;
               three_part->rb = temp;
               phase *= -1;
               break;
     case 6:    temp = three_part->ra;
               three_part->ra =  three_part->rb;
               three_part->rb =  three_part->rc;
               three_part->rc = temp;
               break;
     case 7:   temp = three_part->ra;
               three_part->ra =  three_part->rc;
               three_part->rc = temp;
               phase *= -1;
               break;
    default:   printf("\n\n Error in function id_orbit_interchange():");
	       printf("Not correct interchange number = %d\n",(int)number);
	       exit(1);
  }  // end switch()

   three_part->value *= (double)phase;

       /* if necessary, interchange right and
       ** left orbit identity such that
       ** three_part->left > three_part->right
       */

   three_part->left  =  (ULL_ONE << three_part->la) + (ULL_ONE << three_part->lb) 
                      + (ULL_ONE << three_part->lc); 
   three_part->right =  (ULL_ONE << three_part->ra) + (ULL_ONE << three_part->rb) 
                      + (ULL_ONE << three_part->rc); 

  if(three_part->left < three_part->right) {
    temp              = three_part->la;
    three_part->la    = three_part->ra;
    three_part->ra    = temp;
    temp              = three_part->lb;
    three_part->lb    = three_part->rb;
    three_part->rb    = temp;
    temp              = three_part->lc;
    three_part->lc    = three_part->rc;
    three_part->rc    = temp;
    number            = three_part->left;
    three_part->left  = three_part->right,
    three_part->right = number; 
  }
  // end interchange |right> < |left>


} // End: function  id_orbit_interchange()

     /* 
     ** The function
     **        sort_matrix_elements()             
     ** sorts all matrix elements  in struct three_part()
     ** after increasing three_part()->right. 
     ** Next, for each group of three_part()->right sort the matrix 
     ** elements after increasing three_part()->left. 
     */

static void sort_matrix_elements(int num_elem, THREE_PART *three_part)
{
   int        loop, num;
   ULL        test;
   THREE_PART *ptr_start, *ptr;

  // sort matrix elements into groups after increasing three_part->right

   qsort(three_part, (size_t)num_elem, (size_t)sizeof(THREE_PART),
                      (int(*)(const void *, const void *)) right_config);
   /*
   ** sort matrix elements with same three_part->right
   ** into increasing three_part->left index
   */

   ptr_start = three_part;    // initialization;
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
       qsort(ptr_start, (size_t)num, (size_t)sizeof(THREE_PART),
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
     qsort(ptr_start, (size_t)num, (size_t)sizeof(THREE_PART),
	    (int(*)(const void *, const void *)) left_config);
   }

} /* End: function sort_matrix_elements() */

     /*
     ** The function                         
     **        int right_config()                  
     ** is a utility function for the library function qsort() in order to sort the
     ** three-body matrix elements into groups of increasing value of three_part->right.
     */

static int right_config(const THREE_PART *one, const THREE_PART *two)
{
  if(one->right > two->right)       return +1;
  else  if(one->right < two->right) return -1;
  else                              return  0;
} // End: function right_config()

     /*
     ** The function                         
     **        int left_config()                  
     ** is a utility function for the library function qsort() in order to sort the
     ** three-body matrix elements for the same three_part->right 
     ** into sunbgroups of increasing value of three_part->left.
     */

static int left_config(const THREE_PART *one, const THREE_PART *two)
{
  if(one->left > two->left)       return +1;
  else  if(one->left < two->left) return -1;
  else                            return  0;     
} // End: function left_config()


      /*
      ** The function 
      **      id_compress_matr_elem()
      ** takes num list of sorted three-particle matrix 
      ** element and removes possibly equal elements 
      ** Note: The function must have num > 1;
      ** The final number of elements is returned.
      */

static int id_compress_matr_elem(int num, THREE_PART *three_part)
{
  int   k, m, limit;

  limit = num - 1;  // initialization
  k     = 0;

  while(k < limit) {
    if(  (three_part[k].left  == three_part[k+1].left)
       &&(three_part[k].right == three_part[k+1].right))  {
      if(fabs(three_part[k].value - three_part[k+1].value) 
	                                  > MATRIX_LIMIT) {
	printf("\n\nError in function  id_compress_matr_elem():");
	printf("\nInput data file contains two equal matrix elements");
        printf(" with different values:");
	printf("\n<%d,%d,%d|V|%d,%d,%d> = %10.6lf and %10.6lf\n",
               three_part[k].la,three_part[k].lb,three_part[k].lc,
               three_part[k].ra,three_part[k].rb,three_part[k].rc,
               three_part[k].value,three_part[k+1].value);
	exit(1);
      }
      for(m = k+1; m < limit; m++) {
	three_part[m].la    = three_part[m+1].la;
	three_part[m].lb    = three_part[m+1].lb;
	three_part[m].ra    = three_part[m+1].ra;
	three_part[m].rb    = three_part[m+1].rb;
	three_part[m].left  = three_part[m+1].left;
	three_part[m].right = three_part[m+1].right;
	three_part[m].value = three_part[m+1].value;
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
     **              id_m_pot_diag(...)                    
     ** transfers all m-scheme diagonal three-particle matrix
     ** elements to the struct MATR_OP op_veff->diag
     */

static void id_m_pot_diag(SP_BAS *sp_bas, MATR_OP *op_veff, 
                      int tot_num_elem, THREE_PART *three_part)
{
  int         k, l, m, limit, step;
  ULL         right;
  double      *id_diag;
  THREE_PART  *ptr;

  id_diag    = op_veff->id_diag; // initialization
  limit      = sp_bas->numm_orb - 2;
  for(k = 0; k < limit; k++) {
    for(l = k+1; l <= limit; l++) {
      for(m = l+1; m < sp_bas->numm_orb; m++) {
	right = (ULL_ONE << k) + (ULL_ONE << l) + (ULL_ONE << m);  
	ptr  = three_part;
        step = 0;
	while( !(right == ptr->right)) {
	  ptr++;
	  if((++step >= tot_num_elem))  {
	    printf("\n\n Error in function id_m_pot_diag():");
	    printf("The m-scheme diagonal three-particle matrix ");
	    printf("element <%d,%d,%d|V|%d,%d,%d> is missing\n",k,l,m,k,l,m);
	    exit(1);
	  }
	}
	if(right != ptr->left) {
	  printf("\n\n Error in function id_m_pot_diag():");
	  printf("\nNot a diagonal matrix elements:");
	  printf("\n <%d,%d,%d|V|%d,%d,%d> \n",ptr->la, ptr->lb, ptr->lc,k,l,m);
	}

	*(id_diag++) = ptr->value; 

      } // end m index
    } // end l index
  } // end k index

} /* End: function id_m_pot_diag() */

     /*
     ** The function                                          
     **             id_m_pot_nondiag(,...)                    
     ** transfer all nondiagonal m-scheme three-particle matrix elements
     ** <h < i < j | VEFF | k < l < m> with (h,i,j) > (k,l,m) to 
     ** op_veff->id_nondiag_bas[], one group of elements for each (k,l.m).
     ** op_veff->id_nondiag_table[((2*num_m - k - 3) * k)/2 + (l - 1)] 
     ** element points to the first element of the (k,l,m) group.
     */

static int id_m_pot_nondiag(SP_BAS *sp_bas, MATR_OP *op_veff, 
                                                     THREE_PART *three_part)
{
  int        k, l, m, k_limit, l_limit, m_limit, count, number;
  ULL        right;
  ID_INT     *nondiag, **nondiag_table = NULL;
  THREE_PART *ptr;

  k_limit = sp_bas->numm_orb - 2;  // initialization
  l_limit = sp_bas->numm_orb - 1;
  m_limit = sp_bas->numm_orb;


  nondiag                 = op_veff->id_nondiag_bas;
  nondiag_table           = op_veff->id_nondiag_table;
  op_veff->maxNondiagElem = 0;

  // correction two-part-veff --> three-part-veff 

  number = 0;

  for(k = 0; k < k_limit; k++) {
    for(l = k+1; l < l_limit; l++)  {
      for(m = l+1; m < m_limit; m++, nondiag_table++) {

	*nondiag_table = nondiag; // store pointer

	right = (ULL_ONE << k) + (ULL_ONE << l) +  (ULL_ONE << m);

	ptr = three_part;
	while(ptr->right != right) ptr++;
	ptr++;
	count = 0;       // initialization
	while(ptr->right == right) {
	  nondiag->val = (ptr->value);
	  if(fabs(nondiag->val) > MATRIX_LIMIT) {
	    nondiag->one    =  (ULL_ONE<<ptr->lc) + (ULL_ONE<<ptr->lb)
		             + (ULL_ONE<<ptr->la);
	    nondiag->tr_one =  (ULL_ONE<<(l_limit - ptr->lc)) 
		                 + (ULL_ONE<<(l_limit - ptr->lb))
		                 + (ULL_ONE<<(l_limit - ptr->la));
	    (nondiag)->two  =   ((ULL_ONE<<ptr->la) - ULL_ONE)
		                      +( (ULL_ONE<<ptr->lc) 
                                        -(ULL_ONE<<(ptr->lb + 1)));
	    count++;
	    number++;  
	    nondiag++;
	  } 
	  ptr++; // new nondiag |left>
	} // end  nondiag |left>

	if(!count)  {// no elements for current (k,l,m)
	  *nondiag_table = NULL_PTR;
	}
	else   { // count elements for current (k,l,m)
	  op_veff->maxNondiagElem = MAX(op_veff->maxNondiagElem, count);
	  number++; 
	  (nondiag++)->one = ULL_ZERO;// terminate each (h,i,j)-set with a ZERO
	}
      } // end index m
    } // end index l
  } // end index k

  return number; // return total number of elements

} // End: function id_m_pot_nondiag()
