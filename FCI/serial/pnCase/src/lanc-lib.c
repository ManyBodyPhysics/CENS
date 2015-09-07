
/*******************  The module lanc-lib.c  ******************/

#include "shell.h"

         /* a macro used in function pythag() */

static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

 		    /**** function declarations ****/

void  *mem_alloc(size_t size, char *func, char *var);
       /*
       ** reserves dynamic memory in  heap using malloc(). No initialization
       ** of the elements.   
       **         size  - number of bytes reserved       
       **         func  - name of calling function
       **         var   - name of variable  
       ** Returns a void pointer to the reserved memory location
       ** or if(size == 0) returns NULL pointer. 
       **/

void  *mem_calloc(size_t num_elem, size_t size, char *func, char *var);
       /*
       ** reserves dynamic memory in  heap using calloc(). All elements are
       ** initialized to ZERO
       **         num_elem   - number of elements
       **         size       - number of bytes for each element       
       **         func       - name of calling function
       **         var        - name of variable  
       ** Returns a void pointer to the reserved memory location
       ** or if(num_elem = 0) returns a NULL pointer.
       */

void  *mem_realloc(void *ptr, size_t num_elem, size_t size, char *func, char *var);
      /*
      ** reallocates previously reserved dynamic memory in heap using realloc().
      ** In case of increased memory no initialization of the new elements.   
      **         *ptr  - pointer to the previosly reserved memory 
      **         size  - number of bytes reserved       
      **         func  - name of calling function
      **         var   - name of variable  
      ** Returns a void pointer to the new reserved memory location.
      **/

void add_a_vector(int, double, double *,double *);
       /*
       ** multiply a vector ptr_1[] by a constant, add the result to 
       ** a vector ptr_2[] and store the result in vector ptr_2[].   
       **        ptr_2[] <--- ptr_2[] + factor * ptr_1[]              
       */

double norm_of_a_vector(int, double *);
       /*
       ** calculates the norm of a vector stored in in ptr[] 
       ** with "dim" number of elements. 
       ** The function returns 
       **     result = sqrt( Sum(x_i * x_i))                
       */

void  scale_a_vector( int, double, double *);
       /*
       ** multiply each component of a vector by  factor.
       **        ptr[] <- factor * ptr[]                 
       */

void normalize_a_vector(int, double *);
       /*
       ** normalizes a vector stored in vec[]. The old vector is destroyed and the 
       ** function returns the normalized  vector stored in the same place. 
       */  

double scalar_product_of_vectors(int, double *, double *);
       /*
       ** calculates and returns the scalar product of 
       ** two vectors stored in vec_1[] and vec_2[].
       */

void eigenvalues(int, double *, double *);
       /*
       ** diagonalyzes a general real symmetrix matrix of dimension dim 
       ** stored in h_matrix[] in the format specified in function 
       ** pn_eigenvalue_lanczo_process() in file lanc-pn-process.c
       ** The function returns all eigenvalues in increasing order
       ** in list_eigenvalue[].
       */

void eigenvalue_tred2(double **a, int, double *, double *);
       /*
       ** (C) Copr. 1986-92 Numerical Recipes Software )%. 
       */

void eigenvalue_tqli(double *, double *, int);
       /*
       ** (C) Copr. 1986-92 Numerical Recipes Software )%. 
       */

double pythag(double, double);
       /*
       ** (C) Copr. 1986-92 Numerical Recipes Software )%. 
       */

int eigenvalue_comp(const double *, const double *);
       /*
       ** uses by the library function qsort() to put
       ** eigenvalues in an increasing sequence.
       */

void eigenvectors(int, double *, EIGEN *);
       /*
       ** diagonalyzes a general real symmetrix matrix of dimension dim 
       ** stored in h_matrix[] in the format specified in function 
       ** pn_eigenvalue_lanczo_process() in file lanc-pn-process.c
       ** The function returns all eigenvalues and corresponding
       ** eigenvectors  in increasing order in structure 
       ** EIGEN list_eigen[]
       */

void eigenvector_tred2(double **, int, double *, double *);
       /*
       ** (C) Copr. 1986-92 Numerical Recipes Software )%. 
       */

void eigenvector_tqli(double *, double *, int, double **);
       /*
       ** (C) Copr. 1986-92 Numerical Recipes Software )%. 
       */

int eigenvector_comp(const EIGEN *, const EIGEN *);
       /*
       ** uses by the library function qsort() to put eigenvalues and 
       ** corresponding eigenvectors in an increasing sequence.
       */

double clebsch_gordan(int j1, int j2, int j3, int m1, int m2);
       /*
       ** returns the value of the Clebsch-Gordan coefficient                                  
       ** <j1/2, m1/2, j2/2, m2/2 | j3/2, (m1 + m2)/2> 
       */

double fac_ratio(int m, int n);
       /*
       ** calculates and returns the ratio (n! / m!).           
       */

void read_Lanczos_vector_from_open_file(char *file_name, FILE *file_ptr, int n,
					              int dim, double *vec_ptr);
       /*
       ** reads from open file seqiencially  the Lanczos vectors in the format
       **          int n             -   as local vector number on current file and
       **          int dim           -   as number of components
       **          double vec_ptr[]  -   dim amplitudes 
       */

void run_time(int number, TID *run);
     /* 
     ** calculates start and stop time and returns the difference.                
     ** Input data:                            
     **    int number      - = 1 for start     
     **                        2 for stop      
     **    TID run  -  returns the time difference. 
     */
         
               /**** End: function declarations ****/

               /**** The function defin|<itions  ****/ 

       /*
       ** The function 
       **        void  *mem_alloc()
       ** reserves dynamic memory in  heap using malloc(). No initialization
       ** of the elements.   
       **         size  - number of bytes reserved       
       **         func  - name of calling function
       **         var   - name of variable  
       ** Returns a void pointer to the reserved memory location
       ** or if(size == 0) returns NULL pointer. 
       **/

void  *mem_alloc(size_t size, char *func, char *var)
{
   void  *new_memory;

    /* allow a call to library function malloc() during memory allocation */

#undef     malloc   

   if(size == 0) new_memory = NULL_PTR;
   else {
      if((new_memory = (void *) malloc(size)) == NULL)  {
         printf("\n\nError from function mem_alloc():");
         printf("\nmemory allocation from function %s: %s ", func, var); 
         printf("\nMemory request: %d bytes \n", size);
         exit(1);
      }
   }

#define malloc             "don't call malloc directly anymore!"

   return new_memory;

} /* End: function  mem_alloc() */

       /*
       ** The function 
       **        void  *mem_calloc()
       ** reserves dynamic memory in  heap using calloc(). All elements are
       ** initialized to ZERO
       **         num_elem   - number of elements
       **         size       - number of bytes for each element       
       **         func       - name of calling function
       **         var        - name of variable  
       ** Returns a void pointer to the reserved memory location
       ** or if(num_elem = 0) returns a NULL pointer.
       */

void  *mem_calloc(size_t num_elem, size_t size, char *func, char *var)
{
   void  *new_memory;

     /* allow a call to library function calloc() during memory allocation */

#undef     calloc   

   if(num_elem == 0) new_memory = NULL_PTR;
   else {
      if((new_memory = (void *) calloc(num_elem, size)) == NULL)  {
         printf("\n\nError from function mem_calloc():");
         printf("\nmemory allocation from function %s: %s ", func, var);
         printf("\nMemory request: num_elem = %d size = %d bytes \n", num_elem, size);
         exit(1);
      }
   }

#define calloc             "don't call calloc directly anymore!"

   return new_memory;

} /* End: function  mem_calloc() */

      /*
      ** The function 
      **        void  *mem_realloc()
      ** reallocates previously reserved dynamic memory in heap using realloc().
      ** In case of increased memory no initialization of the new elements.   
      **         *ptr  - pointer to the previosly reserved memory 
      **         size  - number of bytes reserved       
      **         func  - name of calling function
      **         var   - name of variable  
      ** Returns a void pointer to the new reserved memory location.
      **/

void  *mem_realloc(void *ptr, size_t num_elem, size_t size, char *func, char *var)
{
   void  *new_memory;

    /* allow a call to library function realloc() during memory allocation */

#undef     realloc   

   if((new_memory = (void *) realloc(ptr, num_elem * size)) == NULL)  {
      printf("\n\nError from function mem_realloc():");
      printf("\nmemory allocation from function %s: %s ", func, var); 
      printf("\nMemory request: %d bytes \n", size);
      exit(1);
   }

#define realloc             "don't call realloc directly anymore!"

   return new_memory;

} /* End: function  mem_realloc() */

       /*
       ** The function                                               
       **                    add_a_vector()                          
       ** multiply a vector ptr_1[] by a constant, add the result to 
       ** a vector ptr_2[] and store the result in vector ptr_2[].   
       **        ptr_2[] <--- ptr_2[] + factor * ptr_1[]              
       */

void add_a_vector(int dim, double factor, double *ptr_1,double *ptr_2)
{
  do  {
    *(ptr_2++) += factor * (*(ptr_1++));
  } while (--dim);

} /* End: function add_a_vector() */

       /*
       ** The function                                       
       **               norm_of_a_vector()                   
       ** calculates the norm of a vector stored in in ptr[] 
       ** with "dim" number of elements. 
       ** The function returns 
       **     result = sqrt( Sum(x_i * x_i))                
       */

double norm_of_a_vector(int dim, double *ptr)
{
  double
           result;

  result = D_ZERO;
  do  {
    result += (*ptr) * (*ptr);
    ptr++;
  } while (--dim);

  return sqrt(result);

} /* End: function norm_of_a_vector() */

       /*
       ** The function                                   
       **             scale_a_vector()                   
       ** multiply each component of a vector by  factor.
       **        ptr[] <- factor * ptr[]                 
       */

void  scale_a_vector( int dim, double factor, double *ptr)
{
  do  {
    *(ptr++) *= factor;
  } while (--dim);
} /* End: function scale_a_vector */

       /*
       ** The function
       **             normalize a vector()
       ** normalizes a vector stored in vec[]. The old vector is destroyed and the 
       ** function returns the normalized  vector stored in the same place. 
       */  

void normalize_a_vector(int dim, double *vec)
{
   double
             sum;

   sum = norm_of_a_vector(dim, vec);

   if(sum < ZERO_LIMIT)   {
      printf("\n\n Error in function normalize a vector()");
      printf("\nNorm of the vector is too small - norm = %f\n",sum);
      exit(1);
   }
   scale_a_vector(dim, 1.0 / sum, vec);

} /* End: function normalize_a_vector () */ 

       /*
       ** The function                                       
       **        scalar_product_of_vectors()                  
       ** calculates and returns the scalar product of 
       ** two vectors stored in vec_1[] and vec_2[].
       */

double scalar_product_of_vectors(int dim, double *vec_1, double *vec_2)
{
  double  
           result;

  result = D_ZERO;             /* initialization */
  do  {
    result += (*(vec_1++)) * (*(vec_2++));
  } while (--dim);

  return result;

} /* End: function scalar_product_of_vectors() */

      /* 
      ** The function                                    
      **      eigenvalues()                          
      ** diagonalyzes a general real symmetrix matrix of dimension dim 
      ** stored in h_matrix[] in the format specified in function 
      ** pn_eigenvalue_lanczo_process() in file lanc-pn-process.c
      ** The function returns all eigenvalues in increasing order
      ** in list_eigenvalue[].
      */

void eigenvalues(int dim, double *h_matrix, double *list_eigenvalue)
{
  char          
                    *func = {" eigenvalues(): "};
  register int
                    loop, left, right;
  double
                    *d, *e, **a;


    /* temporary memory allocation */

  d = MALLOC(dim + 1, double, func, "d[]");
  e = MALLOC(dim + 1, double, func, "e[]");
  a = MALLOC(dim + 1, double *, func,"a[][]");
  for(loop = 0; loop < dim + 1; loop++) {
    a[loop] = MALLOC(dim + 1, double, func,"a[]");
  }

   /* transfer the Lanczo energy matrix to a[][] */

  for (loop = 0, right = 1; right <= dim; right++)  {
    for (left = 1; left <= right; left++)   {
       a[left][right] = h_matrix[loop];
       a[right][left] = h_matrix[loop++];
    }
  }
    /* diagonalize the matrix a[][] */

  eigenvalue_tred2(a, dim, d, e);       /* transfer to tri-diag. form */
  
  for( loop = dim; loop >= 0; loop--) free(a[loop]);
  free(a);

  eigenvalue_tqli(d, e, dim);   /* diagonalize */

    /* return eigenvalues */

  for(left = 1; left <= dim; left++)   list_eigenvalue[left - 1] = d[left];

    /*
     * Sort eigenvalues and the corresponding eigenvectors 
     * after increasing values of the eigenvalues.         
     */

  qsort((void *) list_eigenvalue, (size_t) dim , (size_t) sizeof(double),
      (int(*)(const void *, const void *))eigenvalue_comp);

    /* release temporary memory */

  free(e);
  free(d);

} /* End: function eigenvalues() */

/* (C) Copr. 1986-92 Numerical Recipes Software )%. */

void eigenvalue_tred2(double **a, int n, double d[], double e[])
{
  register int
                  l,k,j,i;
  double
                  scale,hh,h,g,f;

  for (i=n;i>=2;i--) {
    l=i-1;
    h=scale=0.0;
    if (l > 1) {
      for (k=1;k<=l;k++)
        scale += fabs(a[i][k]);
      if (scale == 0.0)
        e[i]=a[i][l];
      else {
        for (k=1;k<=l;k++) {
          a[i][k] /= scale;
          h += a[i][k]*a[i][k];
        }
        f=a[i][l];
        g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
        e[i]=scale*g;
        h -= f*g;
        a[i][l]=f-g;
        f=0.0;
        for (j=1;j<=l;j++) {
          a[j][i]=a[i][j]/h;
          g=0.0;
          for (k=1;k<=j;k++)
            g += a[j][k]*a[i][k];
          for (k=j+1;k<=l;k++)
            g += a[k][j]*a[i][k];
          e[j]=g/h;
          f += e[j]*a[i][j];
        }
        hh=f/(h+h);
        for (j=1;j<=l;j++) {
          f=a[i][j];
          e[j]=g=e[j]-hh*f;
          for (k=1;k<=j;k++)
            a[j][k] -= (f*e[k]+g*a[i][k]);
        }
      }
    } else
      e[i]=a[i][l];
    d[i]=h;
  }
  d[1]=0.0;
  e[1]=0.0;

       /* Contents of this loop can be omitted if eigenvectors not
          wanted except for statement d[i]=a[i][i]; */

  for (i=1;i<=n;i++) {
    d[i]=a[i][i];
  }
} /* End: function eigenvalue_tred2() */

void eigenvalue_tqli(double d[], double e[], int n)
{
  register int
                 m,l,iter,i;
  double
                 s,r,p,g,f,dd,c,b;

  for (i=2;i<=n;i++) e[i-1]=e[i];
  e[n]=0.0;
  for (l=1;l<=n;l++) {
    iter=0;
    do {
      for (m=l;m<=n-1;m++) {
        dd=fabs(d[m])+fabs(d[m+1]);
        if ((double)(fabs(e[m])+dd) == dd) break;
      }
      if (m != l) {
        if (iter++ == 30) {
          printf("\n\nToo many iterations in tqli.\n");
          exit(1);
        }
        g=(d[l+1]-d[l])/(2.0*e[l]);
        r=pythag(g,1.0);
        g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
        s=c=1.0;
        p=0.0;
        for (i=m-1;i>=l;i--) {
          f=s*e[i];
          b=c*e[i];
          e[i+1]=(r=pythag(f,g));
          if (r == 0.0) {
            d[i+1] -= p;
            e[m]=0.0;
            break;
          }
          s=f/r;
          c=g/r;
          g=d[i+1]-p;
          r=(d[i]-g)*s+2.0*c*b;
          d[i+1]=g+(p=s*r);
          g=c*r-b;
        }
        if (r == 0.0 && i >= l) continue;
        d[l] -= p;
        e[l]=g;
        e[m]=0.0;
      }
    } while (m != l);
  }
} /* End: function eigenvalue_tqli() */

double pythag(double a, double b)
{
  double
           absa,absb;

  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
} /* End: function pythag() */

    /*
     * The function                           
     *            eigenvalue_comp()                
     * uses by the library function qsort() to put
     * eigenvalues in an increasing sequence.
     */

int eigenvalue_comp(const double *one, const double *two)
{
  if(*one < *two)       return -1;
  else if(*one > *two)  return +1;
  else                  return  0;
} /* End: function eigenvalue_comp() */

      /* 
      ** The function                                    
      **      eigenvectors()                          
      ** diagonalyzes a general real symmetrix matrix of dimension dim 
      ** stored in h_matrix[] in the format specified in function 
      ** pn_eigenvalue_lanczo_process() in file lanc-pn-process.c
      ** The function returns all eigenvalues and corresponding
      ** eigenvectors  in increasing order in structure 
      ** EIGEN list_eigen[]
      */

void eigenvectors(int dim, double *h_matrix, EIGEN *list_eigen)
{
  char
                    *func = {"eigenvectors(): "};
  register int
                    loop, left, right;
  double
                   *d, *e, **a;


    /* temporary memory allocation */

  
  d = MALLOC(dim + 1, double, func, "d[]");
  e = MALLOC(dim + 1, double, func, "e[]");
  a = MALLOC(dim + 1, double *, func,"a[][]");
  for(loop = 0; loop < dim + 1; loop++) {
    a[loop] = MALLOC(dim + 1, double, func,"a[]");
  }

   /* transfer the Lanczo energy matrix to a[][] */

  for (loop = 0, right = 1; right <= dim; right++)  {
    for (left = 1; left <= right; left++)   {
       a[left][right] = h_matrix[loop];
       a[right][left] = h_matrix[loop++];
    }
  }
    /* diagonalize the matrix a[][] */

  eigenvector_tred2(a, dim, d, e);       /* transfer to tri-diag. form */
  eigenvector_tqli(d, e, dim, a);        /* diagonalize */

                /* return eigenvalues */

  for(left = 1; left <= dim; left++) list_eigen[left - 1].value = d[left];

    /*
     * Sort eigenvalues and the corresponding eigenvectors 
     * after increasing values of the eigenvalues.         
     */

                /* return eigenvactors */

  for(left = 1; left <= dim; left++) {

                /* transfer eigenvector no left */

    for(right = 1; right <= dim; right++) {
      list_eigen[left - 1].vector[right - 1] = a[right][left];
    }
  }
  qsort((void *) list_eigen, (size_t) dim , (size_t) sizeof(EIGEN),
      (int(*)(const void *, const void *))eigenvector_comp);

    /* release temporary memory */

  for( loop = dim; loop >= 0; loop--) free(a[loop]);
  free(a);
  free(e);
  free(d);

} /* End: function eigenvectors() */



/* (C) Copr. 1986-92 Numerical Recipes Software )%. */


void eigenvector_tred2(double **a, int n, double d[], double e[])
{
  register int
                  l,k,j,i;
  double
                 scale,hh,h,g,f;

  for (i=n;i>=2;i--) {
    l=i-1;
    h=scale=0.0;
    if (l > 1) {
      for (k=1;k<=l;k++)
        scale += fabs(a[i][k]);
      if (scale == 0.0)
        e[i]=a[i][l];
      else {
        for (k=1;k<=l;k++) {
          a[i][k] /= scale;
          h += a[i][k]*a[i][k];
        }
        f=a[i][l];
        g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
        e[i]=scale*g;
        h -= f*g;
        a[i][l]=f-g;
        f=0.0;
        for (j=1;j<=l;j++) {
          a[j][i]=a[i][j]/h;
          g=0.0;
          for (k=1;k<=j;k++)
            g += a[j][k]*a[i][k];
          for (k=j+1;k<=l;k++)
            g += a[k][j]*a[i][k];
          e[j]=g/h;
          f += e[j]*a[i][j];
        }
        hh=f/(h+h);
        for (j=1;j<=l;j++) {
          f=a[i][j];
          e[j]=g=e[j]-hh*f;
          for (k=1;k<=j;k++)
            a[j][k] -= (f*e[k]+g*a[i][k]);
        }
      }
    } else
      e[i]=a[i][l];
    d[i]=h;
  }
  d[1]=0.0;
  e[1]=0.0;

        /* 
	** Contents of this loop can be omitted if eigenvectors
	** not wanted except for statement d[i]=a[i][i];
        */

  for (i=1;i<=n;i++) {
    l=i-1;
    if (d[i]) {
      for (j=1;j<=l;j++) {
        g=0.0;
        for (k=1;k<=l;k++)
          g += a[i][k]*a[k][j];
        for (k=1;k<=l;k++)
          a[k][j] -= g*a[k][i];
      }
    }
    d[i]=a[i][i];
    a[i][i]=1.0;
    for (j=1;j<=l;j++) a[j][i]=a[i][j]=0.0;
  }

} /* End: function eigenvector_tred2() */

void eigenvector_tqli(double d[], double e[], int n, double **z)
{
  register int
                 m,l,iter,i,k;
  double
                 s,r,p,g,f,dd,c,b;

  for (i=2;i<=n;i++) e[i-1]=e[i];
  e[n]=0.0;
  for (l=1;l<=n;l++) {
    iter=0;
    do {
      for (m=l;m<=n-1;m++) {
        dd=fabs(d[m])+fabs(d[m+1]);
        if ((double)(fabs(e[m])+dd) == dd) break;
      }
      if (m != l) {
        if (iter++ == 30) {
          printf("\n\nToo many iterations in tqli.\n");
          exit(1);
        }
        g=(d[l+1]-d[l])/(2.0*e[l]);
        r=pythag(g,1.0);
        g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
        s=c=1.0;
        p=0.0;
        for (i=m-1;i>=l;i--) {
          f=s*e[i];
          b=c*e[i];
          e[i+1]=(r=pythag(f,g));
          if (r == 0.0) {
            d[i+1] -= p;
            e[m]=0.0;
            break;
          }
          s=f/r;
          c=g/r;
          g=d[i+1]-p;
          r=(d[i]-g)*s+2.0*c*b;
          d[i+1]=g+(p=s*r);
          g=c*r-b;

  	     /*
	     ** Next loop can be obmitted if
	     ** eigenvectors are not wanted
	     */
 
          for (k=1;k<=n;k++) {
            f=z[k][i+1];
            z[k][i+1]=s*z[k][i]+c*f;
            z[k][i]=c*z[k][i]-s*f;
          }
        }
        if (r == 0.0 && i >= l) continue;
        d[l] -= p;
        e[l]=g;
        e[m]=0.0;
      }
    } while (m != l);
  }
} /* End: function eigenvector_tqli() */

    /*
    ** The function                           
    **            eigenvector_comp()                
    ** uses by the library function qsort() to put eigenvalues and 
    ** corresponding eigenvectors in an increasing sequence.
    */

int eigenvector_comp(const EIGEN *one, const EIGEN *two)
{
  if(one->value < two->value)       return -1;
  else if(one->value > two->value)  return +1;
  else                              return  0;
} /* End: function eigenvector_comp() */

       /*
       ** The function                                 
       **           clebsch_gordan()                   
       ** returns the value of the Clebsch-Gordan coefficient                                  
       ** <j1/2, m1/2, j2/2, m2/2 | j3/2, (m1 + m2)/2> 
       */

double clebsch_gordan(int j1, int j2, int j3, int m1, int m2)
{
         /*   fac[n] = n! / sqrt{10^n}  */
   
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
       ** Interchange the angular momenta such that the smallest j-value 
       ** is placed in position 2. A parameter is set to remember the interchange 
       **    change  = -1 - j1 is placed in position 2.  
       **            =  0 - no change                       
       **            =  1 - j3 is placed in position 2.  
       **/

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
       ** The function               
       **      double fac_ratio()    
       ** calculates and returns the ratio (n! / m!).           
       */

double fac_ratio(int m, int n)
{
  int i, high;
  double value;

  value = 1.0;      /* initialization */
  if( n == m)    return (value);          /* Return for n equal m */
  high = MAX(n,m); 
  for (i = MIN(n,m)+1; i <= high; i++) value*= i;
  if (n < m) return (1.0/value);       /* Return for n less  m  */
  return (value);                  /* Return for n greater than m */
} /* End: function fac_ratio() */

     /*
     ** The function                                           
     **            read_Lanczos_vector_from_open_file( )                    
     ** reads from open file seqiencially  the Lanczos vectors in the format
     **          int n             -   as local vector number on current file and
     **          int dim           -   as number of components
     **          double vec_ptr[]  -   dim amplitudes 
     */

void read_Lanczos_vector_from_open_file(char *file_name, FILE *file_ptr, int n,
                                                       int dim, double *vec_ptr)
{
/***
  int     k;
*****/
  int 
        vec_no, num_elem;

   if( fread((void *)&vec_no,(size_t) sizeof(int), 1, file_ptr) != 1)  { /* vector number */
      printf("\nError in function read_Lanczos_vector_from_open__file()");
      printf("\nIn reading the vector number from file %s\n",file_name);
      exit(1);
   }   /* end of if-test  */

   if(n != vec_no) {
     printf("\nError in function read_Lanczos_vector_from_open__file()");
     printf("\nread vector number = %d is not the wanted n = %d\n", vec_no, n);
     exit(1);
   }
   if( fread((void *)&num_elem,(size_t) sizeof(int),1,file_ptr)!=1) {/* read dim */
      printf("\nError in function read_Lanczos_vector_from_open__file()");
      printf("\nIn reading dimension from  file %s\n",file_name);
      exit(1);
   }   /* end of if-test  */

   if(dim  != num_elem) {
     printf("\nError in function read_Lanczos_vector_from_open__file()");
     printf("\nread vector dimension = %d is not the wanted dimension = %d\n", num_elem, dim);
     exit(1);
   }
   if(fread((void *)vec_ptr,(size_t) sizeof(double),                /* readvector */
       (size_t) dim, file_ptr) != (size_t) dim)  {
      printf("\nError in function read_Lanczos_vector_from_open__file()");
      printf("\nIn reading vector number %d from file %s\n",n, file_name);
      exit(1);
   } /* end of if-test  */
/*******************   test output ************
  printf("\n\nread  Lanczos vector no %d in file %s\n", n, file_name);
  for(k = 0; k < num_elem; k++) {
    printf("  %10.6f",vec_ptr[k]);
    if((((k + 1)/6) * 6) == (k + 1)) printf("\n");
  }
  printf("\n");
******************    end  ******************/

} /* End: function read_Lanczos_vector_from_open_file() */
    /*
    ** The function                           
    **      void run_time(..)                 
    ** calculates start and stop time and returns the difference.                
    ** Input data:                            
    **    int number      - = 1 for start     
    **                        2 for stop      
    **    TID run  -  returns the time difference. 
    */

void run_time(int number, TID *run)
{
  clock_t
                 num, stop = 0;
  static clock_t
                 start;

  if(number == 1) {
    start = clock();
    return;
  }
  else if(number == 2)  {
    stop = clock();
  }
  num       = ((stop - start) /CLOCK_UNIT);
  run->sec  = (UL) (num % 60);
  run->min  = (UL) (num / 60);
  run->hour = (UL) (run->min/ 60);
  run->min  = (UL) (run->min % 60);
  
} /* End: of function run_time(..) */
