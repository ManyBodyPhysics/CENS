     /* The program module                              
     **                      shell.h                    
     ** for proton - neutron the Lanczo program package.
     */

     /**  Standard ANSI-C include files **/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <string.h> 
#include <time.h>
#include <ctype.h>

#define   UL         unsigned int
#define   UL_ONE     1u
#define   UL_ZERO    0u
#define   D_ZERO     0.0
#define   D_ONE      1.0
#define   TRUE       1
#define   FALSE      0
#define   YES        1
#define   NO         0
#define   EVEN       1
#define   ODD        0
#define   PROTON     0
#define   NEUTRON    1
#define   PLUS       0
#define   MINUS      1
#define   VEFF       0
#define   ANG        1
#define   MAX_M_ORBITS   (8 * sizeof(int))
#define   FINISHED     2
#define   NULL_PTR     (void *) 0

#define CLOCK_UNIT    CLOCKS_PER_SEC    /* Use of special library functions */

          /* max. factorials in Clebsch-Gordan coefficients */

#define FAC_LIM                300
#define ONE_LINE                80
#define ZERO_LIMIT         1.0E-10
#define MATRIX_LIMIT       1.0E-06
#define ENERGY_LIMIT       1.0E-05           /* limit of energy changes through Lanczos iterartion */
#define ANG_LIMIT          1.0E-05           /* limit ang mom admixture in Lanczos basis vectors   */


      /* Dynamical memory allocation procedure */

#define malloc             "don't call malloc directly!"
#define MALLOC(num_elem, type, func, text)\
               (type *)mem_alloc((num_elem) * sizeof(type), func, text) 
#define calloc             "don't call calloc directly!"
#define CALLOC(num_elem, type, func, text)\
               (type *)mem_calloc((num_elem), sizeof(type), func, text) 
#define realloc             "don't call realloc directly!"
#define REALLOC(ptr, num_elem, type, func, text)\
               (type *)mem_realloc(ptr, (num_elem), sizeof(type), func, text) 

     /* Macro definitions for integer arguments only */

#define   MAX(a,b)   ( ((a) > (b)) ? (a) : (b) )
#define   MIN(a,b)   ( ((a) < (b)) ? (a) : (b) )
#define   PARITY(a)  ( (a) % 2) 
#define   PHASE(a)   (1 - 2 * (abs(a) % 2))
#define   MOD(a,b)   ((a) - ((a)/(b)) * b)
#define   SIGN(a,b) ((b)<0 ? -fabs(a) : fabs(a))


typedef   struct {             /* The spherical single-particle orbits */
  int
              osc,   /* osc. quantum number              */
	        l,   /* orbital ang. momemtum            */
	        j,   /* twice total ang. momentum        */
         min_part,   /* min particles in spherical orbit */    
         max_part;   /* max particles in spherical orbit */
  double
                e;   /* single-part. energy              */
}JBAS;
      
typedef struct {   /* Data structure to store non-diagonal SD matrix elements <SD'|OP|SD> */
  double
                            value;   /* value of <SD'|OP|SD>  */
  int
                            final;   /* number of final |SD'> */
} STORE;

typedef   struct  {   /* structure definition for execution time */
  int
              tick,
               sec,
               min,
              hour;
} TID;

typedef   struct { /* m-scheme quantum numbers */
  int  
              orb,   /* j-orbit number - j_max first     */
	      osc,   /* oscillator number                */
	        l,   /* orbital ang. momentum            */
	      par,   /* parity of orbit (+1,-1)          */
	        j,   /* twice total spin                 */
  	        m,   /* twice proj. of total spin        */
	    phase;   /* "Whitehead" phase factor         */   
  double
                e;   /* single-particle energy           */
} MBAS;

typedef  struct {/* information of j-orbits with particle limitations */
  int     
              orb,   /* j_orbit number                        */
              min,   /* lower particle limit                  */
              max;   /* upper particle limit                  */
  UL
             mask;   /* define all m-orbits for given j-orbit */
} MASK;
   
typedef   struct { /* information list of j-orbits with particle limitations */
  int
            num;   /* number of j-orbits with limitations */
  MASK
         *table;   /* list  of j-orbits with limitations  */ 
} PART_LIM;

typedef struct { /* identical two--body matrix elements <ij|OP|kl> **/
  UL
               one,  /* identify orbits i and k     */
               two;  /* two-body pauli phase factor */
    double
               val;  /* value of matrix elements    */
} ID_INT;

typedef struct  {     /* proton-neutron two--body matrix elements <i(p)j(n)|OP|k(p)l(n) */
  UL
               one[2],  /* proton([0] orbit i / neutron ([1]) orbit j            */   
               two[2];  /* proton([0])/neutron([1]) two-body pauli phase factor  */
  double
                 val;  /* value of matrix elements            */
} PN_INT;

typedef struct {       /* index table to store two--body matrix elements */
  double
                    diag_val;  /* value of diagonal matrix element   */
  ID_INT
              *start_nondiag;  /* ptr to nondiagonal matrix elements */
} ID_TABLE;

typedef struct {  /* index table to store proton-neutron two--body matrix elements */
  double
                    diag_val;  /* value of diagonal matrix element   */
  PN_INT
              *start_nondiag;  /* ptr to nondiagonal matrix elements */
} PN_TABLE;

/*******************   old version ****************
typedef struct  {  * pointers to table of two-body matrix elements *
  ID_TABLE
            *id_table[2]; * table of two-body matrix elements                 *
		          *  for protons([0]) and neutrons([1])               *
  PN_TABLE
              *pn_table;  * table of proton-neutron two-body matrix elements  *
} MATR_OP;
******************    end old version *************/

/*******************   new version ****************/

typedef struct  {  /* pointers to table of two-body matrix elements */
  double
	      *id_diag[2],  /* list of diagonal <xx|V|kl> elements,   */ 
	         *pn_diag;  /* for eack (kl) pair, (pp->[0], nn->[1]) */
  ID_INT
          **id_nondiag[2];  /* list of addresses to list of nondiagonal */
                            /* <xx|V|kl> elements for each (kl) pair,   */
                            /* (pp->[0], nn->[1])                       */
  PN_INT
            ***pn_nondiag;  /* list of addresses to list of addresses for */
                            /* each del_mZ and del_pZ list of nondiagonal */
                            /* <xx|V|KL> elements for each (kl) pair      */
} MATR_OP;
/******************   end new version *************/

    /*
    ** Data structure to store the eigenvalues 
    ** and the corresponding eigenvectors.     
    */

typedef   struct  {
    double
              value,
              *vector;
} EIGEN;

    /* The input data for the shell model calculation */

typedef    struct  {
  char
     type_calc[ONE_LINE],  /* text for type of calculation process         */
         title[ONE_LINE],  /* Input title for the calculation              */
       in_data[ONE_LINE],  /* input data                                   */ 
      out_data[ONE_LINE],  /* resulting out_data                           */
       h_final[ONE_LINE],  /* data to restart the Lanczos process          */
       veff_pp[ONE_LINE],  /* proton-proton matrix elements in J-scheme    */
       veff_nn[ONE_LINE],  /* neutron-neutron matrix elements in J-scheme  */
       veff_pn[ONE_LINE],  /* proton-neutron matrix elements in J-scheme   */
     eigen_vec[ONE_LINE],  /* final eigenvector file                       */
     start_vec[ONE_LINE];  /* start vector file                            */
   int
                 control,   /* identification of input data read                */
               part_type,   /* particle coding from particle_structure()        */
                       P,  /* parity  states                                   */
                      MJ,  /* twice proj. of total ang.mom.(2*M)               */
                 j_value,  /* = EVEN(+1) or ODD(0) for even(odd) tot_J         */  
                 part[2],  /* particle number,proton[0]/neutron[1]             */ 
                num_j[2],  /* number of spherical orbits, proton[0]/neutron[1] */
          max_iterations,  /* maximum number  of lanczo iterations             */
                  states,  /* number of eigenstates wanted                     */
                  tot_2J,  /* ang. mom. for start vectors                      */
           num_start_vec,  /* number of start vectors                          */
                 *vec_no,  /* list of start vectors                            */
             mem_nondiag,  /* max number of nondiag STORE matr elem in memory  */
            file_nondiag;  /* curr max num of nondiag STORE matr elem on files */
  JBAS
                *jbas[2];  /* ptr. to spherical orbits,proton[0]/neutron[1]    */ 
} INPUT;

     /* The shell model data of proton-neutron  many-body system */

typedef    struct  {
  char
         title[ONE_LINE],  /* Input title for the calculation                      */
     type_calc[ONE_LINE],  /* type of calculation process                          */
      out_data[ONE_LINE];  /* resulting out_data                                   */
  int
               part_type,  /* particle coding from particle_structure()            */
                       P,  /* parity  of states (+1, -1)                           */
                      MJ,  /* twice proj. of total ang.mom.(2*M)                   */
                 part[2],  /* particle number,proton[0]/neutron[1]                 */ 
                num_j[2],  /* number of spherical orbits, proton[0]/neutron[1]     */
             same_pn_orb;  /* = YES if proton/neutron orbits are equal, els NO     */
  double 
       *list_eigenvalues;  /* table for final eigenvalues from Lanczos iteration   */
  JBAS
                *jbas[2];  /* ptr to spherical orbits,proton[0]/neutron[1]         */ 
} PN_SHELL;

        /* Data information common to all proton/neutron group */

typedef     struct {
  char 
    store_veff[ONE_LINE],  /* store diag/nondiag SD veff matrix elements   */
    store_ang[ONE_LINE];  /* store diag/nondiag SD J**2 matrix elements   */
  int
                part[2], /* number of particles  for |SD> proton[0]/neutron[1]  */
               m_orb[2], /* number of m-orbits  for |SD> proton[0]/neutron[1]   */
                   parZ,  /* possible parities of proton |SD(Z)>                */
                          /* coding: -1: parZ  -1; both parZ = 0; +1: parZ = +1 */
                      P,  /* parity  states                                     */
                     MJ,  /* twice proj. of total ang.mom.(2*M)                 */
                 num_gr, /* total number of groups                              */
            tot_numSD_Z, /* total number of |SD(Z)>                             */
            tot_numSD_N, /* total number of |SD(N)>                             */
              tot_numSD, /* total number of |SD(Z),SD(N)>                       */
               maxSD[2], /* max |SD()> within a group - proton [0]/ neutron [1] */ 
                *occ[2], /* lists to store occ orb in all |SD(Z)> and SD(N)>    */
                         /* in a group                                          */ 
            mem_nondiag, /* max number of nondiag STORE matr elem in memory      */
           file_nondiag; /* curr max num of nondiag STORE matr elem on files     */
  MBAS
               *mbas[2];  /* ptr to m-scheme orbits, proton[0]/neutron[1]        */
  PART_LIM
            part_lim[2]; /* data of j-orbits with particle limitations           */
  MATR_OP
                 op_int,  /* table of m-scheme effective two-particle            */
                          /* matrix elements for pp, nn and pn                   */
                 op_ang;  /* table of m-scheme angular momentum two-particle     */
                          /* matrix elements for pp, nn and pn                   */
  char
               *mem_pos,  /* dynamic memory to store different types of data     */
              *mem_free;  /* position to first free byte reserved dynamic memory */
} GR_BASIS;                           
   
        /* Data information  for a proton/neutron group */

typedef   struct   {
  int                    /* for OP = VEFF and ANG                                */
                   m[2], /* ang.mom. projection for |SD> proton[0]/neutron[1]    */
                 par[2], /* parity P_Z(= +1,-1) for |SD> for |SD>                */
                         /* proton[0]/neutron[1]                                 */
               numSD[2], /* number of |SD>      for |SD> proton[0]/neutron[1]    */
                 max_gr, /* max group number for nondiag matrix elements         */
              start_amp; /* start position for Lanczos vector amplitudes         */
  UL
                 *SD[2]; /* ptr to list of |SD>     for |SD(Z)>|SD(N)>           */
} GROUP;

typedef      struct   {
  char 
             type_calc[ONE_LINE],  /* type of calculation process                */
              out_data[ONE_LINE],  /* resulting out_data                         */
            lanc_store[ONE_LINE],  /* file for lanczos vectors                   */
             eigen_vec[ONE_LINE],  /* final eigenvector file                     */
             start_vec[ONE_LINE],  /* start vector file                          */
               h_final[ONE_LINE];  /* storage of h_matrix[] prepared for next run*/
  int
                   num_start_vec,  /* number of start vectors                    */
                         *vec_no,  /* list of start vectors                      */
                  max_iterations,  /* maximum number  of lanczo iterations       */
                  num_iterations,  /* actual number of lanczo iterations         */
                          states,  /* number of eigenstates wanted / found       */
                          tot_2J,  /* ang. mom. for start vectors                */
                        num_totJ,  /* number of possible J_values in the model   */
                     mem_nondiag,  /* max number of nondiag STORE matr elem in memory*/
                    file_nondiag,  /* curr max num of nondiag STORE matr elem on file*/
                        run_code;  /* type of termination for the Lanczos process */
  double
	            *eigenvalues;  /* table of final eigenvalues                  */
} ITERATION;

        /* structure for initial and final vectors in a lanczos process */

typedef struct {
  double     
                 *one,
                 *two;
} LANC;

           /* function declaration */

/* file: lanc-main.c */
    int main(int , char **);

/* file: lanc-input.c */
    void input_process(INPUT *);

/* file: lanc-pn-sp-basis.c */
    void pn_single_particle_basis(PN_SHELL *, GR_BASIS *);

/* file: lanc-pn-SD-basis.c */
    int pn_slater_determinant(PN_SHELL *, GR_BASIS *, GROUP **);

/* file: lanc-pn-pot.c */
    void proton_neutron_effective_interaction(INPUT *, PN_SHELL *, GR_BASIS *); 

/* file: lanc-pn-store.c */
    void pn_store_SD_matr_elem(int, GR_BASIS *, GROUP *);

/* file: lanc-pn-process.c */
void pn_eigenvalue_Lanczos_process(GR_BASIS *, GROUP *,ITERATION *);

/* file: lanc-pn-matrix-vector.c */
void matrix_vector_multiplication(GR_BASIS *, GROUP *, int, LANC); 


/* file: lanc-lib.c */

void  *mem_alloc(size_t, char *, char *);
void  *mem_calloc(size_t, size_t, char *, char *);
void  *mem_realloc(void *, size_t, size_t, char *, char *);
void add_a_vector(int, double, double *,double *);
double norm_of_a_vector(int, double *);
void  scale_a_vector( int, double, double *);
void normalize_a_vector(int, double *);
double scalar_product_of_vectors(int, double *, double *);
void eigenvalues(int, double *, double *);
void eigenvalue_tred2(double **a, int, double *, double *);
void eigenvalue_tqli(double *, double *, int);
double pythag(double, double);
int eigenvalue_comp(const double *, const double *);
void eigenvectors(int, double *, EIGEN *);
void eigenvector_tred2(double **, int, double *, double *);
void eigenvector_tqli(double *, double *, int, double **);
int eigenvector_comp(const EIGEN *, const EIGEN *);
double clebsch_gordan(int, int, int, int, int);
double fac_ratio(int, int);
void read_Lanczos_vector_from_open_file(char *, FILE *, int n, int, double *);
void run_time(int, TID *);
