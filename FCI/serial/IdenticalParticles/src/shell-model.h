
     /* The program module                                   
     **        shell-model.h
     ** the shell model calculation.
     */

     //  Standard ANSI-C include files

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>

#include <string.h> 
#include <time.h>
#include <ctype.h>
#include <sys/time.h>

       // definition of filenames 

#define RESULT_OUTPUT            "-result.dat"
#define SD_STORE_VEFF            "-SD-veffJ"
#define SD_STORE_ANG             "-SD-J**2"
#define SD_STORE_CM             "-SD-CM"
#define STORE_MEMORY            0
#define STORE_FILE              1



#define DIAG                     "-diag-elem.dat"
#define NONDIAG                  "-nondiag-elem.dat"
#define LANCZO_VECTORS           "-lanc-vectors.dat"
#define EIGEN_VECTORS            "-eigen-vectors.dat"
#define CM_VECTORS               "-CM-vectors.dat"

   //  Definition of interactions

#define VEFF_INT        0
#define VEFF_INT3       1
#define ANG_INT         2
#define CM_INT          3

// definition of sizes for memory and file allocations

#define PERM_MEM        500 // storage of <SD'|veff|SD> in units of Mbytes
#define WIND_MEM        100 // temp storage of <SD'|veff|SD in units of Mbytes

#define LANC_FILES      10  // lanczos iteration files in units of Gbytes
#define MATR_ELEM_FILES 10  // file storage of <SD'|veff|SD> in units of Gbytes
 
#define M_BYTES         1000000 // unit trasfer: Mbytes --> bytes   
#define MAX_FILE_SIZE   1000    // in units of Mbytes

#define   D_ZERO     0.0
#define   D_ONE      1.0
#define   TRUE       1
#define   FALSE      0
#define   YES        1
#define   NO         0
#define   EVEN       1
#define   ODD        0
#define   NULL_PTR   (void *) 0
#define CLOCK_UNIT   CLOCKS_PER_SEC    // Use of special library functions   

          // max. factorials in Clebsch-Gordan coefficients   

#define FAC_LIM                300
#define ONE_LINE                80
#define ZERO_LIMIT         1.0E-10
#define MATRIX_LIMIT       1.0E-06
#define   VEFF             0
#define   ANG              1
#define   CM               2
#define ENERGY_LIMIT       1.0E-06    // limit of energy changes through Lanczos iterartion   
#define ANG_LIMIT          1.0E-07    // limit ang mom admixture in Lanczos basis vectors     

#define   UL           unsigned int
#define  DOUBLE_INT    YES             // 64 bits version (YES) 32 bits version (NO)  

#if DOUBLE_INT
    #define   ULL            unsigned long long int
    #define   ULL_ONE        1ULL
    #define   ULL_ZERO       0ULL
    #define   MAX_M_ORBITS   (8 * sizeof(ULL))
#else
    #define   ULL            unsigned int
    #define   ULL_ONE        1U
    #define   ULL_ZERO       0U
    #define   MAX_M_ORBITS   (8 * sizeof(UL))

#endif


      // Dynamical memory allocation procedure   

#define malloc             "don't call malloc directly!"
#define MALLOC(num_elem, type, func, text)\
               (type *)mem_alloc((num_elem) * sizeof(type), func, text) 
#define calloc             "don't call calloc directly!"
#define CALLOC(num_elem, type, func, text)\
               (type *)mem_calloc((num_elem), sizeof(type), func, text) 
#define realloc             "don't call realloc directly!"
#define REALLOC(ptr, num_elem, type, func, text)\
               (type *)mem_realloc(ptr, (num_elem), sizeof(type), func, text) 

     // Macro definitions for integer arguments only   

#define   MAX(a,b)   ( ((a) > (b)) ? (a) : (b) )
#define   MIN(a,b)   ( ((a) < (b)) ? (a) : (b) )
#define   PARITY(a)  ( (a) % 2) 
#define   PHASE(a)   (1 - 2 * (abs(a) % 2))
#define   MOD(a,b)   ((a) - ((a)/(b)) * b)
#define   SIGN(a,b) ((b)<0 ? -fabs(a) : fabs(a))



typedef   struct  {   // structure definition for execution time   
  unsigned long long int
                          tick,
                           sec,
                           min,
                          hour;
} TID;

typedef   struct {                        // The spherical single-particle orbits   
  int
              osc,   // osc. quantum number
                l,   // orbital ang. momemtum 
                j,   // twice total ang. momentum 
         min_part,   // min allowed particles in spherical orbit
         max_part;   // max allowed particles in spherical orbit
  double
                e;   // single-part. energy
} JBAS;
      
typedef      struct {
  int  eigenVecNo,  
           tot_2J, 
       senOverlap;
} OVERLAP;

    // The basic data for the shell model calculation   

typedef    struct  {
  char
     type_calc[ONE_LINE],  // text for type of calculation process
         title[ONE_LINE],  // Input title for the calculation
     start_vec[ONE_LINE],  // start vector file                                    
      out_data[ONE_LINE],  // out-data   
  veffMatrElem[ONE_LINE],  // Effective interaction matrix elements
  centerOfMass[ONE_LINE],  // center_of_Mass matrix elements

             *permMemPtr,  // permanent memory to store diag/nondiag matrix elements
             *windMemPtr;  // window memory to store diag/nondiag matrix elements 
  int
                typeVeff,  // type of input effective interaction
                       P,  // parity  states
                      MJ,  // twice proj. of total ang.mom.(2*M)
                  J_type,  // = EVEN(+1) or ODD(0) for even(odd) tot_J
                    part,  // particle number
                 calc_CM,  // calculate Center_of_Mass effect(yes = 1, no = 0) 
                numj_orb,  // number of input single particle j-orbits
          max_iterations,  // maximum number  of lanczo iterations
                  states,  // number of eigenstates wanted
                 memCase,  // = 0 - only memory for diag VEFF, = 1 also memory for nondiag VEFF 
            file_nondiag,  // Max. total file storage for nondiag matrix elements (Mbytes)
        sect_file_no_ang,  // number of last section file for J**2 matrix elements
         sect_file_no_CM,  // number of last section file for CM matrix elements
       sect_file_no_veff,  // number of last section file for nondiag veff matrix elements
             file_lanczo,  // Max. total file storage for lanczos iteration vectors (Mbytes)
               sen_value,  // choice of seniority calc - see input datafile
           *list_sen_val,  // list of 2J-values for seniority diag
              numOverlap;  // number of seniority overlaps to calculate
  UL
             permMemSize,  // mem (bytes) to store perm. VEFF diag/nondiag on mem      
             windMemSize;  // mem (bytes) to store all types diag/nondiag on files 
  JBAS
                   *jbas;  // list of single-particle j-orbits
} SHELL;

typedef   struct {                        // m-scheme quantum numbers   
  int  
              orb,   // j-orbit number - j_max first 
              osc,   // oscillator number
                l,   // orbital ang. momentum
              par,   // parity of orbit (+1,-1)
                j,   // twice total spin 
                m,   // twice proj. of total m-value
            phase;   // "Whitehead" phase factor
  double
                e;   // single-particle energy
} MBAS;

         /* structure to specify particle reduction
	 ** or maskng in j-orbitals
	 */

typedef    struct {
  int 
         num,   // number of reduced j-orbits
        *lim;   // upper and lower particle limits for j-orbits
  ULL
       *list;   // list of j-orbits mask structure
} MASK;

         /* complete data list for spherical-j  and
	 ** m-scheme single-particle orbits_orbits  
	 ** Possible limitations on paricle number
	 ** in j-orbits are included
	 */

typedef    struct {
  int
         numj_orb,   // number of j-orbits 
         numm_orb;   // number of m-orbits
  JBAS
            *jbas;   // list of j-orbits
  MBAS
            *mbas;   // list of m-orbits
  MASK
             mask;   // structure to specify reduction/masking
                     // of number of particles in j-orbits  
} SP_BAS;

      // data structure for the basis |SD> states 

typedef   struct {
   int
                            P, // parity  states
                           MJ, // twice proj. of total ang.mom.(2*M)
                     numm_orb, // number of m-orbits
                       J_type, // = EVEN(+1) or ODD(0) for even(odd) tot_J
                         part, // number of identical particles
                     num_totJ, // number of different J-values in the model 
                     numSD[2], // total number of |SD()>, asym/nosym = 0, sym = 1
                    tot_dimSD; // total number of basis |SD()>
   ULL
                          *SD;  // |SD()> config, stored as |asym_SD> (+ |sym_SD>)             
} SD_BAS;

typedef     struct   {   // data structure for Center_of_Mass process
    int            N,
             minPart,
             maxPart;
} OSC_STRUCT;

   typedef   struct {
           int
                       a,    /* left orbit  - particle one */
                       b,    /* left orbit  - particle two */
                       c,    /* right orbit - particle one */
                       d;    /* right orbit - particle two */
           ULL
                    left,    /* left two-particle config   */
                   right;    /* right part-particle config */
   } ORB_ID;

typedef struct { // identify orbits in two- or three-particle case 
                 // for elements in m-scheme table: 
                 //  H*|SD_init> --> value*|SD_final>
  ULL
               one, // identify occ particle orbits
            tr_one, // identify time-reversed occ particle orbits
               two; // Pauli phase   
    double
               val;  // value of matrix elements      
} ID_INT;


typedef struct  {  // pointers to table of two-body matrix elements   
  int
         num_nondiag_elem,  // num reserved nondiag elem
           maxNondiagElem;  // max num <..|op|ij>,
  double
         	 *id_diag;  // list of diagonal <ij|V|ij> elements,      
  ID_INT
          *id_nondiag_bas,  // pointer to start of nondiag ID_INT elements   
                            // <xx|V|kl> elements for each (kl) pair,     
       **id_nondiag_table;  // kl_table of pointer to nondiag matrix elements
} MATR_OP;

typedef struct {   // Data structure to write/read  <SD'|OP|SD> to/from file
  float
                   value;   // value of <SD'|OP|SD>
  ULL
                   final;   // number of final |SD'>
} STORE_F;

typedef    struct {
  char
               title[ONE_LINE],  //  Input title for the calculation
           sect_file[ONE_LINE],  // filename for the nondiag section files
              result[ONE_LINE],  // filename for the final result
                   *permMemPtr,  // memory to store <SD'|OP|SD> in memory
               *freePermMemPtr,  // free memory to store <SD|OP|SD> in memory
                   *windMemPtr;  // memory to store diag <SD|VEFF|SD> on files
  int
                            MJ,  // twice proj. of total ang.mom.(2*M)
                       memCase,  // = 0 - permMen, = 1 - mainly windMen
                   diagMemCalc,  // YES if <SD|OP|SD> stored in memory, otherwise NO
                       typeInt,  // type of interaction:
                                 //  (VEFF_INT, VEFF_INT3, ANG_INT, CM_INT)
               lastSD_mem_asym,  // last |initSD_asym>  of mem stored <  |V|initSD_asym>
                lastSD_mem_sym,  // last |initSD_sym>  of mem stored <  |V|initSD_sym>
             storedSD_mem_elem,  // number of stored mem <..|V|initS
              lastSD_file_asym,  // last |initSD_asym>  of file stored <  |V|initSD_asym>
               lastSD_file_sym,  // last |initSD_sym>  of filestored <  |V|initSD_sym>
            storedSD_file_elem,  // number of stored file <..|V|initS
                 tot_file_size,  // total file size for precalculated matr elem (Mbytes)
                sect_file_size,  // current available sect file size (bytes)
                  sect_file_no;  // number of last section file 
  UL
                   permMemSize,  // mem (bytes) to store diag/nondiag in memory
               freePermMemSize,  // free mem (bytes) to store diag/nondiag in memory
                   windMemSize;  // mem (bytes) to store diag/nondiag on files     

  MATR_OP
                            op;  // two- or three-particle effective interaction
} SD_STORE;

typedef  struct {
  char 
              title[ONE_LINE],  // title name for the currentshell model calculation
          type_calc[ONE_LINE];  // see choices in input data file
  int 
                     typeVeff,  // type of input effective interaction
                     run_code,  // identify the termination of Lanczos process
                storedLancVec,  // number of |lancVec> store on file
               max_iterations,  // maximum number  of lanczo iterations
                      calc_CM,  // calculate Center_of_Mass effect(yes = 1, no = 0) 
                       states,  // number of eigenstates wanted
               num_iterations,  // actual number of lanczo iterations
                      memCase,  // = 0 - permMen, = 1 - mainly windMen
                  file_lanczo;  // max file storage  available for lanczos vectors (Mb)
  double
                        *vec1, // ponters to two lanczos vectors
                        *vec2,
                    *h_matrix, // pointer to the lanczos eigenvector matrix
                   *CM_matrix, // pointer to the center_of_Mass matrix
                 *delta_eigen; // data for convergence testing
  SD_STORE
                  sdStoreVeff, // storage of <SD'|VEFF|SD> and two-part int
                   sdStoreAng, // storage of <SD'|ANGF|SD> and two-part int
                    sdStoreCM; //storage of <SD'|CM|SD> and two-part int
} LANC_PROC;

//  structure of pointers for memory in the Lanczos process

typedef   struct {
  double      *vec1,
              *vec2,
          *h_matrix,
       *delta_eigen;
} LANC_MEM;

    /*
    ** Data structure to store the eigenvalues 
    ** and the corresponding eigenvectors.     
    */

typedef   struct  {
    double
              value,
              *vector;
} EIGEN;

    /*
    **  structure for eigenvector data 
    **  from Lanczos iteration procedure
    */

typedef   struct {
     int       vecNum,
                  dim,
             numj_occ;
     double  eigenVal,
               angMom,
               val_CM;
} EIGEN_DATA;


typedef     struct {
  int
            typeSD,
            startSD,
              endSD,
              numSD,
             lastSD,
            memSize;
char
            *memPtr;
}  HEAD;

     
                  // Entrance function definitions for the Lanczos process

void input_process(char *input_file, SHELL *model);   // file: shell-model-input.c 
void single_particle_basis(SHELL *model, SP_BAS *sp_bas);  // file id-sp-basis.c
void id_slater_determinant(SHELL *model, SP_BAS *sp_bas, SD_BAS *sd_bas); 
void id_ang_mom_interaction(SP_BAS *sp_bas, MATR_OP *op_ang);
void id_two_part_intJ(char *filename, SP_BAS *sp_bas, MATR_OP *op_int);
void id_two_part_intM(char *filename, SP_BAS *sp_bas, MATR_OP *op_int);
void id_three_part_veffM(SHELL *model, SP_BAS *sp_bas, MATR_OP *op_veff);
void storeSD_MatrElemInMem(SP_BAS *spBas, SD_BAS *sdBas, SD_STORE *sdStore);
void storeSD_MatrElemOnFile(SP_BAS *spBas, SD_BAS *sdBas, SD_STORE *sdStore);
void twoPartNondiagSD_Calc(HEAD *head, SP_BAS *sp_bas,
			SD_BAS *sd_bas, MATR_OP *op_int);
void threePartNondiagSD_Calc(HEAD *head, SP_BAS *sp_bas,
			       SD_BAS *sd_bas, MATR_OP *op_veff);
void id_eigenvalue_Lanczos_search(SHELL *model, LANC_PROC *lanc_proc, 
                                    SP_BAS *sp_bas, SD_BAS *sd_bas);
void  id_matrix_vector_calc(SP_BAS *spBas, SD_BAS *sdBas,
	      SD_STORE *sdStore, double *initVec, double *finalVec);
double id_lanczo_matrix_element(SP_BAS *sp_bas, SD_BAS *sd_bas,
	      SD_STORE *sdStore, double *initVec, double *finalVec);
void id_lanczos_result(SP_BAS *sp_bas, SD_BAS *sd_bas, LANC_PROC *lancProc);

            // file: shell-model-lib.c

void   *mem_alloc(size_t, char *, char *);
void   *mem_calloc(size_t, size_t, char *, char *);
void   *mem_realloc(void *, size_t, size_t, char *, char *);
int    read_text_string(FILE *in_data, char *text_string);
int    read_int_number(FILE *in_data, int number, int *vector);
int    read_float_number(FILE *in_data, int number, double *vector);
int    read_data_string(FILE *in_data, char *model_string);
double clebsch_gordan(int j1, int j2, int j3, int m1, int m2);
TID    run_clock(double *);
void   hpcwall(double *);
void   run_time(int number, TID *run);
TID    time_step(int, int);
TID    cum_time(int clock_num, int type);
void   add_a_vector(double factor, int dim, double *init_vec, double *final_vec);
void   scale_a_vector(int dim, double factor, double *ptr);
void   append_basis_vector_to_file(char *file_namen, int n, int dim, double *vector);
void   read_basis_vector_from_open_file(char *file_name, FILE *file_ptr,
                                  int n, int dim, double *vec_ptr, float *temp_mem);
void writeLanczosEigenVec(char *fileName,EIGEN_DATA eigenData,double *vector);
void readLanczosEigenVecFromOpenFile(char *file_name,FILE *file_ptr,int vecNum,
				     EIGEN_DATA *eigenData,double *vec_ptr,float *temp_mem);
void readShellmodelEigenVec(char *filename, int vecNum, int dim, 
			    EIGEN_DATA *eigenData, double *vec_ptr,float *temp_mem);
void   eigenvalues(int dim, double *h_matrix, double *eigenvalue);
int    eigenvalue_comp(const double *one, const double *two);
void   eigenvectors(int dim, double *h_matrix, EIGEN *eigen);
int    eigenvector_comp(const EIGEN *one, const EIGEN *two);

void   eigenvalue_tred2(double **a, int n, double d[], double e[]);
void   eigenvalue_tqli(double d[], double e[], int n);
void   eigenvector_tred2(double **a, int n, double d[], double e[]);
void   eigenvector_tqli(double d[], double e[], int n, double **z);
