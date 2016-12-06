     /*
     ** The program module                                   
     **        pnshell-model.h
     ** for the proton-neutron shell model 
     ** calculation.
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
#include "bitarray.h"


#include "mpi.h"     // mpi include file



       // definition of filenames   

#define RESULT_OUTPUT            "-result.dat"
#define SD_STORE_VEFF            "-SD-veffJ"
#define LANCZO_VECTORS           "-lanc-vectors.dat"
#define SD_STORE_ANG             "-SD-J**2"
#define SD_STORE_CM             "-SD-CM"
#define STORE_MEMORY            0
#define STORE_FILE              1
#define MATR_ELEM_FILES         10  // default max file size 
#define LANC_FILES              10  // default max file size


#define DIAG                     "-diag-elem.dat"
#define NONDIAG                  "-nondiag-elem.dat"
#define EIGEN_VECTORS            "-eigen-vectors.dat"
#define CM_VECTORS               "-CM-vectors.dat"

   //  Definition of program process

#define VEFF_LANC      0
#define DIMENSION      1
#define VEFF_CALC      2
#define CM_CALC        3
#define ANG_CALC       4

#define VEFF2J_INT     0
#define VEFF2M_INT     1
#define VEFF3_INT      2
#define CM_INT         3

#define MASTER          0
#define ALLOCATE        0
#define FREE            1


// definition of sizes for memory and file allocations

#define LANC_FILES      10  // default lanczos iteration files in units of Gbytes
#define MATR_ELEM_FILES 10  // default file storage of <SD'|veff|SD> in units of Gbytes
 
#define M_BYTES         1000000 // unit trasfer: Mbytes --> bytes   
#define MAX_FILE_SIZE   1000    // in units of Mbytes
#define MEM_BLOCK_SIZE  10      // Blocksize for <SD'|OP|SD> in memory , unit Mbytes


#define   D_ZERO     0.0
#define   D_ONE      1.0
#define   TRUE       1
#define   FALSE      0
#define   YES        1
#define   NO         0
#define   CON        2
#define   PROTON     0
#define   NEUTRON    1
#define   PLUS       0
#define   MINUS      1
#define   NULL_PTR   (void *) 0
#define CLOCK_UNIT   CLOCKS_PER_SEC    // Use of special library functions   

          // max. factorials in Clebsch-Gordan coefficients   

#define FAC_LIM                300
#define ONE_LINE                80
#define ZERO_LIMIT         1.0E-10
#define MATRIX_LIMIT       1.0E-06
#define ENERGY_LIMIT       1.0E-06    // limit of energy changes through Lanczos iterartion   
#define ANG_LIMIT          1.0E-07    // limit ang mom admixture in Lanczos basis vectors     

#define   UL           unsigned int
#define   UL_ONE       1UL
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

//   GLOBAL MPI VARIABLE

int Rank;       // global parameter: the process Rank number
int NumProcs;   // global parameter: total number of processes

char   FILE_NAME[ONE_LINE],
       FILE_OUTPUT[ONE_LINE];        // global fileName for test output




typedef   struct  {   // structure definition for execution time   
  unsigned long long 
                          tick,
                           sec,
                           min,
                          hour;
} TID;


    // The basic input data for the shell model calculation
    // single-particle j-orbits in struct J_BAS   


typedef   struct { // The input spherical single-particle orbits   
  int
                N,   // N = 2*n +l = 2*osc +l
              osc,   // n osc. quantum number
                l,   // orbital ang. momemtum 
                j,   // twice total ang. momentum 
         min_part,   // min allowed particles in spherical orbit
         max_part;   // max allowed particles in spherical orbit
  double
                e;   // single-part. energy
} JBAS;


typedef    struct  {
  char
            title[ONE_LINE], // Input title for the calculation
  veffMatrElem_pp[ONE_LINE], // Effective pp-V_effective
  veffMatrElem_nn[ONE_LINE], // Effective nn-V_effective
  veffMatrElem_pn[ONE_LINE], // Effective pn-V_effective
  centerOfMass_pp[ONE_LINE], // center_of_Mass pp_matrix elements
  centerOfMass_nn[ONE_LINE]; // center_of_Mass nn_matrix elements
  int
                typeProc,  // Choice calculation process
                typeVeff,  // type of input effective interaction
                       P,  // parity  states
                      MJ,  // twice proj. of total ang.mom.(2*M)
                 part[2],  // particle number: proton[0]/neutron[1] 
             numj_orb[2],  // number of single-particle j_orbits 
                 calc_CM,  // calculate Center_of_Mass effect(yes = 1, no = 0) 
                N_low[2],  // first additional P/N harmonic osc shell
               N_high[2],  // last additional  P/N harmonic osc shell
                  delN_Z,  // max oscillator excitation for proton config
                  delN_N,  // max oscillator excitation for neutron config

                    delN,  // max energy excitation of basis |SD> states
          max_iterations,  // maximum number  of lanczo iterations
                  states;  // number of eigenstates wanted
  ULL
             file_lanczo,  // Max.total file size for lanc vectors, unit bytes
             permMemSize;  // number of <pn_SD'|OP|pn_SD> in permanent memory, unit STORE_F 
  JBAS 
                *jbas[2];  // list of single-particle j-orbits proton[0]/neutron[1]
  double
               oscEnergy;  // Harmonic oscillator energy

} SHELL;

typedef   struct {   // m-scheme quantum numbers   
  int 
              orb,   // j-orbit number - j_max first 
                N,   // N = 2*n +l = 2*osc +l
              osc,   // oscillator number  n 
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
	 ** Possible limitations on particle number
	 ** in j-orbits are included
	 */

typedef    struct {
  int
              part,   // number of particle
          numj_orb,   // number of j-orbits 
          numm_orb;   // number of m-orbits 
  JBAS 
             *jbas;  // list of single-particle j-orbits 
  MBAS
             *mbas;   // list of m-orbits
  MASK
              mask;   // structure to specify reduction/masking
                      // of number of particles in j-orbits  
} SP_BAS;

  typedef   struct {
           int
                       a,    /* left orbit  - particle one */
                       b,    /* left orbit  - particle two */
                       c,    /* right orbit - particle one */
                       d;    /* right orbit - particle two */
           UL 
                    left,    /* left two-particle config   */
                   right;    /* right part-particle config */
   } ORB_ID;

        /*
	** identify orbits in two- or three-particle 
	** case for elements in m-scheme table for 
	** identical particles  
	*/

typedef struct {
  ULL          one, // identify occ particle orbits
               two; // Pauli phase   

  double    val[2];  // two-part VEFF and CM matr elem      
} ID_INT;

        /*
	** identify orbits in two- or three-particle 
	** case for elements in m-scheme table for 
	** proton neutron case 
	*/

typedef struct { // identify orbits in two- or three-particle 
                 // case for elements in m-scheme table: 
  ULL       one[2], // identify occ particle orbits
            two[2]; // Pauli phase   
  double    val[2];  // value of matrix elements      
} PN_INT;


typedef struct  {  // pointers to table of two-body matrix elements   
  int    num_nondiag_elem,  // num reserved nondiag elem
           maxNondiagElem;  //r max num <..|op|ij>,
  double     	 *id_diag,  // list of diagonal <ij|V|ij> elements,
              *id_CM_diag;  // list of diagonal <ij|CM_int|ij> elements,
  ID_INT 
          *id_nondiag_bas,  // pointer to start of nondiag ID_INT elements
                            // <xx|V|kl> elements for each (kl) pair
       ** id_nondiag_table; // kl_table of pointer to nondiag matrix elements 
} MATR_ID;

typedef struct  {  // pointers to table of two-body matrix elements   
  int     numNondiagElem;  // num reserved nondiag elem
  //     maxNondiagElem;  // max num <..|op|ij>,
  double         *pn_diag,  // list of diagonal <ij|V|ij> elements,
              *pn_CM_diag;  // list of diagonal <ij|CM_int|ij> elements,
  PN_INT       ***kl_list,  // ptr to list of kl_table for each del_m and del_p
             *nondiag_bas;  // pointer to start of nondiag ID_INT elements   
                            // <xx|V|kl> elements for each (kl) pair,     
} MATR_PN;

typedef struct {   // Data structure to write/read  <SD'|OP|SD> to/from file
  float
                   value;   // value of <SD'|OP|SD>
  ULL
                   final;   // number of final |SD'>
} STORE_F;

typedef   struct {
  int        startSDnum, // absolute position of first |SD> in current block
               totNumSD, // total number of |SD> in current block
               endSDnum; // absolute position of last |SD> in current block
} SD_BLOCK;


typedef     struct {
  char
               title[ONE_LINE], // Input title for the calculation
      pnSectFilename[ONE_LINE];

  int
             typeProc, // Choice calculation process
             typeVeff, // VEFF2, VEFF3, CM
                    P, // parity  states
                   MJ, // twice proj. of total ang.mom.(2*M)
              calc_CM, // calculate Center_of_Mass effect(yes = 1, no = 0)
                 parZ, // possible parities of proton |SD(Z)>
                       // coding: -1: parZ  -1; both parZ = 0; +1: parZ = +1
               num_gr, // total number of groups
            maxOscN_Z, // Max 0scillator N for the proton config
            maxOscN_N, // Max 0scillator N for the neutron config
              refOscN, // Summed Oscillator N for the lowest P + N config
              maxOscN, // Max oscillator N for the total P/N config
             N_low[2], // first additional P/N harmonic osc shell
            N_high[2], // last additional  P/N harmonic osc shell
               delN_Z, // max oscillator excitation for proton config
               delN_N, // max oscillator excitation for neutron config
                 delN, // max oscillator excitation for the total P/N basis |SD> states
           totNumSD_Z, // total number of |SD(Z)>
      maxGroupNumSD_Z, //  
           totNumSD_N, // total number of |SD(N)
      maxGroupNumSD_N, //  
       RankTotNumSD_N, // total number of |SD(N)> in each Rank
          totNumSD_ZN, // total number of |SD(Z),SD(N)>
   totNumSD_ZN_oscLim, // same number with osc. limitation
      RankTotNumSD_ZN, // total number of |SD(Z),SD(N)> for each Rank
RankTotNumSD_ZN_oscLim,// same number with osc. limitation
             maxSD[2], // max |SD()> within a group - proton [0]/ neutron [1]
                       // in a group
         pnSectFileNo; // file no for nondiag pn matr elem
  ULL
          permMemSize, // size to store <pn_SD'|OP|pn:SD>  in memory, unit STORE_F
       permMemSizeANG, // size to store <SD'|ANG|SD>  in memory, unit STORE_F
     curr_permMemSize, // free memory space
 pnNondiagTotFileSize, // total file size for nondiagMatrElem in STORE_F
       pnSectFileSize, // current pnSectFile size in STORE_F
  curr_pnSectFileSize, // available STORE_F elem in current pnSectFile
    totMemNondiagElem, // number of stored <SD'|OP|SD> in memory
   totFileNondiagElem; // number of stored <SD'|OP|SD> on files
 
  float
         *diagMemElem, // total memory for pnDiagElem
      *diagMemElemANG; // total memory for pnDiagElem

  STORE_F
            *blockMem, // block mem Ptr to <SD'|OP|ST>
             *permMem, // perm mem Ptr to  <SD'|OP|SD>
          *permMemANG, // perm mem Ptr to  <SD'|ANG|SD>


        *curr_permMem; // available perm mem 
  FILE
      *pnSectFilePtr;  // file pointer to current pn section file

  SP_BAS 
           *spBas[2];  // complete set of single-particle basis
  MATR_ID
               op_pp,  // two-particle interaction VEFF or ANG
               op_nn;

  MATR_PN
               op_pn;

  SD_BLOCK
             *sdBlock;  // block struc for lancVec> in each process 
  double
            oscEnergy;  // Harmonic oscillator energy
} GR_BAS;                          
   
// Data information  for a single proton/neutron group

typedef   struct   {
  int                    // for OP = VEFF and ANG
                  grNum, // current group number
                   m[2], // ang.mom. projection for |SD> proton[0]/neutron[1]
                 par[2], // parity P_Z(= +1,-1) for |SD> for |SD>
                         // proton[0]/neutron[1] 
               numSD[3], // number of |SD> proton[0]/neutron[1]
            numAmplStep, // num init/final ampl steps in gurr group
               currSD_Z, // current start |initSD_proton>
               currSD_N, // current start |initSD_neotron>
               numSD_ZN, // total number of |SD> in curr group 
              currSD_ZN, // current start |initSD_proton/neutron>
        numSD_ZN_oscLim, // same number with osc. limitation
           RankNumSD_ZN, // total number of |SD(Z),SD(N)> for each Rank
    RankNumSD_ZN_oscLim, // same number with osc. limitation
                *occ[2], // lists to store occ orb in all |SD(Z)> and SD(N)>
    numSingleDataBlocks, // number of SingleGroup pn data blocks
     numMultiDataBlocks, // number of multiGroup pn data blocks
      sizeLastDatablock, // number of last SingleGroup data block
        singleDataBlock; // number of STORE:F elem in last MultiGroup data block
  float
           *diagMemElem; // pointer to diagMemElem in curr group
  ULL
               startAmp, // start position for Lanczos vector amplitudes
          storedSD_Elem, // total  STORE_F nondiag matr elem in curr group
           blockNumElem, // num of STORE_F elem in a single block
       lastMultiNumElem, // num of STORE_F elem in last MultiGroup
                 *SD[2]; // ptr to list of |SD>     for |SD(Z)>|SD(N)>
  MASK
                mask[2]; // structure to specify reduction/masking
                         // of number of particles in proton[0] and
                         // neutron[1] j-orbits  
  BARR_ELTYPE
              *bitArray; // Array of bits 
} GROUP;

typedef  struct {
  char              
    sectFilename[ONE_LINE]; 
  int
                  typePart,
                   typeInt,
                   partNum,
                   nummOrb, 
                     numSD,
                    currSD,
                      *occ,
             storedSD_Elem,
               lastNumElem,
                sectFileNo,
             numDataBlocks;

  ULL
               tempMemSize, // unit STORE_F
          currSectFileSize,
               totFileSize,
                       *SD;
STORE_F
                  *tempMem;
FILE
              *sectFilePtr;

  MATR_ID 
                    *id_op;
  SP_BAS            *spBas;

} ID_STORE;

typedef  struct {
  char        title[ONE_LINE];  // path to scratch directory
  int                run_code,  // identify the termination of Lanczos process
           totNumSD_ZN_oscLim,  // total number of |SD(Z),SD(N)>
                storedLancVec,  // number of |lancVec> store on file
               max_iterations,  // maximum number  of lanczo iterations
                       states,  // number of eigenstates wanted
               num_iterations,  // actual number of lanczo iterations
          pnLanczoTotFileSize; // total file size for LancVec in float

  double                *vec1, // ponters to two lanczos vectors
                        *vec2,
                    *h_matrix, // pointer to the lanczos eigenvector matrix
                 *delta_eigen; // data for convergence testing
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
          numj_occ[2];
     double  eigenVal,
               angMom,
                valCM;
} EIGEN_DATA;


typedef     struct {
  int
            typeSD,
            startSD,
              endSD,
         numSD_Elem,
             lastSD,
            memSize;
STORE_F     *memPtr;
}  HEAD;


        /*
	** data structure for calculating 
	** and storing nondiagpp- and nn-
	** matrix elements in functions
	**    pnCalcSD_NondiagMatrElem()
	** and 
	**    idTWOPartNondiagSD_Caklc() 
	*/  

  typedef    struct {
    int    numMemElem,
         numTotalElem;
  } ID_DATA;


// Entrance function definitions for the Lanczos process


// file: pn-shell-model-input.c 
int input_process(char *input_file, SHELL *model);

// file: pn--sp-basis.c
void pn_single_particle_basis(GR_BAS *grBas);

// file: pn-SD-basis.c
int pn_slater_determinant(GR_BAS *grBas, GROUP **group);

// file id-two-part-intJ.c
void id_two_part_intJ(char *filename,SP_BAS *spBas, int calc_CM,
                                                MATR_ID *op_int);

// file id-two-part-angJ.c
void id_two_part_angJ(SP_BAS *spBas,MATR_ID *op_ang);

// file pn-two-part-intJ.c
void pn_two_part_intJ(char *filename, int parZ, SP_BAS *spBasZ,
                     SP_BAS *spBasN,int calc_CM,MATR_PN *op_pn);
// file:  pn-two-part-iangJ.c
void pn_two_part_angJ(int parZ, SP_BAS *spBasZ,SP_BAS *spBasN, 
		                             MATR_PN *op_ang);
// file PAR-pnStoreSD_MatrElem.c
void pnStoreSD_MatrElem(int typeCalc,int calcCM,GR_BAS *grBas,
                                                 GROUP *group);

// file: SINGLE-idNondiagCalc.c
int idTwoPartNondiagSD_Calc(int type, int typeCM,GR_BAS *grBas,
                             GROUP *group,STORE_F *idMemBlockPtr);
int idThreePartNondiagSD_Calc(ID_STORE *store, HEAD *head);

void pnMatrixVectorCalc(int typeInt,GR_BAS *grBas, GROUP *group,
			double *initVec,double *finalVec);

// file: SINGLE-pnLancProcess.c
void pnEigenvalueLanczoSearch(GR_BAS *grBas,GROUP *group,
			             LANC_PROC *lancProc);
// file: SINGLE-pnLancResult.c 
void pnLancResult(GR_BAS *grBas, LANC_PROC *lancProc);

            // file: shell-model-lib.c

void   *mem_alloc(size_t, char *, char *);
void   *mem_calloc(size_t, size_t, char *, char *);
void   *mem_realloc(void *, size_t, size_t, char *, char *);
int   
 read_text_string(FILE *in_data, char *text_string);
int    read_int_number(FILE *in_data, int number, int *vector);
int    read_float_number(FILE *in_data, int number, double *vector);
int    read_data_string(FILE *in_data, char *model_string);
double clebsch_gordan(int j1, int j2, int j3, int m1, int m2);
TID    run_clock(double *);
void   hpcwall(double *);
void   run_time(int number, TID *run);
TID    time_step(int, int);
TID    wallClock(int clock_num, int type);
TID    cpuClock(int clockNum, int type);
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
