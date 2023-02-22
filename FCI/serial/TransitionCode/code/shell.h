    /*
     * The program module                              
     *                      shell.h                    
     * for proton - neutron the Lanczo transition program package.
     */

     /**  Standard ANSI-C include files **/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h> 
#include <time.h>


#define CLOCK_UNIT    CLOCKS_PER_SEC    /* Use of special library functions */

#define   UL         unsigned int
#define   UL_ONE     1u
#define   UL_ZERO    0u
#define   D_ZERO     0.0
#define   D_ONE      1.0
#define   TRUE       1
#define   FALSE      0
#define   YES        1
#define   NO         0
#define   FINISHED   2
#define   NULL_PTR   (void *) 0

          /* max. factorials in Clebsch-Gordan coefficients */

#define FAC_LIM                300
#define ONE_LINE                80
#define MAX_PARTICLES          100
#define ZERO_LIMIT         1.0E-10
#define MATRIX_LIMIT       1.0E-05

     /* Macro definitions for integer arguments only */

#define   MAX(a,b)   ( ((a) > (b)) ? (a) : (b) )
#define   MIN(a,b)   ( ((a) < (b)) ? (a) : (b) )
#define   PARITY(a)  ( (a) % 2) 
#define   PHASE(a)   (1 - 2 * (abs(a) % 2))
#define   MOD(a,b)   ((a) - ((a)/(b)) * b)
#define   SIGN(a,b) ((b)<0 ? -fabs(a) : fabs(a))

    /* Definition of constants in memory allocations */

#define  INPUT_PROTON_JBAS                     0
#define  INPUT_NEUTRON_JBAS                    1
#define  PROTON_M_ORBITS                       2
#define  NEUTRON_M_ORBITS                      3
#define  PROTON_NEUTRON_GROUPS                 4
#define  PROTON_CONFIGURATION                  5
#define  NEUTRON_CONFIGURATION                 6
#define  INFO                                  7
#define  IDENTICAL_CONFIGURATION               8
#define  TRANSITION_MATRIX_ELEMENTS            9
#define  M1_E2_E3_MATRIX_ELEMENTS             15
#define  E2_TRANSITION                        16
#define  E3_TRANSITION                        17
#define  M1_TRANSITION                        18
#define  GROUP_CHANGE                         19

#define  DOUBLE_INT    YES    // 64 bits version (YES) 32 bits version (NO)  

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




    /* structure definition for execution time */

   struct  tid {
      int     tick,
               sec,
               min,
              hour;
   };

    /* Basic input data to define the shell model problem. */

    /* The spherical single-particle orbits */

   struct  jbas {
      int     osc,   /* osc. quantum number              */
	        l,   /* orbital ang. momemtum            */
	        j,   /* twice total ang. momentum        */
         min_part,   /* min particles in spherical orbit */    
         max_part;   /* max particles in spherical orbit */
      char status;   /* = n: no coupling to this orbit, normal = a */
      double    e;   /* single-part. energy                        */
   };
     /* Specifications for one-particle transitions */

   struct trans      {
      int   num_init,  /* position of initial state   */
              j_init,  /* 2J for the initial state    */
           num_final,  /* position of the final state */
             j_final;  /* 2J for the final state      */
   };

     /* The shell model data of the many-body system */

   struct shell {
      char                P,  /* parity  states                        */
                  J_type[5],  /* = even (odd) for even(odd) J-state    */
             calc[ONE_LINE];  /* e2, e3 or m1 transition               */
      int                 Z,  /* proton number                         */
                          N,  /* neutron number                        */
                         MJ,  /* twice proj. of total ang.mom.(2*M)    */
                 same_basis,  /* if YES, init. and final basis are same*/
                     numj_Z,  /* num. of proton spher. orbits          */
                     numj_N,  /* num. of neutron spher. ornits         */
                     memory,  /* memory for <SD_n|H|SD_m> in iter.proc.*/
                  num_trans;  /* number of one-particle transitions    */
      struct  jbas   *jbasZ,  /* ptr. to proton spher. orbits          */
                     *jbasN;  /* ptr. to neutron spher. orbits         */
     double        e2_n_eff,  // neutron effective E2 charge
                   e2_p_eff,  // proton  effective E2 charge
                m1_n_eff_gl,  // neutron  effective gl charge
                m1_p_eff_gl,  // proton  effective gl charge
                m1_n_eff_gs,  // neutron  effective gs charge
                m1_p_eff_gs;  // proton  effective gs charge
      struct trans   *trans;  /* ptr to transition matrix elements     */
   };
    
              /* Storage of the m-scheme single-particle orbits. */

              /* m-scheme quantum numbers */

   struct mbas  {
      int     orb, /* number acc. to input           */
	      osc, /* oscillator number              */
	        l, /* orbital ang. momentum          */
	      par, /* parity of orbit (+1,-1)        */
	        j, /* twice total spin               */
  	        m, /* twice proj. of total spin      */
	    phase, /* "Whitehead" phase factor       */   
	   status; /* active(no active) orbit = 1(0) */
      double    e; /* single-particle energy         */
   };
              /*  storage  of the m-orbit   */

   struct m_orbit {
      int            num;  /* number of m-orbits  */
      struct mbas  *mbas;  /* storage of m-orbits */
      int        numMask,   // number of reduced j-orbits
                    *lim;   // upper and lower particle limits for j-orbits
      ULL          *list;   // list of j-orbits mask structure
   };

          /* Storage of FILE pointers for the calculation. */

   struct files    {
      FILE          *in_data,   /* shell model in_data              */
                   *out_data,   /* resulting out_data               */
                     *vector;   /* storage of lanczo vec.           */
   };

        /*
         * Data structure for organization of SD used in both
         * identical particle and proton-neutron lanczo procedure.
         */
   struct         info_com {
      int   type,  /* proton only(= 0), neutron only(= 1), both ( >= 2 ).*/
             par,  /* parity(+1,-1,0 (both)) of SD for the total system.*/ 
              MJ,  /* twice proj. of total ang.mom. */
             num,    /* number of SD for the total system.*/
/*************************
       num_eigen,   dimension of eigenvectors for the total system.
***************************/
          J_even,  /* = YES(NO) for even(odd) J-value for eigenstates */
          tr_sym;  /* = YES(NO) for the case of time-reversal symmetry */
   };
        /* 
         * Data structure for organization of  SD       
         * used in the identical particle lanczo procedure.
         */

   struct     info_id  {   /* data for the total system of  SD.*/
      int      num[2], /* num of asym SD (num[0]) and sym SD (num[1]).*/
                 part, /* number of identical particles.*/
                  orb; /* number of single-particle orbits.*/
   };        
        /* 
         * Data structure for organization of  SD       
         * used in the proton-neutron Lanczo procedure.
         */

   struct     info_pn  {   /* data for the total system of  SD.*/
      int      num_gr,   /* number of groups with both +1 and -1 parity.*/
               Z_part,   /* number of proton particles.*/
                Z_num,   /* total number of proton SD.*/
             max_Z_SD,   /* max. number of proton SD in a group.*/  
                Z_par,   /* possible parities of the proton SD.*/
                         /* coding: -par = -1; both par. = 0; +par = +1.*/
                Z_orb,   /* number of proton single-particle orbits.*/
            Zm_max[2],   /* max. m-values for + and - parity */
            Zm_min[2],   /* min. m-values for + and - parity */
               N_part,   /* number of neutron particles.*/
                N_num,   /* total number of neutron SD.*/
             max_N_SD,   /* max. number of neutron SD in a group.*/  
                N_par,   /* possible parities of the neutron SD.*/
                         /* coding: -par = -1; both par. = 0; +par = +1.*/
                N_orb,   /* number of neutron single-particle orbits.*/
            Nm_max[2],   /* max. m-values for + and - parity */
            Nm_min[2];   /* min. m-values for + and - parity */
   };

   struct       group {
      int        Z_m,      /* ang.mom. projection M_Z for proton SD*/
               Z_par,      /* parity P_Z for proton SD.*/
             Z_start,      /* start num. of proton SD in SD_Z[].*/
             Z_numSD,      /* number of proton  SD.*/
                 N_m,      /* ang.mom. projection M_N for neutron SD.*/
               N_par,      /* parity P_N for neutron SD.*/
             N_start,      /* start num. of neutron SD in SD_N[].*/
             N_numSD,      /* number of neutron SD.*/
          ampl_start;      /* start number of ampl. in lanczo vectors.*/
   };       
            /* data structure for temporary storage of SD */
  
   struct    sd    {
      ULL     config;
      int   m_value,
             parity;
   };

     /* structure to store start positions of proton/neutron SD i each group.*/

   struct sd_pos  {
      int     p_start,
               p_stop,
              n_start,
              n_stop;
            
   };  

        /* Storage of the off-diag one-particle matrix elements */

   struct trans_nondiag {
      double  val;             /* value of the matrix elements       */
      ULL      one,             /*  a ONE in orbit j                  */
           tr_one,             /*  a ONE in time-reversed orbit of j */
              two;             /*  ONE's between l and j             */
/*********************************
      int   phase;              (-)**(j_f - j_i)
****************************/
   };
        /* Index tables to one-particle off-diag. matr. elem */

   struct table_l  {
       double                    diag; /* diag. single-particle matr. elem.*/
       struct trans_nondiag  *nondiag; /* ptr to the off-diag. matr.elem. */
   };

      	/* function definition for trans-1.c */

void basic_input_data(void);
void type_of_calculation(void);
void memory_allocation(int,int);
int max_non_diag_elem(int);
void matrix_elements_in_m_scheme(int, int);
double m_scheme_e2_diag(int, int, struct mbas *);
double m_scheme_e2_nondiag(int, int, int, struct mbas *);
double rad_int2(int, int, int, int);
double m_scheme_m1_diag(int, int, struct mbas *);
double m_scheme_m1_nondiag(int, int, int, struct mbas *);
double m_scheme_e3_nondiag(int, int, int, struct mbas *);
double rad_int_lambda(int, int, int, int ,int);

        /* function definition for id-transition.c  */

void id_calculation(void);
void id_slater_determinant(void);
void number_of_identical_SD(int, int);
void number_of_SD_tr_sym(int, int, struct m_orbit);
void number_of_SD_no_sym(int, struct m_orbit);
void identical_SD_configuration(int, int);
void identical_SD_tr_sym(int, int, struct m_orbit);
void identical_SD_no_sym(int, struct m_orbit);
void id_particle_transitions(int);
void id_diag_trans(int, int, struct table_l *, ULL *, double *);
void read_id_eigenvector(char *, int, int, int, double *);
void id_trans_contribution(int, int, int, int, int, int, struct table_l *, double *, double *);
void diag_contr(int, int, double *, double *, double *);
void id_nondiag_asymI(int, int, int, int, int, struct table_l *, double *, double *);
void id_nondiag_symI(int, int, int, int, int, struct table_l *, double *, double *);
void id_nondiag_asymII(int, int, int, int, struct table_l *, double *, double *);
void id_nondiag_symII(int, int, int, int, struct table_l *, double *, double *);
void id_nondiag_nosymIII(int, int, int, int, int, struct table_l *, double *, double *);
void id_nondiag_nosymIV(int, int, int, int, struct table_l *, double *, double *);
double overlap(int, int, double *, double *);


        /* function definition for pn-transition.c  */

void pn_calculation(void);
void pn_slater_determinant(void);
void parity_check(int);
void number_of_groups(int);
void group_m_limit(int, int, int, struct m_orbit, int *);
int number_of_proton_SD(int);
void proton_SD_configuration(int);
int number_of_neutron_SD(int);
void neutron_SD_configuration(int);
void pn_particle_transitions(int);
void pn_diag_matrix_elements(void);
void read_pn_eigenvector(char *,int, int, double *);
void group_change(int, int, int, int, int, int, int  *);
void pn_nondiag_iteration(double *, double *);
void nondiag_proton_contribution(int, int, ULL *, int, ULL *, int, double *, double *);
void nondiag_neutron_contribution(int, int, ULL *, int, ULL *, int, double *, double *);
double pn_overlap(int, double *, double *);

      	/* function definition for trans-lib.c */

void  *vector(int, int);
void free_vector(void *);
void  *vector_realloc(void *, int);
void  **matrix(int, int, int);
void free_matrix( void **, int);
void read_orbit(int, struct jbas *);
int j_comp(const struct jbas *, const struct jbas *);
int read_line(char *, char *);
void read_numbers_into_vector(char *, int, int *);
void removeBlanks(char *string);
void single_m_states(int, struct jbas *, struct m_orbit *);
void part_orb_limitations(int num_j, struct jbas *jbas, 
                                    struct m_orbit *orbit);
void new_single_m_states(int, struct jbas *, struct m_orbit *);
int mbas_parity(const struct mbas *, const struct mbas *);
void read_transitions(void);
int trans_comp(const struct trans *, const struct trans *);
int plus_parity(const struct sd *, const struct sd *);
int minus_parity(const struct sd *, const struct sd *);
int plus_m_value(const struct sd *, const struct sd *);
int minus_m_value(const struct sd *, const struct sd *);
int config_comp(const struct sd *, const struct sd *);
int binomial(int, int);
void run_time(int, struct tid *);
double clebsch_gordan(int, int, int, int, int);
double fac_ratio(int, int);
