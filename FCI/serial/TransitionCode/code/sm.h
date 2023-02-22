          /*
	   * The program module                           
	   *                sm.h                     
	   * contains all global variable definitions for the one-particle 
	   * transformation program package.
	   */

	     /* Storage of file names for the calculation.*/

char 
     File_output[ONE_LINE],       /* output data                          */
     File_eigen_init[ONE_LINE],   /* storage file of initial eigenvectors */
     File_eigen_final[ONE_LINE];  /* storage file of final eigenvectors   */


          /* Basic input data to define the shell model problem. */

struct shell Model[2];

          /* Storage of the m-scheme single-particle orbits for protons and neutrons.*/

struct m_orbit      Z_orbit,    /* proton particle type */
		    N_orbit;    /* neutron particle type */

          /* Storage of FILE pointers  for the calculation.*/

struct files    Files;

          /* 
           * Permanent storage table containing all informations 
           * of the organization of the complete set of SD.
           */

struct info_com	    Info_com[2];         /* for all cases */
struct info_id      Info_id[2];         /* for the identical particle case.*/
struct info_pn      Info_pn[2];         /* for the proton-neutron case.*/

struct group        *Group[2];

int                 *Group_change_p,     /* identify change of group for OP(p,p') */
                    *Group_change_n;     /* identify change of group for OP(n,n') */

          /* 
           * Permanent storage for the  proton and neutron SD for both 
           * initial(case = 0) and the final case(case = 1).
           */

ULL		 *SD_Z[2],  /* SD for protons in proton/neutron case       */
                 *SD_N[2],  /* SD for neutrons in proton/neutron case      */
              *SD_asym[2],  /* asym. SD identical particle case            */
           *SD_tr_asym[2],  /* time-rev. asym. SD identical particle case  */
               *SD_sym[2];  /* sym. SD identical particle case,            */
           
	/* Identification of passive m-orbits */

ULL      Passive_orbit;

double      *SD_diag,   /* Storage of diagonal  matrix <SD|H|SD>  */
          *Eigen_vec,   /* Storage for initial/final eigen vector */
        *Trans_contr;   /* Storage for OP*|Eigen_vec >            */

       /* 
        * Storage of the one-particle matrix elements and index table 
        * for proton ( = 0) and neutron ( = 1) particles.
        */

struct table_l        *Table_pn_l[2];         /* index tables */

struct trans_nondiag  *Trans_pn_nondiag[2];   /* nondiag. matrix elements */

       /* Storage of the transition matrix elements between eigenstates. */

double     *Transitions;              





