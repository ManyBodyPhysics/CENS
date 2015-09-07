

  /*
  **          Shell-model proton-neutron code 
  **                 August 2012
  ** Shell model energy spectrum program based on an m-scheme
  ** basis of the proton-neutron slater determinants 
  **       |SD(M_Z), SD(M_N): M = M_Z + M_N> 
  **
  */

#include "PAR-pnShellModel.h"

// **************  local function declarations

static void  pn_LancIterateCalc(SHELL *model);
    /* 
    ** performs a Lanczos iteration process for a system
    ** of protons and neutrons to obtain
    ** a limited number of shell model eigenstates.
    ** The result is stored on file
    */

static void  shellModelBasis(SHELL *model,GR_BAS *grBas,GROUP **group);
   /*
   ** calculates 
   **   1. The complete set of j-scheme 
   **      and m-scheme single-particle states  
   **   2. The complete many-particle slater 
   **      determinants in m-scheme
   **      |SD(M,P)> = |SD(Z:M_Z,P_Z), SD(N:M_N,P_N)
   **      separated into groups for fixed M_Z and M_P
   */

static void blockSD_struct(GR_BAS *grBas);
   /*
   ** calculates and return in SD_BLOCK sdBlock 
   ** number of |SD_init> for each node
   */

static void pn_print_group_data(char *file_name, GR_BAS  *gr_basis,
                                                      GROUP  *group);
   /*
   ** writes to output file all calculated group data.
   */

static void addID_1to2SinglePartTermVeff(int num_part, SP_BAS *spBas,
					 int calc_CM, MATR_ID *op_int);
     /*
     ** adds contributions from single-particle terms to the
     ** effective two-particle matrix elements  <k.l |VEFF| k.l>
     ** and if(calc_CM == 1) <k.l |CM_INT| k.l>
     */

//********  End local function declarations ******
 
           //  *** The function definitions  ***

int main(int argc, char *argv[])
{
  char     inputFile[ONE_LINE]; 
  SHELL    model;
  TID      wallTime, cpuTime;
  int      calcCM;



  //  MPI initialization

  MPI_Init(&argc, &argv);      
  MPI_Comm_size(MPI_COMM_WORLD,&NumProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &Rank);

  // Read and open input filename containing basic input data

  if(argc >= 2) { // as program parameters
    strcpy(inputFile, argv[1]); // save in_data filename
  }
  else  { // or as separate typed in filename
    printf("\n Type input filename for basic shellModelData = ");
    scanf("%s",inputFile);
  } // end of data fileName input

  // Time control for the total shell model calculation

  wallClock(1,0); // initialization 
  cpuClock(1,0);
  wallClock(1,1); // start time
  cpuClock(1,0);

  if(input_process(inputFile, &model) == YES) {    //  shell-model-input.c

    pn_LancIterateCalc(&model);

    // Stop and read total process time 

    wallClock(1,2); // stop time
    cpuClock(1,2);
    wallTime = wallClock(1,3);  // read time
    cpuTime  = cpuClock(1,3);

    if(Rank == MASTER) {
      char  filename[ONE_LINE]; 
      FILE  *file_ptr;
      
      sprintf(filename,"%s%s",model.title,RESULT_OUTPUT); 
      if((file_ptr = fopen(filename,"a"))== NULL) {
	printf("\n\nRank%d: Error in function lanc-main.c():",Rank);
	printf("\nWrong file = %s for the output data\n\n", filename);
	MPI_Abort(MPI_COMM_WORLD,Rank);
      }
      fprintf(file_ptr, "\n\nTotal shell model process wall time used");  
      fprintf(file_ptr, " %llu hour %llu min %llu sec",
	      wallTime.hour, wallTime.min, wallTime.sec);
      fprintf(file_ptr, "\nTotal shell model process cpu time used");  
      fprintf(file_ptr, " %llu hour %llu min %llu sec\n",
	      cpuTime.hour, cpuTime.min, cpuTime.sec);
      fclose(file_ptr);

    } //end MASTER

  }  //end test  total input Data

  printf("\n\n");

  MPI_Finalize();

  return 0; // sucsessfull termination

} // End: function main()
 
    /* 
    ** The function
    **     pn_LancIterateCalc()
    ** performs a Lanczos iteration process for a system
    ** of protons and neutrons to obtain
    ** a limited number of shell model eigenstates.
    ** The result is stored on file
    */
  
static void  pn_LancIterateCalc(SHELL *model)
{
  char          *func = {" pn_LancIterateCalc(): "};
  int           calcCM;  // center-of Mass in use ??
  GR_BAS        grBas;
  SP_BAS        spBas[2];
  GROUP         *group = NULL;
  LANC_PROC     lancProc;

   calcCM = 0;

  // block structure of lanczo vectors for each process
 
    grBas.sdBlock 
      = MALLOC(NumProcs,SD_BLOCK, func,"grBas->sdBlock[]");

 // Pointer ref of single-particle proton/neutron data

  grBas.spBas[0] = &spBas[0];  
  grBas.spBas[1] = &spBas[1];

  //  Memory to store interaction SD_MatrElem

  grBas.permMemSize  = model->permMemSize;
 
   /*
   ** calculate 
   **   1. The complete set of j-scheme 
   **      and m-scheme pn-single-particle states  
   **   2. The complete pn-many-particle slater 
   **      determinants in m-scheme 
   **      |SD(M,P)> = |SD(Z:M_Z,P_Z), SD(N:M_N,P_N)
   */

  shellModelBasis(model, &grBas, &group);


  /***** Calculation of shell model dimension only ******/

  if(model->typeProc == DIMENSION) return;

  // block structure of lanczo vectors  // Output ??
 
  blockSD_struct(&grBas);

     /*
     ** Read from input all veff matrix elements and store 
     ** the complete set m-scheme veff in op_veff.
     ** Add single-particle energies to op_veff
     ** Calculate and store diagonal and nondiagonal slater
     ** determinant matrix elements <SD'|VEFF|SD> in memory 
     ** or on files according to specifications in 
     ** MPI-shell-model.h
     */

  switch(model->typeVeff) {
      case 0:
              /* 
	      ** For model->part[PROTON] >= 2: 
	      ** Read from input all J_coupled  two-particle 
	      ** proton-proton matrix element and store in 
	      ** grBas->op_VEFFpp the complete set of m-scheme
	      ** two-particle matrix elements.
	      */

	if(model->part[0] >= 2) {
	  id_two_part_intJ(model->veffMatrElem_pp,grBas.spBas[0], 
                                      grBas.calc_CM,&grBas.op_pp);

	  addID_1to2SinglePartTermVeff(model->part[0],grBas.spBas[0],
	  		       calcCM,&grBas.op_pp);

        }     
              /* 
	      ** For model->part[NEUTRON] >= 2: 
	      ** Read from input all J_coupled  two-particle 
	      ** neutron-neutron matrix element and store in 
	      ** sdStore->op_VEFFpp the complete set of m-scheme
	      ** two-particle matrix elements.
	      */

	if(model->part[1] >= 2) {
	  id_two_part_intJ(model->veffMatrElem_nn,grBas.spBas[1],
			              grBas.calc_CM,&grBas.op_nn);

	  addID_1to2SinglePartTermVeff(model->part[1],grBas.spBas[1],
	  		       calcCM,&grBas.op_nn);

	}
              /* 
	      ** Read from input all J_coupled  two-particle 
	      ** proton-neutron matrix element and store in 
	      ** sdStore->op_VEFFpn the complete set of m-scheme
	      ** two-particle matrix elements.
	      */

	pn_two_part_intJ(model->veffMatrElem_pn,grBas.parZ,grBas.spBas[0],
                                grBas.spBas[1],grBas.calc_CM,&grBas.op_pn);
	break;
      case 1:
              /* 
	      ** For model->part[PROTON] >= 2: 
	      ** Read from input all m_scheme two-particle 
	      ** proton-proton matrix element and store the 
	      ** result in sdStore->op_VEFFpp
	      */

	if(model->part[0] >= 2) {
	  // id_two_part_intM(model->veffMatrElem_pp,grBas.spBas[0],
          //                        grBas->calcCM,&grBas->op_VEFFpp);
	}
              /* 
	      ** add proton single-particle energies to the matrix
	      ** elements for the current particle number
	      */

	// add1to2SinglePartTermVeff(model->part[0],grBas.spBas[0],
        //                                          &grBas->op_VEFFpp);

              /* 
	      ** For model->part[NEUTRON] >= 2: 
	      ** Read from input all m_scheme two-particle 
	      ** neutron-neutron matrix element and store the 
	      ** result in sdStore->op_VEFFpp 
	      */

	if(model->part[1] >= 2) {
	  // id_two_part_intM(model->veffMatrElem_nn, grBas.spBas[1],
          //                               grBas->calcCM,&grbas->_nn);
	}
              /* 
	      ** add neutron single-particle energies to the matrix
	      ** elements for the current particle number
	      */

	// addID_1to2SinglePartTermVeff(model->part[1],grBas.spBas[1],
        //                                          &grBas->op_VEFFnn);

              /* 
	      ** Read from input all m_scheme two-particle 
	      ** proton-neutron matrix element and store the 
	      ** result in sdStore->op_VEFFpn
	      */

//     not implemented
//     pn_two_part_intM(model->veffMatrElem_pn, grBas.parZ,grBas.spBas[0],
//                                         grBas.spBas[1], sdStore.op_VEFFpn);
	break;
//      case 2: 

//      not implemented
//	break;

      default: 
         printf("\n\nRank%d: Error in function main(id_LancIterateCalc)",
                                                                   Rank);
         printf("\nWrong Veff interation value - typeVeff = %d\n\n",
                                                         model->typeVeff);
	 MPI_Abort(MPI_COMM_WORLD,Rank);

  }  // end switch()

     /*
     ** Calculate and store in memory all 
     ** diagonal proton-neutron matrix elements
     **     <SD(Z)g,p;SD(N)g,n|OP|SD(Z)g,p;SD(N)g,n>
     ** and as many as possible nondiag proton-neutron
     ** matrix elements 
     **     <SD(Z)g',p';SD(N)g',n'|OP|SD(Z)g,p;SD(N)g,n>
     ** Next, calculate nondiagonal proton-neutron 
     ** matrix elements
     **     <SD(Z)g',p';SD(N)g',n'|OP|SD(Z)g,p;SD(N)g,n> 
     ** The remaining may be storen on file
     **         not yet implemented
     */

  calcCM = NO;
  pnStoreSD_MatrElem(VEFF_CALC,calcCM,&grBas,group);

  // data preparation for lanczo process

  sprintf(lancProc.title,"%s",model->title);

  lancProc.pnLanczoTotFileSize = model->file_lanczo;
  lancProc.max_iterations      = model->max_iterations;
  lancProc.states              = model->states;
  lancProc.totNumSD_ZN_oscLim  = grBas.totNumSD_ZN_oscLim;

  // run Lanczo process


  pnEigenvalueLanczoSearch(&grBas,group,&lancProc);

  if(Rank == MASTER) {
     pnLancResult(&grBas,&lancProc);
  }

  } // End function pn_LancIterateCalc()

   /*
   ** The function
   **       shellModelBasis()
   ** calculates 
   **   1. The complete set of j-scheme 
   **      and m-scheme single-particle states  
   **   2. The complete many-particle slater 
   **      determinants in m-scheme
   **      |SD(M,P)> = |SD(Z:M_Z,P_Z), SD(N:M_N,P_N)
   **      separated into groups for fixed M_Z and M_P
   **   3. Calculate the data for the MPI block structure
   **      separated into groups for fixed M_Z and M_P
   */

static void  shellModelBasis(SHELL *model,GR_BAS *grBas,
                                           GROUP **group)
{
  char     filename[ONE_LINE];
  FILE     *filePtr;

     /*
     ** Store in grBas->spBas[] all single-particle
     ** input data from SHELL model
     */

  sprintf(grBas->title,"%s",model->title);

  // Total group basis data 

  sprintf(FILE_OUTPUT,"%s%s",grBas->title,RESULT_OUTPUT);

  grBas->typeProc            = model->typeProc;
  grBas->typeVeff            = model->typeVeff;
  grBas->P                   = model->P;
  grBas->MJ                  = model->MJ;
  grBas->calc_CM             = model->calc_CM;
  grBas->spBas[0]->part      = model->part[0];
  grBas->spBas[0]->numj_orb  = model->numj_orb[0];
  grBas->spBas[0]->jbas      = model->jbas[0];
  grBas->N_low[0]            = model->N_low[0]; 
  grBas->N_high[0]           = model->N_high[0]; 
  grBas->N_low[1]            = model->N_low[1]; 
  grBas->N_high[1]           = model->N_high[1];
  grBas->delN_Z              = model->delN_Z;
  grBas->delN_N              = model->delN_N;
  grBas->delN                = model->delN; 
  grBas->oscEnergy           = model->oscEnergy;
  grBas->spBas[1]->part      = model->part[1];
  grBas->spBas[1]->numj_orb  = model->numj_orb[1];
  grBas->spBas[1]->jbas      = model->jbas[1];

    /*
    ** the complete single-particle basis 
    ** in j- and m-scheme for both protons 
    ** and neutrons
    */ 

  pn_single_particle_basis(grBas);

    /*
    ** the complete set of many-particle basis states
    **    |SD(M,P)> = |SD(Z:M_Z,P_Z), SD(N:M_N,P
    ** separated into groups for fixed M_Z and M_P
    */

  if(pn_slater_determinant(grBas,group) == FALSE)  {
    printf("Rank%d: Error in function shellModelBasis():",Rank);
    printf("\nNumber of basis |SD(Z),SD(N)> states are ZERO!!!\n\n");
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }
 
  // PrintOut of total number of basis Slater determinant

  if(Rank == MASTER) {

    sprintf(filename,"%s%s",grBas->title, RESULT_OUTPUT);

    if((filePtr = fopen(filename,"a")) == NULL) {
      printf("\n\nRank%d: Error in function pn_shell_model_calc()",
                                                             Rank);
      printf("\nWrong file = %s for output data\n\n", filename);
      MPI_Abort(MPI_COMM_WORLD,Rank);
    } 

    fprintf(filePtr,"\n\nNumber of proton |SD_Z> configurations = %d", 
	                                        grBas->totNumSD_Z); 
    fprintf(filePtr,"\nNumber of neutron |SD_N> configurations = %d",
                                               grBas->totNumSD_N); 
    fprintf(filePtr,"\nTotal number of   |SD_ZN> configurations = %d",
                                                     grBas->totNumSD_ZN); 
    fprintf(filePtr,"\n and with oscillator restriction         = %d",
                                              grBas->totNumSD_ZN_oscLim); 
    fclose(filePtr);

    // print Group structure  of basis states 

    // pn_print_group_data(filename, grBas, *group);
    
      // block structure of lanczo vectors
 
    blockSD_struct(grBas);

  } // end Rank = MASTER process

}  // End: function  shellModelBasis()

   /*
   ** The function
   **     blockSD_struct()
   ** calculates and return in SD_BLOCK sdBlock 
   ** number of |SD_init> for each MPI node
   */

void blockSD_struct(GR_BAS *grBas)
{
 static char   filename[ONE_LINE];
 int           k,initSize, blockSize,remainSD,
               startSDnum,endSDnum;
  FILE         *file_ptr;  

  // Test: NumProcs < grBas->tot_dimSD 

  if((int)grBas->totNumSD_ZN_oscLim < NumProcs) {
    sprintf(filename,"%s%s",grBas->title,RESULT_OUTPUT);
    if((file_ptr = fopen(filename,"a"))== NULL) {
      printf("\n\nRank(%d): Error in function  blockSD_struct():",
                                                            Rank);
      printf("\nWrong file = %s for the output data\n", filename);
      MPI_Abort(MPI_COMM_WORLD,Rank);
    }
    fprintf(file_ptr,"\n\nRank%d: MPI - error in function blockSD_struct():",
                                                                       Rank);
    fprintf(file_ptr,"\nNumber of processes = %d", NumProcs);
    fprintf(file_ptr," exceeds number of basis states = %d\n",
                                   grBas->totNumSD_ZN_oscLim);
    fclose(file_ptr);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }

  initSize  = grBas->totNumSD_ZN_oscLim / NumProcs;
  remainSD  = grBas->totNumSD_ZN_oscLim % NumProcs;

  startSDnum = 0; // initialization
  endSDnum   = 0;

  for(k = 0; k < NumProcs; k++) {
    blockSize = initSize + ((remainSD > 0) ? +1 : 0);
    remainSD--;
    grBas->sdBlock[k].startSDnum = startSDnum;
    endSDnum                    += blockSize; 
    grBas->sdBlock[k].endSDnum   = endSDnum;
    grBas->sdBlock[k].totNumSD   = blockSize;
    startSDnum                  += blockSize;
  }



/**************  OUTPUT TEST DATA ONLY  ***************

  //******* Block data output  for MASTER only 

  if(Rank == MASTER ) {
    sprintf(filename,"%s%s",grBas->title,RESULT_OUTPUT);
   
    if((file_ptr = fopen(filename,"a"))== NULL) {
      printf("\n\nRank%d: Error in function id_lanczo_calc.c():",Rank);
      printf("\nWrong file = %s for the output data\n", filename);
      MPI_Abort(MPI_COMM_WORLD,Rank);
    }
    fprintf(file_ptr, "\n\nBlock data for MASTER node with Rank = %d",MASTER);
    fprintf(file_ptr, "\n\nTotal processes = %d",NumProcs);
    fprintf(file_ptr, "\nsdBlock[%d].startSDnum = %d",
                         k,grBas->sdBlock[k].startSDnum);
    fprintf(file_ptr, "\nsdBlock[%d].endSDnum   = %d",
                          k,grBas->sdBlock[k].endSDnum);
    fprintf(file_ptr, "\nsdBlock[%d].totNumSD  = %d",
                         k,grBas->sdBlock[k].totNumSD);
  
    fprintf(file_ptr,"\n");
    fflush(file_ptr);
    fclose(file_ptr);
    
  } //********* end block output for MASTER only

  ************* END OUTPUT TEST DATA ONLY ***********/

  return;

} // End: function blockSD_struct()

   /*
   ** The function 
   **         print group_data()
   ** writes to output file all calculated group data.
   */

static void pn_print_group_data(char *file_name, GR_BAS  *gr_basis, GROUP  *group)
{
  int        k; 
  FILE       *file_ptr;

  if( (file_ptr= fopen(file_name,"a")) == NULL) { /* open out_data file */
    printf("\n\nRank%d: Error in function pn_print_group_data()\n\n",Rank);
     MPI_Abort(MPI_COMM_WORLD,Rank);
  }
                     /******   Group heading  ****/

  fprintf(file_ptr,"\n\n     ******* Group data *******");
  fprintf(file_ptr,"\n  k  parZ  parN  2*M_Z  2*M_N   numSD_Z   num_SD_N      tot_SD");
  fprintf(file_ptr,"    sect");

  for(k = 0; k < gr_basis->num_gr; k++)  {
    fprintf(file_ptr,"\n%3d    %c    %c    %4d   %4d    %6d     %6d    %8d",
	                             k, ((group[k].par[0] == +1) ? '+' : '-'),
                                          ((group[k].par[1] == +1) ? '+' : '-'),
            group[k].m[0], group[k].m[1], group[k].numSD[0],group[k].numSD[1],
	    group[k].numSD_ZN_oscLim);
  } 

  fprintf(file_ptr,"\n");
  fclose(file_ptr);

} // End: function pn_print_group_data()

     /*
     ** The function
     **       addID_1to2SinglePartTermVeff()
     ** adds contributions from single-particle terms to the
     ** effective two-particle matrix elements  <k.l |VEFF| k.l>
     ** and if(calc_CM == 1) <k.l |CM_INT| k.l>
     */

static void addID_1to2SinglePartTermVeff(int num_part, SP_BAS *spBas,
				       int calc_CM, MATR_ID *op_int)
{
   int         k, l, limit;
   double      factor1to2,*id_diag,*id_CM_diag;
 
   calc_CM = 0;

  id_diag = op_int->id_diag;
  if(calc_CM == 1) {
    id_CM_diag = op_int->id_CM_diag;
  }
  else {
    id_CM_diag  = (double *)NULL_PTR;
  }

  limit   = spBas->numm_orb - 1;  /* max. k - values  */
  factor1to2 = 1.0 /((double) (num_part - 1));

  for(k = 0; k < limit; k++) { 
    for(l = k+1; l <= limit; l++) {
      *(id_diag++) += (spBas->mbas[k].e + spBas->mbas[l].e)*factor1to2;

      if(calc_CM == 1) {
      *(id_CM_diag++) -= (double)(
                  2*spBas->mbas[k].osc + spBas->mbas[k].l + 1.5
		+ 2*spBas->mbas[l].osc + spBas->mbas[l].l + 1.5)*factor1to2;
      }
    } // end l-loop
  } // end k-loop

} // End: function addID_1to2SinglePartTermVeff()
