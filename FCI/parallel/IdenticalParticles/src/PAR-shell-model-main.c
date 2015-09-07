  /*
  **          Version April  2006
  ** Shell model energy spectrum program based on an m-scheme
  ** basis of the slater determinants |SD()> 
  **
  **  MPI commands for mpr2
  ** mpdboot -r rsh -v # initialization of the mpi2 system
  ** Two file necessary in home directory
  **  .mpd.conf - contains a secret word
  **  mpd.hosts contains  all node-names included - 
  ** mpicc  file.c -o file  # complie and link
  ** mpiexec -n 10 ./file  #run 10 proc   
  ** mpdallexit # terminate mpi2 system
  */

#include "PAR-shell-model.h"

        // *** local function declaration ***

static void id_LancIterateCalc(SHELL *model);
    /*
    ** performs a Lanczos iteration process to obtain
    ** a limited number of shell model eigenstates.
    ** The result is stored on file
    */
 
static void  shellModelBasis(SD_BLOCK *sdBlock,SHELL *model,
                             SP_BAS *spBas, SD_BAS *sdBas);
   /*
   ** calculates 
   **   1. The complete set of j-scheme 
   **      and m-scheme single-particle states  
   **   2. The complete many-particle slater 
   **      determinants in m-scheme |SD(M)>
   **   3. Calculate the data for the MPI 
   **      block structure
   */
 
void blockSD_struct(SD_BLOCK *sdBlock,SHELL *model,SD_BAS *sdBas);
   /*
   ** calculates and return in SD_BLOCK sdBlock 
   ** number of |SD_init> for each node
   */

static void process_memory(SD_BLOCK sdBlock, SHELL *model,
                                            SD_BAS *sdBas);
   /*
   ** check computer limitations:
   **      1. lanc_sect_files < 1Gbytes
   **      2. total file storage for basic 
   **         lanczo vectors
   ** and calculate max lanczo vectors allowed in a single
   ** lanc_sec_file. Reserve permanent memory to store 
   ** diag/nondiag <SD'|VEFF|SD>
   */

static void lanczoIterationLimit(SHELL *model, SD_BAS *sdBas);
   /*
   ** may reduce maximum lanczos iteration due to
   ** limitation in lanczos file storage
   */

static void veffInteraction(SD_BLOCK sdBlock, SHELL *model, SP_BAS *spBas,
                                       SD_BAS *sdBas, SD_STORE *sdStore);

   /*
   ** calculates and stores diagonal and nondiagonal 
   ** slater determinant matrix elements of veffJ <SD'|VEFFJ|SD>
   ** according to specification in the input data.
   */

static void add1to2SinglePartTermVeff(int num_part, SP_BAS *spBas, 
                                                  MATR_OP *op_veff);
   /*
   ** adds contributions from single-particle terms to the effective
   ** two-particle matrix elements  <k.l |int| k.l>
   */

static void center_of_MassInteraction(SD_BLOCK sdBlock, SHELL *model,
                                      SP_BAS *spBas, SD_BAS *sdBas, 
                                                   SD_STORE *sdStore);
   /*
   ** reads from input all J_coupled two-particle center_of_Mass
   ** matrix element, calculate the complete set of m-scheme one-
   ** and two-particle matrix elements and store the result in op_CM.
   ** Calculate and store diagonal and nondiagonal slater determinant
   ** matrix elements of center_of_Mass <SD'|CM|SD> according to
   ** specifications in shell-model.h
   */

static void add1To2SinglePartTermCM(int num_part, SP_BAS *spBas, 
                                                 MATR_OP *op_CM);
   /*
   ** adds contributions from single-particle terms to
   ** center_of_Mass two-particle matrix elements
   **   <k.l |CM| k.l>
   */

static void add1To3SinglePartTermVeff(int num_part, SP_BAS *spBas,
				                  MATR_OP *op_veff);
   /*
   ** adds contributions from single-particle terms to the effective
   ** three-particle matrix elements  <k.l,m | VEFF | k.l.m>
   */

void  id_lanczo_data(SHELL *model, LANC_PROC *lancProc);
   /*
   ** sets up the data strucure LANC_PROC *lancProc 
   ** for the Lanczos process 
   */
           // *** End: local function declaration  ***

           //  *** The function definitions  ***

int main(int argc, char *argv[])
{
  char   input_file[ONE_LINE],
         *func = {" main(int argc, char *argv[]) :  "};
  FILE   *file_ptr;
  SHELL  model;
  TID   wallTime, cpuTime;


  MPI_Init(&argc, &argv);      
  MPI_Comm_size(MPI_COMM_WORLD,&NumProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &Rank);

  time_step(10,1);  // time controll for total shell model calc 

  // Read file names and open file containing basic input data

  if(argc >= 2) { // as program parameters
    strcpy(input_file, argv[1]); // save in_data input_file
  }
  else  { // or as separate typed in filename
    printf("\n Type in file name for basic shell model data = ");
    scanf("%s",input_file);
  } // end of data-file-name input

  // Time control for the total shell model calculation

  wallClock(1,0); // initialization 
  cpuClock(1,0);
  wallClock(1,1); // start time
  cpuClock(1,0);

  input_process(input_file, &model); // module shell-model-input.c

   /********************  test output  ********/

  if(Rank == MASTER) {
    FILE    *file_ptr;

    if((file_ptr = fopen(model.title,"a"))== NULL) {
      printf("\n\nRank(%d): Error in function main()():",Rank);
      printf("\nWrong file = %s for the output data\n", model.title);
      MPI_Abort(MPI_COMM_WORLD,Rank);
    }
    fprintf(file_ptr,"\n\n     NumProcs = %d\n\n\n",NumProcs);
    fflush(file_ptr);
    fclose(file_ptr);
  } // end MASTER
  /****************   end test outpt ********/


  strcpy(TEST_FILE_NAME, model.title);    // parameter for test output 


  // all files identifies with rank number

  if(  (!strcmp(model.type_calc,"random-start"))
     ||(!strcmp(model.type_calc,"dimension"))) {

       /*
       ** A Lanczos iteration shell model calculation
       ** based on random initial state.
       */

    id_LancIterateCalc(&model);

  }
  else {
    printf("\n\n The calculation process: %s not implemented\n\n",
                                                model.type_calc);
   MPI_Abort(MPI_COMM_WORLD,Rank);
  }

  // Stop and read total process time 

  wallClock(1,2); // stop tile
  cpuClock(1,2);
  wallTime = wallClock(1,3);  // read time
  cpuTime  = cpuClock(1,3);

  if(Rank == MASTER) {
    if((file_ptr = fopen(model.title,"a"))== NULL) {
      printf("\n\nError in function lanc-main.c():");
      printf("\nWrong file = %s for the output data\n", model.title);
     MPI_Abort(MPI_COMM_WORLD,Rank);
    }
    fprintf(file_ptr, "\n\nTotal shell model process wall time used");  
    fprintf(file_ptr, " %llu hour %llu min %llu sec",
	  wallTime.hour, wallTime.min, wallTime.sec);
    fprintf(file_ptr, "\nTotal shell model process cpu time used");  
    fprintf(file_ptr, " %llu hour %llu min %llu sec\n",
	  cpuTime.hour, cpuTime.min, cpuTime.sec);
    fclose(file_ptr);
  }  // end MASTER

  //   MPI_Finalize();

  MPI_Finalize();




  return 0; // sucsessfull termination

} // End: function main()

    /* 
    ** The function
    **     id_LancIterateCalc()
    ** performs a Lanczos iteration process to obtain
    ** a limited number of shell model eigenstates.
    ** The result is stored on file
    */
  
static void  id_LancIterateCalc(SHELL *model)
{
  char         *func = {" id_LancIterateCalc(): "};
  SP_BAS       spBas;
  SD_BAS       sdBas;
  LANC_PROC    lancProc;
  SD_BLOCK     *sdBlock;

  
  sdBlock = MALLOC(NumProcs, SD_BLOCK, func,"model->mem_ptr[]");

        /*
	** single-particle and many-particle
	**  shell model basis and block 
	** structure for possible MPI calculation
	*/

  shellModelBasis(sdBlock, model, &spBas, &sdBas);

        /*
	** limitation for lanczo iteration
	** and memory allocation for diag
	** <SD|VEFF|SD> and nondiag <SD|VEFF|SD>
	** matrix elements
	*/

  process_memory(sdBlock[Rank], model, &sdBas);

     /*
     ** If center_of_Mass two-particle matrix elements are included
     ** in the input file (model->calc_CM = 0, 1), calculate and store 
     ** in op_CM[] the complete set of m-scheme one- and two-particle
     ** matrix elements  of the center_of_Mass operator CM. 
     ** Calculate and store on files diagonal and nondiagonal slater
     ** determinant matrix elements <SD'|CM|SD>
     */ 

     /*
     ** Read from input all veff matrix elements and store 
     ** the complete set m-scheme veff in op_veff.
     ** Add single-particle energies to op_veff
     ** Calculate and store diagonal and nondiagonal slater
     ** determinant matrix elements <SD'|VEFF|SD> in memory 
     ** or on files according to specifications in 
     ** shell-model.h
     */


 switch(model->typeVeff) {
      case 1:
              /* 
	      ** Read from input all J_coupled  two-particle 
	      ** matrix element and store in op->veff the 
	      ** complete set of m-scheme two-particle matrix 
	      ** elements.
	      */

	      id_two_part_intJ(model->veffMatrElem, &spBas, 
                                  &lancProc.sdStoreVeff.op);
              /* 
	      ** add single-particle energies to the matrix
	      ** elements for the current particle number
	      */

	      add1to2SinglePartTermVeff(sdBas.part, &spBas,
		  	           &lancProc.sdStoreVeff.op);

	      /*
	      ** calculates and stores diagonal and nondiagonal 
	      ** slater determinant matrix elements of  <SD'|VEFFJ|SD>
	      ** according to specification in the input data.
	      */


	      lancProc.sdStoreVeff.typeInt = VEFF_INT; // twoPart veff

	      veffInteraction(sdBlock[Rank],model, &spBas, &sdBas,
                                            &lancProc.sdStoreVeff);
              break;
      case 2:
              /* 
	      ** Read from input all m-scheme two-particle matrix
	      ** element and and store the result in op_veff
	      */
 
              id_two_part_intM(model->veffMatrElem, &spBas, 
                                 &lancProc.sdStoreVeff.op);

	      /* 
	      ** add single-particle energies to the matrix
	      ** elements for the current particle number
	      */

	      add1to2SinglePartTermVeff(sdBas.part, &spBas,
                                 &lancProc.sdStoreVeff.op);

	      /*
	      ** calculates and stores diagonal and nondiagonal 
	      ** slater determinant matrix elements of  <SD'|VEFFM|SD>
	      ** according to specification in the input data.
	      */

	      lancProc.sdStoreVeff.typeInt = VEFF_INT; // twoPart veff

    	      veffInteraction(sdBlock[Rank],model, &spBas, &sdBas,
                                            &lancProc.sdStoreVeff);
	      break;
  case 3:

	      /* 
	      ** Read from input all m-scheme two-part and 
	      ** three-particle matrix element and store 
	      ** the result in op_veff
	      */

	      id_three_part_veffM(model, &spBas, &lancProc.sdStoreVeff.op);

	      /* 
	      ** add single-particle energies to the matrix
	      ** elements for the current particle number
	      */

	      add1To3SinglePartTermVeff(sdBas.part, &spBas, &lancProc.sdStoreVeff.op);

	      /*
	      ** calculates and stores diagonal and nondiagonal 
	      ** slater determinant matrix elements of  <SD'|VEFF3|SD>
	      ** according to specification in the input data.
	      */

	      lancProc.sdStoreVeff.typeInt = VEFF_INT3; // threePart veff

     	      veffInteraction(sdBlock[Rank],model, &spBas, &sdBas,
                                            &lancProc.sdStoreVeff);
	      break;
      default: 
              printf("\n\nError in function main(id_LancIterateCalc):");
	      printf("\nWrong Veff interation value - typeVeff = %d\n",
                                                         model->typeVeff);
	     MPI_Abort(MPI_COMM_WORLD,Rank);
  }  // end switch()

     /*
     ** Set up the data structure LANC_PROC *lancProc 
     ** for the Lanczos process
     */

  id_lanczo_data(model, &lancProc);

  // Lanczos iteration process: module  lanc-id-process.c

  id_eigenvalue_Lanczos_search(sdBlock, model, &lancProc,
                                          &spBas, &sdBas); 
  // Print final resultin MASTER process only

  if(Rank == MASTER) {

    id_lanczos_result(&spBas, &sdBas, &lancProc);

  }

} // End: function  id_LancIterateCalc()

   /*
   ** The function
   **       shellModelBasis()
   ** calculates 
   **   1. The complete set of j-scheme 
   **      and m-scheme single-particle states  
   **   2. The complete many-particle slater 
   **      determinants in m-scheme |SD(M)>
   **   3. Calculate the data for the MPI 
   **      block structure
   */

static void  shellModelBasis(SD_BLOCK *sdBlock, SHELL *model,
                                  SP_BAS *spBas, SD_BAS *sdBas)
{
  // the complete single-particle basis in j- and m-scheme

  single_particle_basis(model, spBas);

  // calculate the complete set of basis |SD> states:

  id_slater_determinant(model, spBas, sdBas);

 // block structure for each node 

  blockSD_struct(sdBlock, model, sdBas);

}  // End: function  shellModelBasis()

   /*
   ** The function
   **     blockSD_struct()
   ** calculates and return in SD_BLOCK sdBlock 
   ** number of |SD_init> for each node
   */

void blockSD_struct(SD_BLOCK *sdBlock, SHELL *model, SD_BAS *sdBas)
{
 int           k, num_block,remainSD, startSDnum,
               endSDnum, num;
  FILE         *file_ptr;  

  // Test: NumProcs < sb_bas.tot_dimSD 

  for(k = 0; k < NumProcs; k++) {
    sdBlock[k].MJ = sdBas->MJ;  // > 0 if no asymmetric  symmetry
  }
  
  if((int)sdBas->tot_dimSD < NumProcs) {
    if((file_ptr = fopen(model->title,"a"))== NULL) {
      printf("\n\nRank(%d): Error in function  blockSD_struct():",Rank);
      printf("\nWrong file = %s for the output data\n", model->title);
     MPI_Abort(MPI_COMM_WORLD,Rank);
    }
    fprintf(file_ptr,"\n\nMPI - error in function blockSD_struct():");
    fprintf(file_ptr,"\nNumber of processes = %d", NumProcs);
    fprintf(file_ptr," exceeds number of basis states = %d\n",sdBas->tot_dimSD);
    fclose(file_ptr);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }
  num_block = sdBas->tot_dimSD / NumProcs;
  remainSD  = sdBas->tot_dimSD % NumProcs;

  startSDnum = 0; // initialization
  endSDnum   = 0;

  for(k = 0; k < NumProcs; k++) {
    num = num_block + ((remainSD > 0) ? +1 : 0);
    remainSD--;
    sdBlock[k].startSDnum = startSDnum;
    endSDnum              += num; 
    sdBlock[k].endSDnum   = endSDnum;
    sdBlock[k].tot_numSD  = num;

    if(sdBlock[k].endSDnum <= sdBas->numSD[0]) {
      sdBlock[k].numSD[0] = num;
      sdBlock[k].numSD[1] = 0;
    }
    else if(sdBlock[k].startSDnum >= sdBas->numSD[0]) {
      sdBlock[k].numSD[0] = 0;
      sdBlock[k].numSD[1] = num;
    }
    else {
      sdBlock[k].numSD[0] = sdBas->numSD[0] - startSDnum;
      sdBlock[k].numSD[1] = num - sdBlock[k].numSD[0];
    }
    startSDnum += num;
  } // end loop through all blocks

   /******** Block data output  in MASTER only  ****/


  if(Rank  == MASTER) {
   
    if((file_ptr = fopen(model->title,"a"))== NULL) {
      printf("\n\nError in function id_lanczo_calc.c():");
      printf("\nWrong file = %s for the output data\n", model->title);
      MPI_Abort(MPI_COMM_WORLD,Rank);
    }

    fprintf(file_ptr, "\n\nBlock data for process with NumProcs = %d\n",NumProcs);
    fprintf(file_ptr, "\nBlock data for process with Rank = %d",Rank);
    fprintf(file_ptr, "\nsdBlock[0].startSDnum = %d",  sdBlock[0].startSDnum);
    fprintf(file_ptr, "\nsdBlock[0].endSDnum   = %d",  sdBlock[0].endSDnum);
    fprintf(file_ptr, "\nsdBlock[0].num_asym   = %d",  sdBlock[0].numSD[0]);
    fprintf(file_ptr, "\nsdBlock[0].num_sym    = %d",  sdBlock[0].numSD[1]);
    fprintf(file_ptr, "\nsdBlock[0].tot_numSD  = %d\n",sdBlock[0].tot_numSD);
    fflush(file_ptr);
    fclose(file_ptr);
  
  } //end block output  in MASTER

  return;

} // End: function blockSD_struct()

   /*
   ** The function 
   **    process_memory(model)
   ** check computer limitations:
   **      1. lanc_sect_files < 1Gbytes
   **      2. total file storage for basic 
   **         lanczo vectors
   ** and calculate max lanczo vectors allowed in a single
   ** lanc_sec_file. Reserve memory to store 
   ** diag/nondiag <SD'|VEFF|SD>
   */

static void process_memory(SD_BLOCK sdBlock, SHELL *model,
                                             SD_BAS *sdBas)
{
  char *func = {"process_memory(): "}; 

  lanczoIterationLimit(model, sdBas);

     /*
     ** Reservation of permanent memory to be used throughout
     ** the whole calculation to store:
     **     diag  <SD|VEFF|SD>
     ** Necessary minimum memory: 1 x sizeof(|lancVec>) 
     ** data units float
     */

  model->diagMemPtr = MALLOC((sdBlock.tot_numSD * sizeof(float)),
                                   char, func,"model->diagMemPtr[]");

  // necessary memory for nondiag <SD'|VEFF|SD>

  if(model->tempMemSize <= (UL)(sdBlock.tot_numSD * sizeof(STORE_F))) {
    printf("\n\nRank = %d: Error in function process_memory()", Rank);
    printf("\ntempMem = %u is too small",model->tempMemSize);
    printf("\n for nondiag  <SD'|VEFF|SD> ");
    printf("specified in shell_model.input\n");
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }
    model->tempMemPtr = MALLOC(model->tempMemSize,char, func,
			               "model->tempMemPtrA");
} // End: function: process_memory()

      /*
      ** The function 
      **       lanczosIterationLimit()
      ** may reduce maximum lanczos iteration due to
      ** limitation in lanczos file storage
      */

static void lanczoIterationLimit(SHELL *model, SD_BAS *sdBas)
{
  char   filename[ONE_LINE];
  int    vec_size, fileSizeMin;
  FILE   *file_ptr;

  vec_size    = sdBas->tot_dimSD * sizeof(float)+sizeof(int);  // in byte
  fileSizeMin = (model->max_iterations * vec_size)/M_BYTES + 1; // in Mbyte

  if(model->file_lanczo < fileSizeMin)  {

    // reduce number of lanczos iterations
  
    model->max_iterations = model->file_lanczo/(vec_size/M_BYTES) + 1;

    if((file_ptr = fopen(model->title,"a")) == NULL)   {
      printf("\n\nError in function lanczoIterationLimit():");
      printf("\nWrong file = %s to open the output data file\n",filename);
      MPI_Abort(MPI_COMM_WORLD,Rank);
    }
    fprintf(file_ptr,"\n\nReduction in number of lanczo iterations to %d\n",
	                                               model->max_iterations);
    fclose(file_ptr);
  }		     

} // End function  lanczoIterationLimits()

     /*
     ** The function 
     **        veffInteraction()
     ** calculates and stores diagonal and nondiagonal 
     ** slater determinant matrix elements of  <SD'|VEFFJ|SD>
     ** according to specification in the input data.
     */

static void veffInteraction(SD_BLOCK sdBlock, SHELL *model, SP_BAS *spBas, 
                  SD_BAS *sdBas, SD_STORE *sdStore)
{

  // sdStore VEFF data initialization



  strcpy(sdStore->result, model->title);
  strcpy(sdStore->scratchFile,model->scratchFile);

  sdStore->tot_file_size = model->file_nondiag;
  sdStore->MJ            = model->MJ;
  sdStore->storedSD_elem =  0;
  sdStore->diagMemPtr    = model->diagMemPtr;
  sdStore->tempMemPtr    = model->tempMemPtr;
  sdStore->tempMemSize   = model->tempMemSize;

  // calculate and store <SD'|VEFF|SD> 

  storeSD_MatrElem(spBas, sdBas, sdStore);

} // End: function veffInteraction()

     /*
     ** The function
     **       add1to2SinglePartTermVeff()
     ** adds contributions from single-particle terms to the
     ** effective two-particle matrix elements  <k.l |int| k.l>
     */

static void add1to2SinglePartTermVeff(int num_part, SP_BAS *spBas,
                                                   MATR_OP *op_veff)
{
   int         k, l, limit;
   double      factor1to2, *id_diag;

  id_diag = op_veff->id_diag;
  limit  = spBas->numm_orb - 1;  /* max. k - values  */
  factor1to2 = 1.0 /((double) (num_part - 1));
  for(k = 0; k < limit; k++) { 
    for(l = k+1; l <= limit; l++, id_diag++) {   
      *id_diag += (spBas->mbas[k].e + spBas->mbas[l].e) * factor1to2;
    } /* end l-loop */
  } /* end k-loop */

} // End: function add1to2SinglePartTermVeff()

     /*
     ** The function 
     **       center_of_MassInteraction()
     ** reads from input all J_coupled two-particle center_of_Mass
     ** matrix element, calculate the complete set of m-scheme one-
     ** and two-particle matrix elements and store the result in op_CM.
     ** Calculate and store diagonal and nondiagonal slater determinant
     ** matrix elements of center_of_Mass <SD'|CM|SD> according to
     ** specifications in shell-model.h
     */

static void center_of_MassInteraction(SD_BLOCK sdBlock, SHELL *model,
                                      SP_BAS *spBas, SD_BAS *sdBas, 
                                                   SD_STORE *sdStore)
{
       /* 
       ** Read from input all J_coupled  two-particle CM matrix element.
       ** Calculate and store the complete set of m-scheme two-particle
       ** matrix elements and store the result in op_CM
       ** Calculate and store diagonal and nondiagonal slater determinant
       ** matrix elements of center_of_Mass <SD'|CM|SD> according to
       ** specifications in shell-model.h
       */

  id_two_part_intJ(model->centerOfMass, spBas, &sdStore->op);

       /* 
       ** add single-particle veff matrix elements
       ** for the current particle number
       */

  add1To2SinglePartTermCM(sdBas->part, spBas, &sdStore->op);
 
  // sdStore_CM data initialization

  sprintf(sdStore->result,"%s", model->title); 
  sprintf(sdStore->scratchFile,"%s%s%d%s" , model->scratchFile,
                                        "Rank",Rank,SD_STORE_CM); 

  sdStore->tot_file_size = model->file_nondiag;
  sdStore->MJ            = model->MJ;
  sdStore->storedSD_elem =   0;

  sdStore->tempMemPtr  = model->tempMemPtr;
  sdStore->tempMemSize = model->tempMemSize;


  // calculate and store <SD'|CM|SD> on file 

  storeSD_MatrElem(spBas, sdBas, sdStore);

  model->file_nondiag -= sdStore->tot_file_size;  // remaining file size

} // End: function  center_of_MassInteraction()

     /*
     ** The function
     **       add1To2SinglePartTermCM()
     ** adds contributions from single-particle terms to
     ** center_of_Mass two-particle matrix elements
     **   <k.l |CM| k.l>
     */

static void add1To2SinglePartTermCM(int num_part, SP_BAS *spBas,
                                                 MATR_OP *op_CM)
{
   int         k, l, limit;
   double      num, *id_diag;

  id_diag = op_CM->id_diag;
  limit  = spBas->numm_orb - 1;  /* max. k - values  */
  num    = 1.0 /((double) (num_part - 1));
  
  for(k = 0; k < limit; k++) { 
    for(l = k+1; l <= limit; l++, id_diag++) {   
      *id_diag += (double)(2*spBas->mbas[k].osc + spBas->mbas[k].l + 1.5
		+ 2*spBas->mbas[l].osc + spBas->mbas[l].l + 1.5)* num;

    } /* end l-loop */
  } /* end k-loop */

} // End: function  add1To2SinglePartTermCM()

     /*
     ** The function
     **       add1To3SingleParttermVeff()
     ** adds contributions from single-particle terms to the effective
     ** three-particle matrix elements  <k.l,m | VEFF | k.l.m>
     */

static void add1To3SinglePartTermVeff(int num_part, SP_BAS *spBas,
                                                    MATR_OP *op_veff)
{
   int         k, l, m, limit;
   double      *id_diag, factor1to3;

  id_diag = op_veff->id_diag;
  limit  = spBas->numm_orb - 2;  // max. k - values
  factor1to3 = 2.0 /((double)((num_part - 1)*(num_part - 2)));
  for(k = 0; k < limit; k++) { 
    for(l = k+1; l <= limit; l++) {  
      for(m = l+1; m < spBas->numm_orb; m++) {
	*id_diag += (  spBas->mbas[k].e + spBas->mbas[l].e
                     + spBas->mbas[m].e ) * factor1to3;
        id_diag++;
      } // end m-loop
    } // end l-loop
  } // end k-loop

} // End: function  add1To3SinglePartTermVeff()

     /*
     ** The function
     **         id_lanczo_data()
     ** sets up the data strucure LANC_PROC *lancProc 
     ** for the Lanczos process 
     */

void  id_lanczo_data(SHELL *model,LANC_PROC *lancProc)
{
  lancProc->calc_CM         = model->calc_CM;

  strcpy(lancProc->title, model->title);
  strcpy(lancProc->scratchFile, model->scratchFile);
  strcpy(lancProc->type_calc, model->type_calc);
       
  lancProc->max_iterations     = model->max_iterations;
  lancProc->states             = model->states;

  lancProc->file_lanczo = model->file_lanczo; // file storage for lanc vec

} // End: function: id_lanczo_data()
