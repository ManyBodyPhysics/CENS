  /*
  **          Version october 2005
  ** Shell model energy spectrum program based on an m-scheme
  ** basis of the slater determinants |SD()> 
  */

#include "shell-model.h"

        // *** local function declaration ***

static void id_LancIterateCalc(SHELL *model);
    /*
    ** performs a Lanczos iteration process to obtain
    ** a limited number of shell model eigenstates.
    */
 
static void  shellModelBasis(SHELL *model, SP_BAS *spBas, SD_BAS *sdBas);
   /*
   ** calculates 
   **   1. The complete set of j-scheme 
   **      and m-scheme single-particle states  
   **   2. The complete many-particle slater 
   **      determinants in m-scheme |SD(M)>
   */

static void process_memory(SHELL *model, SD_BAS *sdBas);
   /*
   ** check computer limitations:
   **      1. lanc_sect_files < 1Gbytes
   **      2. total file storage for basic 
   **         lanczo vectors
   ** and calculate max lanczo vectors
   ** allowed in a single lanc_sec_file
   */

static void lanczoIterationLimit(SHELL *model, SD_BAS *sdBas);
   /*
   ** may reduce maximum lanczos iteration due to
   ** limitation in lanczos file storage
   */

static void veffInteraction(SHELL *model, SP_BAS *spBas,
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

static void center_of_MassInteraction(SHELL *model, SP_BAS *spBas,
                                  SD_BAS *sdBas, SD_STORE *sdStore);
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
   ** adds contributions from single-particle terms to center_of_Mass
   ** two-particle matrix elements  <k.l |CM| k.l>
   */

static void add1To3SinglePartTermVeff(int num_part, SP_BAS *spBas,
				                  MATR_OP *op_veff);
   /*
   ** adds contributions from single-particle terms to the effective
   ** three-particle matrix elements  <k.l,m | VEFF | k.l.m>
   */

void  id_lanczo_data(SHELL *model, SD_BAS *sdBas,
  		            LANC_PROC *lanc_proc);
   /*
   ** sets up the data strucure LANC_PROC *lanc_proc 
   ** for the Lanczos process 
   */
           // *** End: local function declaration  ***

           //  *** The function definitions  ***

int main(int argc, char *argv[])
{
  char   input_file[ONE_LINE],  filename[ONE_LINE]; 
  FILE   *file_ptr;
  SHELL  model;
  TID   ex_time;

  time_step(10,1);  // time controll for total shell model calc 

  // Read file names and open file containing basic input data

  if(argc >= 2) { // as program parameters
    strcpy(input_file, argv[1]); // save in_data filename
  }
  else  { // or as separate typed in filename
    printf("\n Type in file name for basic shell model data = ");
    scanf("%s",input_file);
  } // end of data-file-name input

    /*
    ** Necessary default dparameter are written to structure SHELL model
    ** The remaining data for the shell model problem and the lanczos
    ** procedure are read and stored in SHELL model, checked and written
    ** to output file.
    */

  model.permMemSize  = (UL)(PERM_MEM * M_BYTES); // default mem in byte
  model.windMemSize  = (UL)(WIND_MEM * M_BYTES);

  model.file_nondiag = LANC_FILES * 1000; // default files in Mbyte
  model.file_lanczo  = MATR_ELEM_FILES * 1000;

  input_process(input_file, &model); // module shell-model-input.c

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
    exit(1);
  }

  ex_time = time_step(10,2);  // read clock for total process

  sprintf(filename,"%s%s",model.title,RESULT_OUTPUT); // data to output file

  if((file_ptr = fopen(filename,"a"))== NULL) {
    printf("\n\nError in function lanc-main.c():");
    printf("\nWrong file = %s for the output data\n", filename);
    exit(1);
  }
  fprintf(file_ptr, "\n\nTotal shell model process time used");  
  fprintf(file_ptr, " %lld hour %lld min %lld sec\n",
	  ex_time.hour, ex_time.min, ex_time.sec);
  fclose(file_ptr);

  printf("\n\n");

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
  SP_BAS       spBas;
  SD_BAS       sdBas;
  LANC_PROC    lancProc;

  // single-particle and many-particle shell model basis

  shellModelBasis(model, &spBas, &sdBas);

  // Computer limitation and memory allocation

  process_memory(model, &sdBas);

     /*
     ** If center_of_Mass two-particle matrix elements are included
     ** in the input file (model->calc_CM = 0, 1), calculate and store 
     ** in op_CM[] the complete set of m-scheme one- and two-particle
     ** matrix elements  of the center_of_Mass operator CM. 
     ** Calculate and store on files diagonal and nondiagonal slater
     ** determinant matrix elements <SD'|CM|SD>
     */ 

  if(model->calc_CM == 1) {

    lancProc.sdStoreCM.typeInt = CM_INT;

    center_of_MassInteraction(model, &spBas, &sdBas, 
			           &lancProc.sdStoreCM);
  }
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
	      ** Calculate and store the complete
	      ** set of m-scheme two-particle matrix
	      ** elements and store the result in op_veff
	      */
 
	      lancProc.sdStoreVeff.typeInt = VEFF_INT; // twoPart veff

    	      veffInteraction(model, &spBas, &sdBas, &lancProc.sdStoreVeff);
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
	      ** Calculate and store the complete
	      ** set of m-scheme two-particle matrix
	      ** elements and store the result in op_veff
	      */
 
	      lancProc.sdStoreVeff.typeInt = VEFF_INT; // twoPart veff

	      veffInteraction(model, &spBas, &sdBas,
                               &lancProc.sdStoreVeff);
	      break;
      case 3: 
	      /* 
	      ** Read from input all m-scheme three-particle matrix
	      ** element and and store the result in op_veff
	      */

	      id_three_part_veffM(model, &spBas,
			      &lancProc.sdStoreVeff.op);


	      /* 
	      ** add single-particle energies to the matrix
	      ** elements for the current particle number
	      */

	      add1To3SinglePartTermVeff(sdBas.part, &spBas,
                                   &lancProc.sdStoreVeff.op);
	      /*
	      ** Calculate and store the complete
	      ** set of m-scheme two-particle matrix
	      ** elements and store the result in op_veff
	      */

	      lancProc.sdStoreVeff.typeInt = VEFF_INT3; // threePart veff

	      veffInteraction(model, &spBas, &sdBas,
                               &lancProc.sdStoreVeff);
	      break;
      default: 
              printf("\n\nError in function main(id_LancIterateCalc):");
	      printf("\nWrong Veff interation value - typeVeff = %d\n",
                                                         model->typeVeff);
	      exit(1);
  }  // end switch()

     /*
     ** Set up the data structure LANC_PROC *lancProc 
     ** for the Lanczos process
     */

  id_lanczo_data(model, &sdBas, &lancProc);

      // Lanczos iteration process: module  lanc-id-process.c

  id_eigenvalue_Lanczos_search(model, &lancProc, &spBas, &sdBas); 

  // Print final result:

  id_lanczos_result(&spBas, &sdBas, &lancProc);


} // End: function  id_LancIterateCalc()

   /*
   ** The function
   **       shellModelBasis()
   ** calculates 
   **   1. The complete set of j-scheme 
   **      and m-scheme single-particle states  
   **   2. The complete many-particle slater 
   **      determinants in m-scheme |SD(M)>
   */

static void  shellModelBasis(SHELL *model, SP_BAS *spBas, SD_BAS *sdBas)
{
  // the complete single-particle basis in j- and m-scheme

  single_particle_basis(model, spBas);

  // calculate the complete set of basis |SD> states:

  id_slater_determinant(model, spBas, sdBas);

}  // End: function  shellModelBasis()

   /*
   ** The function 
   **    process_memory(model)
   ** check computer limitations:
   **      1. lanc_sect_files < 1Gbytes
   **      2. total file storage for basic 
   **         lanczo vectors
   ** and calculate max lanczo vectors allowed in a single
   ** lanc_sec_file. Reserve permanent memory to store 
   ** diag/nondiag <SD'|VEFF|SD>
   */

static void process_memory(SHELL *model, SD_BAS *sdBas)
{
  char *func = {"process_memory(): "}; 

  lanczoIterationLimit(model, sdBas);

     /*
     ** Reservation of permanent memory to be used throughout
     ** the whole calculation to store:
     **     diag/nondiag  <SD'|VEFF|SD>
     ** Necessary minimum memory: 1 x sizeof(|lancVec>) 
     ** data units float
     */

  if(model->windMemSize <= (UL)(sdBas->tot_dimSD*sizeof(float))) {
    printf("\n\n Error in function process_memory()");
    printf("\n in File shell-modelmain.c");
    printf("\n WIND_MEM = %d is too small", model->windMemSize);
    printf("Specified in shell_model.h\n");
    exit(1);
  }
  
  if(model->permMemSize <= (UL)(sdBas->tot_dimSD*sizeof(float))) {
    model->permMemSize = (UL)(sdBas->tot_dimSD*sizeof(float));
    model->permMemPtr  = MALLOC(model->permMemSize,char, func,
                                       "model->permMemPtr[]");
    model->memCase    = 0;
    model->windMemPtr = MALLOC(model->windMemSize,char, func,
                                       "model->windMemPtr[]");
  }
  else {
    model->memCase = 1;
    model->permMemPtr = MALLOC(model->permMemSize,char, func,
                                       "model->pernMemPtr[]");
    model->windMemPtr = MALLOC(model->windMemSize,char, func,
                                       "model->windMemPtr[]");
  }
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
    sprintf(filename,"%s%s",model->title,RESULT_OUTPUT);
    if((file_ptr = fopen(filename,"a")) == NULL)   {
      printf("\n\nError in function lanczoIterationLimit():");
      printf("\nWrong file = %s to open the output data file\n",filename);
      exit(1);
    }
    fprintf(file_ptr,"\n\nReduction in number of lanczo iterations to %d\n",
	                                               model->max_iterations);
    fclose(file_ptr);
  }		     

} // End function  lanczoIterationLimits()

     /*
     ** The function 
     **        veffIntcenter_of_Meraction()
     ** calculates and stores diagonal and nondiagonal 
     ** slater determinant matrix elements of  <SD'|VEFFJ|SD>
     ** according to specification in the input data.
     */

static void veffInteraction(SHELL *model, SP_BAS *spBas, 
                             SD_BAS *sdBas, SD_STORE *sdStore)       
{

  // sdStore VEFF data initialization

  sprintf(sdStore->result,"%s%s", model->title, RESULT_OUTPUT); 
  sprintf(sdStore->title,"%s%s" , model->title, SD_STORE_VEFF); 

  sdStore->tot_file_size = model->file_nondiag;
  sdStore->diagMemCalc   = YES;
  sdStore->memCase       = model->memCase;
  sdStore->MJ            = model->MJ;
 
  sdStore->lastSD_mem_asym    = - 1;
  sdStore->lastSD_mem_sym     = - 1;  
  sdStore->storedSD_mem_elem  =   0;

  sdStore->lastSD_file_asym   = - 1;
  sdStore->lastSD_file_sym    = - 1;  
  sdStore->storedSD_file_elem =   0;

  sdStore->permMemPtr      = model->permMemPtr;
  sdStore->permMemSize     = model->permMemSize;
  sdStore->freePermMemPtr  = sdStore->permMemPtr;
  sdStore->freePermMemSize = sdStore->permMemSize;
  sdStore->windMemPtr      = model->windMemPtr;
  sdStore->windMemSize     = model->windMemSize;

  // calculate and store <SD'|VEFF|SD> 

   storeSD_MatrElemInMem(spBas, sdBas, sdStore);

   storeSD_MatrElemOnFile(spBas, sdBas, sdStore);

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

static void center_of_MassInteraction(SHELL *model, SP_BAS *spBas,
                                   SD_BAS *sdBas,SD_STORE *sdStore) 
{
       /* 
       ** Read from input all J_coupled  two-particle CM matrix element.
       ** Calculate and store the complete set of m-scheme two-particle
       ** matrix elements and store the result in op_CM
       */

  id_two_part_intJ(model->centerOfMass, spBas, &sdStore->op); 

       /* 
       ** add single-particle veff matrix elements
       ** for the current particle number
       */

  add1To2SinglePartTermCM(sdBas->part, spBas, &sdStore->op);
 
  // sdStore_CM data initialization

  sprintf(sdStore->result,"%s%s", model->title, RESULT_OUTPUT); 
  sprintf(sdStore->title,"%s%s" , model->title, SD_STORE_CM); 

  sdStore->tot_file_size = model->file_nondiag;
  sdStore->diagMemCalc   = NO;
  sdStore->memCase       = 0;
  sdStore->MJ            = model->MJ;
 
  sdStore->lastSD_mem_asym    = - 1;
  sdStore->lastSD_mem_sym     = - 1;  
  sdStore->storedSD_mem_elem  =   0;

  sdStore->lastSD_file_asym   = - 1;
  sdStore->lastSD_file_sym    = - 1;  
  sdStore->storedSD_file_elem =   0;

  sdStore->permMemPtr      = model->permMemPtr;
  sdStore->permMemSize     = model->permMemSize;
  sdStore->freePermMemPtr  = sdStore->permMemPtr;
  sdStore->freePermMemSize = sdStore->permMemSize;
  sdStore->windMemPtr      = model->windMemPtr;
  sdStore->windMemSize     = model->windMemSize;

  // calculate and store <SD'|CM|SD> on file 


  storeSD_MatrElemOnFile(spBas, sdBas, sdStore);

  model->file_nondiag -= sdStore->tot_file_size; // remaining file size

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
      *id_diag += (double)(
                  2*spBas->mbas[k].osc + spBas->mbas[k].l + 1.5
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
	*(id_diag++) += (  spBas->mbas[k].e + spBas->mbas[l].e
                     + spBas->mbas[m].e ) * factor1to3;
      } // end m-loop
    } // end l-loop
  } // end k-loop

} // End: function  add1To3SinglePartTermVeff()

     /*
     ** The function
     **         id_lanczo_data()
     ** sets up the data strucure LANC_PROC *lanc_proc 
     ** for the Lanczos process 
     */

void  id_lanczo_data(SHELL *model, SD_BAS *sdBas,
		     LANC_PROC *lanc_proc)
{
  lanc_proc->memCase         = model->memCase;
  lanc_proc->calc_CM         = model->calc_CM;

  strcpy(lanc_proc->title, model->title); 
  strcpy(lanc_proc->type_calc, model->type_calc);
       
  lanc_proc->max_iterations     = model->max_iterations;
  lanc_proc->states             = model->states;

  lanc_proc->file_lanczo = model->file_lanczo; // file storage for lanc vec

} // End: function: id_lanczo_data()
