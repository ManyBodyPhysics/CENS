/*******************  The module PAR-lanc-id-process.c  ******************/

#include "PAR-shell-model.h"

     /*
     ** The entrance function
     **         id_eigenvalue_Lanczos_search()
     ** performs a Lanczos iteration for a system of identical particle to a
     **  convergence criterium is reached or a maximum iteration limit. A 
     ** temporary file is needed to store the basis Lanczos vectors. The
     **  corresponding energy matrix is diagonalized.
     ** The result is written to output file and eigenvectors in a |SD()>
     ** basis are stored iteration->eigen_vec.
     */

  /**** local function declarations ****/

void localLanczoMem(int type, int dimSD, LANC_PROC *lancProc);
     /*
     ** reserves dynamic memory if(type == 0)
     ** to be used in the Lanczos iteration 
     ** process to obtain shell model eigenvalues
     ** and the corresponding eigenvectors
     ** if(type != 0) release previous 
     ** reserved memory
     */

static void initialization_process(char *type_calc, SD_BAS *sdBas, double *initVec);
     /*
     ** returns a start vector for a Lanczos iteration 
     ** process in vec[].
     */

static void random_start_vector(SD_BAS *sdBas, double *init_vec);
     /*
     ** calculates and return a normalized start vector where
     **  amplitudes are selected through a random procedure
     */

static void storeLancVec(int lancVecNum, int dim, 
		   LANC_PROC *lancProc, double *vec);
     /*
     ** opens a sect file and stores a lanczos iteration
     ** vector,If current sect file has reacded 
     ** MAX_FILE_SIZE current sect file is closed and 
     ** a new is opened
     */

static void appendLancVecToFile(char *file_name, int n, int dim, 
                                                double *vector);
     /*
     ** opens the file "file_name" and store 
     **     int      n             - local vector number on current file
     **     int      dim           - dim of vector[]
     **     double   vector[]      - list of vector elementsa
     ** Note that the first vector on file is numbered n = ZERO and the
     ** vector is transferred to type float and stored on file
     */

static void lancIterationProc(SD_BLOCK *sdBlock, SP_BAS *spBas, 
                            SD_BAS *sdBas,LANC_PROC *lancProc);
     /*
     ** takes a start lanczos vector stored in |vec1> and 
     ** performs a Lanczos iteration process 
     ** to convergence is obtained
     */

static int orthogonalizationProc(SD_BLOCK sdBlock, double diagElem,
            LANC_PROC *lancProc, double *initVec, double *finalVec, int n);
     /*
    ** runs in parallel over all processors, orthogonalize
    ** |finalVec> against |initVec> and all previous Lanczos 
    ** vectors stored on file. 
    ** Each process treats sdBas[Rank.tot_numSD amplitudes of 
    ** |lancVec> stored locally from 
    ** lancVec> ---> lancVec + sbBlock[Rank} + sbBlock[Rank].startSDnum
    ** Finally, |lancVec> is normalized.
    ** if(norm > matrix_limit) the function return TRUE. Otherwise FALSE
    ** Finally, MASTER calculates and store <n|H|n+1> 
    */

static void readLancVec(int lancVecNum, int dim, 
		 LANC_PROC *lancProc, double *vec);
     /*
     ** opens a sect file and reads a lanczos iteration
     ** vector,If current sect file has reacded 
     ** MAX_FILE_SIZE current sect file is closed and 
     ** a new is opened
     */

static void readLancVecFromOpenFile(char *fileName, FILE **filePtr, int fileNum, 
			                   int dim, float *floatMem);
     /*
     ** reads from open file the data
     **     int      n             - local vector number on current file
     **     int      dim           - dim of vector[]
     **     float   floatMem[]     - list of vector elements
     ** Note first vector on file is numbered n = ZERO
     */

static int eigen_value_convergence(LANC_PROC *lancProc);
     /*
     ** calculates and stores the changes in eigenvalues of the model.states lowest
     ** eigenstates. The changes for the last four iterations are saved. If any 
     ** changes are larger than accuracy the function returns FALSE,
     ** otherwise TRUE.
     */

static void saveFinalData(SD_BLOCK sdBlock, SHELL *model, SP_BAS *spBas,
                                     SD_BAS *sdBas,LANC_PROC *lancProc);
     /*
     ** diagonalizes the energy matrix in Lanczos basis and
     ** transforms eigenvectors into a |SD> basis
     ** Save Eigenvalue, ang-mom and eigenvector on file
     */


static void lancToSD_Bas(LANC_PROC *lancProc, int totDimSD, double *initVec, 
                                     double *finalVec, int dimH, double *lancAmp);
     /*
     ** The function                                        
     **          lanczoToSD_Basis()                       
     ** transform an eigenvector from Lanczos vector basis 
     ** into SD basis and store the result in iteration->final_vec.   
     */

static double scalar_product_of_id_vectors(SD_BLOCK sdBlock, 
                                 double *vec1, double *vec2);
     /*
     ** calculates and returns the scalar product of 
     ** two vectors stored vec1 and vec2.
     ** Possible time reversal symmetry is included
     */

static double SquareNorm_of_a_vector(SD_BLOCK sdBlock,double *vec);
     /*
     ** calculates the squareNorm of a vector stored in 
     ** vec[] with iteration->tot_dimSD elements. 
     ** The function returns 
     **     result = sqrt( Sum(x_i * x_i))
     ** Possible time-reversal symmetry is included   
     */

static void CM_energy_matr_elem(SD_BLOCK sdBlock,LANC_PROC *lancProc,
                           SP_BAS *spBas,SD_BAS *sdBas,double *init_vec, 
                                   double *final_vec, int stored_lanc_vec);
     /*
     ** calculates nondiagonal center_of_Mass matrix 
     ** elements 
     **      <lanc_vec(k)|CM_INTfinal_vec(n)>
     ** for k = 0, 1,,...., < n
     */

static void  angularMomentumInteraction(SD_BLOCK sdBlock, SHELL *model,
                    SP_BAS *spBas,SD_BAS *sdBas, SD_STORE *sdStore);
   /*
   ** calculates and stores on file the complete set of m-scheme
   ** one- and two-particle matrix elements  of the angular
   ** momentum operator J**2. Calculate and store diagonal 
   ** and nondiagonal slater determinant matrix elements 
   ** <SD'|J**2|SD> up to the limits specified in the 
   ** input data.
   */
   
static void add_id_single_particle_ang(int num_part, SP_BAS *spBas, MATR_OP *op_ang);
   /*
   ** adds contributions from single-particle terms of J**2 to 
   ** the two-particle matrix elements  <k.l | J**2 | k.l> 
   */

 static void appendEigenVecToFile(LANC_PROC *lancProc,
                              EIGEN_DATA *h_eigenData,
			 double *j_occ,double *vector);
   /*
   ** opens the file "file_name" and store 
   **    h_eigenData.vecNum   - local vector number
   **                           on current file
   **               .dim      - dim of vector[]
   **               .eigenVal - eigenvalue amplitudes
   **               .angMom   - < |J**2|>
   **               .val_CM   - < |CM| > if val_CM = 1
   **               .numj_occ - num j-orbitals
   **
   **    j_occ[]  - list of single-particle occ prob
   **    vector[] - eigenvector ampl in float
   ** Note that the first vector on file is numbered 
   ** vecNum = ZERO and the amplitudes  is transferred 
   ** to type float and stored on file
   */

static void id_j_occupation(SD_BLOCK sdBlock, SP_BAS *spBas, SD_BAS *sdBas, 
			    double *j_occ, double *eigenVec);
   /*
   ** calculate the j-occupation numbers  of an identicle particle shell model
   ** eigenstate from the lanczo procedure and returned the result in j_occ[]. 
   */

static void writeTimeControl(LANC_PROC *lancProc,
                             TID wallTime, TID cpuTime,
                             TID wallTimeA, TID cpuTimeA);
   /*
   ** collects time controll data from
   ** all processes to MASTER an writes
   ** the result to output file 
   */

               /**** End: function declarations ****/

               /**** The function definitions  ****/ 
     /*
     ** The entrance function
     **         id_eigenvalue_Lanczos_search()
     ** performs a Lanczos iteration for a system of identical particle 
     ** over a set MPI connnected processes to a  convergence criterium 
     ** or a maximum iteration limit is reached. 
     ** The corresponding energy matrix is diagonalized in all nodes.
     ** The result is written to output file and eigenvectors in a |SD()>
     ** basis are stored throughout the MPI processes
     */

void id_eigenvalue_Lanczos_search(SD_BLOCK *sdBlock, SHELL *model,
               LANC_PROC *lancProc, SP_BAS *spBas, SD_BAS *sdBas)
{
  TID    wallTime, cpuTime;

  // reserve local memory for Lanczos process

   localLanczoMem(ALLOCATE, sdBas->tot_dimSD, lancProc);

  if(Rank == MASTER) {

    // Calculate a complete |lancStartVec>

    initialization_process(lancProc->type_calc, sdBas,lancProc->vec1);

  } // end  Rank == MASTER

  // MASTER sends |startVec> to all processes

  MPI_Bcast(lancProc->vec1, sdBas->tot_dimSD, MPI_DOUBLE, 
                                   MASTER, MPI_COMM_WORLD);
     /* 
     ** all processes store its blockpart of lancVec> on a local file
     ** |init_lancVec> = lancProc->vec1 calculated - append to file 
     */

  lancProc->storedLancVec = 0; // initialization
  storeLancVec(lancProc->storedLancVec, sdBlock[Rank].tot_numSD,
         lancProc,lancProc->vec1 + sdBlock[Rank].startSDnum);


  lancProc->storedLancVec++;


 // wall time control for total lanczos iteration 

  wallClock(4, 0); // initialization
  cpuClock(4, 0);

  wallClock(4, 1);  // start clock
  cpuClock(4, 1);

  lancIterationProc(sdBlock,spBas,sdBas,lancProc);


  wallClock(4, 2);  // stop clock
  cpuClock(4, 2);

  /*****************   print out all lanczo processes  ***********/

  if(Rank == MASTER) {
    FILE   *filePtr;

    if((filePtr = fopen(lancProc->title,"a")) == NULL)   {
      printf("\n\nRank%d: Error in function id_eigenvalue_Lanczos_search():",
                                                                        Rank);
      printf("\nWrong file = %s to open the output data file\n",lancProc->title);
       MPI_Abort(MPI_COMM_WORLD,Rank);
    }
    wallTime = wallClock(4,3);   // wall time for total lanczo process
    cpuTime  = cpuClock(4,3);   // CPU time  for total lanczo process 

    fprintf(filePtr,"\nRank%d: Total time used for Lanczos iteration proc",
                                                                 Rank);
    fprintf(filePtr,"\nWall time: %llu hour %llu min %llu sec",
                         wallTime.hour, wallTime.min, wallTime.sec);
    fprintf(filePtr,"\nCPU  time: %llu hour %llu min %llu sec\n",
                            cpuTime.hour, cpuTime.min, cpuTime.sec);
    fflush(stdout);
  
    fclose(filePtr);

  } // end print

/**********   end print out for all lanczo process  ********/

  // calculate and save final data


  saveFinalData(sdBlock[Rank],model, spBas, sdBas, lancProc);



  // release local memory for lanczo process

  localLanczoMem(FREE, sdBas->tot_dimSD, lancProc);

} // End: function  id_eigenvalue_Lanczos_search()

     /* 
     ** The function 
     **        localLanczoMem()
     ** reserves dynamic memory if(type == 0)
     ** to be used in the Lanczos iteration 
     ** process to obtain shell model eigenvalues
     ** and the corresponding eigenvectors
     ** if(type != 0) release previous 
     ** reserved memory
     */

void localLanczoMem(int type, int dimSD, LANC_PROC *lancProc) 
{
  char   *func = {"localLanczoMem(): "};  

  if(type == 0) {

    lancProc->vec1 = MALLOC(dimSD, double, func,"lancProc->vec1[]"); 
    lancProc->vec2 = MALLOC(dimSD, double, func,"lancProc->vec2[]"); 

    // local memory for the Lanczos energy matrix

    lancProc->h_matrix 
      = CALLOC(((lancProc->max_iterations + 1) 
            *(lancProc->max_iterations + 2))/2, 
             double, func, "h_matrix[]");
    
    // memory for center-of-Mass  matrix

    if(lancProc->calc_CM == 0) {
      lancProc->CM_matrix = NULL;
    }
    else {
      lancProc->CM_matrix 
	= CALLOC(( (lancProc->max_iterations + 1) 
		   *(lancProc->max_iterations + 2))/2,
		         double, func, "CM_matrix[]");
    } 

    if(Rank == MASTER) {

      // local memory to save dynamical changes in eigenvalues

      lancProc->delta_eigen = MALLOC(5 * lancProc->states, double, 
                                         func, "delta_eigen[]");
    } // end Master node
  } // end reserve memory type = 0

  else {  // type = 1

    if(Rank == MASTER) {
      free(lancProc->delta_eigen); // remove lanczo process memory 
      lancProc->delta_eigen = NULL;
    } // end MASTER node 

    free(lancProc->h_matrix);
    lancProc->h_matrix= NULL;

    if(lancProc->calc_CM == 1) {
      free(lancProc->CM_matrix);
      lancProc->CM_matrix = NULL;
    }
    free(lancProc->vec2);
    lancProc->vec2 = NULL;
    free(lancProc->vec1);
    lancProc->vec1 = NULL;
  }

} // End function localLanczoMem()

      /*
      ** The function 
      **   initialization_process();
      ** returns a start vector for a Lanczos iteration 
      ** process in vec[].
      */

static void initialization_process(char *type_calc, SD_BAS *sdBas, double *initVec)
{
  switch(strlen(type_calc)) {
      case 12 :                                   /* random-start */
               random_start_vector(sdBas, initVec);
               break;
      case 15:                                   /* random-continue */
               printf("\n\nError in function initialization_process():");
               printf("\ntype_calc = %s", type_calc);
               printf("\n Not implemented!!!!!!!!!\n");
                MPI_Abort(MPI_COMM_WORLD,Rank);
               break;
      case 16 :                                   /* fixed-J-continue */
               printf("\n\nError in function initialization_process():");
               printf("\ntype_calc = %s\n",type_calc);
               printf("\n Not implemented!!!!!!!!!\n");
                MPI_Abort(MPI_COMM_WORLD,Rank);
               break;
      default :
               printf("\n\nError in function initialization_process():");
               printf("\nWrong type_calc = %s",type_calc);
               printf("\n Not implemented!!!!!!!!!\n");
                MPI_Abort(MPI_COMM_WORLD,Rank);
               break;
  } /* end switch() */

} /* End: function initialization_process() */

      /*
      ** The function                           
      **        random_start_vector()                
      ** calculates and return a normalized start vector in vec[] 
      ** where amplitudes are selected through a random procedure
      */

static void random_start_vector(SD_BAS *sdBas, double *initVec)
{
  int
            phase, loop, dim;
  double
           *ptrVec, value, norm;

  dim = sdBas->tot_dimSD;

/***********************
  srand((unsigned long) time(NULL)); 
  phase = (rand() & TRUE) ? +1 : -1;
********************/

  phase = 1.0;

  ptrVec = initVec;
  norm = D_ZERO;

  if(sdBas->MJ == 0) { // case: time-reversal symmetry
 
    for(loop = 0; loop < sdBas->numSD[0]; loop++) {
      *ptrVec = (value = 1.0+((double)loop)/dim,value*phase);
      ptrVec++;
      norm +=  value * value;

/*********************
      phase = (rand() & TRUE) ? +1 : -1;
***********************/
 
      phase = -phase;

    }
    norm *= 2.0;     /* contr from time-reversal comp */

    for( ; loop < dim; loop++) {
      *ptrVec = (value = 1.0+((double)loop)/dim,value*phase);
      ptrVec++;
      norm +=  value * value;

/*********************
       phase = (rand() & TRUE) ? +1 : -1;
***********************/

      phase = - phase;

    }
    // normalization of lanczo vector

    if(fabs(norm) < ZERO_LIMIT)  {
      printf("\n\nError in function  random_start_vector():");
      printf("\nNorm of the initial random lanczos vector is below");
      printf(" accepted limit: norm = %E",norm);
       MPI_Abort(MPI_COMM_WORLD,Rank);
    }
    norm = 1.0 / sqrt(norm);
    ptrVec = initVec;
    for(loop = 0; loop < dim; loop++)  {
      *ptrVec *= norm;
      ptrVec++;
    }

  }  // end case: time-reversal symmetry

  else  {                            /* case: no time-reversal symmetry */

    for(loop = 0; loop < dim; loop++) {
      *ptrVec = (value = 1.0+((double)loop)/dim,value*phase);
      ptrVec++;
      norm +=  value * value;

/**************************
      phase = (rand() & TRUE) ? +1 : -1;
************/

      phase = -phase;

    }
    // normalization of lanczo vector

    if(fabs(norm) < ZERO_LIMIT)  {
      printf("\n\nError in function  random_start_vector():");
      printf("\nNorm of the initial random lanczos vector is below");
      printf(" accepted limit: norm = %E",norm);
       MPI_Abort(MPI_COMM_WORLD,Rank);
    }
    norm = 1.0 / sqrt(norm);

    ptrVec  = initVec;
    for(loop = 0; loop < dim; loop++)  {
      *ptrVec *= norm;
      ptrVec++;
    }
  } // end case: no time-reversal symmetry

} // End: function random_start_vector()


    /*
    ** The function 
    **     storeLancVec()
    ** opens a sect file and stores a lanczos iteration
    ** vector, If current sect file has reacded 
    ** MAX_FILE_SIZE current sect file is closed and 
    ** a new is opened
    */

static void storeLancVec(int lancVecNum, int dim, 
                         LANC_PROC *lancProc, double *vec)
{
  static char   basFile[ONE_LINE], sectFile[ONE_LINE];
  static int    maxLancVecPerFile, sectFileNo;
  FILE          *filePtr;

  if(lancVecNum == 0) { // initialization

    sprintf(basFile,"%s%s%d%s",lancProc->scratchFile,"Rank",
                                        Rank,LANCZO_VECTORS);

    // maximum lanczos vector in a single file

    maxLancVecPerFile 
          = (  (MAX_FILE_SIZE * M_BYTES)
	     / (dim * sizeof(float)+ 2*sizeof(int)));
    sectFileNo = 1;  // initialization
    sprintf(sectFile,"%s%d", basFile, sectFileNo);

    // create a new section file to store lanczo iteration vectors

    if((filePtr = fopen(sectFile,"wb")) == NULL) {
      printf("\n\nError in function storeLancVec():");
      printf("\nNot allowed to create the file %s\n", sectFile);
      MPI_Abort(MPI_COMM_WORLD,Rank);
    }
    fclose(filePtr);
  } // end lancVecNum = 0

  if((lancVecNum/maxLancVecPerFile) >= sectFileNo) {

    // create a new section file

    sectFileNo++;
    sprintf(sectFile,"%s%d", basFile, sectFileNo);
    if((filePtr = fopen(sectFile,"wb")) == NULL)  { 
      printf("\n\nRank(%d): Error in function storeLancVec():", Rank);
      printf("\nNot allowed to create the file %s\n", sectFile);
       MPI_Abort(MPI_COMM_WORLD,Rank);
    }
    fclose(filePtr);
  }
  appendLancVecToFile(sectFile, lancVecNum, dim, vec); 

} // End: function storeLancVec()

     /*
     ** The function                                           
     **            appendLancVecToFile( )                    
     ** opens the file "file_name" and store 
     **     int      n             - local vector number on current file
     **     int      dim           - dim of vector[]
     **     double   vector[]      - list of vector elementsa
     ** Note that the first vector on file is numbered n = ZERO and the
     ** vector is transferred to type float and stored on file
     */

static void appendLancVecToFile(char *filename, int n, int dim, double *vector) 
{
  char    *func = {" appendLancIterationVecToFile(): "};
  int     k;
  float   *tempMem, *finalPtr;
  double  *initPtr;
  FILE    *filePtr;

  if((filePtr = fopen(filename,"ab+")) == NULL) {
    printf("\nRank(%d): Error in function appendLancIterationVecToFile();",
                                                                   Rank);
    printf("\nWrong file = %s to store a vector no %d\n",filename, n);
     MPI_Abort(MPI_COMM_WORLD,Rank);
  }

  // vector number

  if(fwrite((const void *)&n,(size_t) sizeof(int), 1, filePtr) != 1) {
    printf("\nRank(%d): Error in function appendLancVecToFile()",Rank);
    printf("\nIn writing vector number =  %d", n);
    printf("\nto file %s\n",filename);
     MPI_Abort(MPI_COMM_WORLD,Rank);
  }   // end of if-test

  // dimension 

  if(fwrite((const void *)&dim,(size_t) sizeof(int),1,filePtr)!=1) { 
    printf("\nRank(%d): Error in function appendLancVecToFile()",Rank);
    printf("\nIn writing numSD =  %d", dim);
    printf("\nto file %s\n",filename);
     MPI_Abort(MPI_COMM_WORLD,Rank);
  }  // end of if-test

  tempMem = MALLOC(dim, float, func, "tempMem");
  initPtr = vector;
  finalPtr = tempMem;
  for(k = 0; k < dim; k++)  {
    *finalPtr = (float)(*initPtr);
    initPtr++;
    finalPtr++;
  }
  if(fwrite((const void *)tempMem,(size_t) sizeof(float), 
                              (size_t) dim, filePtr) != dim) { 
    printf("\n\nRank(%d): Error in function appendLancVecToFile():",
                                                              Rank);
    printf("\nIn writing %d diag matrix elements to file %s\n\n",
                                                  dim, filename);
     MPI_Abort(MPI_COMM_WORLD,Rank);
  }
  free(tempMem);
 
  fclose(filePtr);

} // End: function appendLancVecToFile()

      /*
      ** The function 
      **    lancIterationProc()
      ** takes a start lanczos vector stored in |vec1> and 
      ** performs a Lanczos iteration process 
      ** to convergence is obtained
      */

static void lancIterationProc(SD_BLOCK *sdBlock, SP_BAS *spBas,
                              SD_BAS *sdBas, LANC_PROC *lancProc)
{
  char      *func = {"lancIterationProc(): "}, sectFile[ONE_LINE]; 
  int       iterate, k, control,*rcounts, *displs;
  double    result, diagElem, *initVec, *finalVec,
            *finalPtr, *temp;

  // local data for the MPI_Allgatherv process

  rcounts = MALLOC(NumProcs, int, func,"rcounts[]");
  displs  = MALLOC(NumProcs, int, func,"displs[]");

  for(k = 0; k < NumProcs; k++) {
    rcounts[k] = sdBlock[k].tot_numSD;
    displs[k]  = sdBlock[k].startSDnum;
  } 
     /*  lancProc->run_code:
     **    = 0 - reached maximum lanczo iterations
     **    = 1 - Lanczos process has converged
     **    = 2 - no more Lanczos vectors
     */ 

  lancProc->run_code = 0; // initialization
  initVec            = lancProc->vec1; 
  finalVec           = lancProc->vec2;

  for(iterate = 1; iterate < lancProc->max_iterations; iterate++) {

    // time control for a single lanczo iteration

    wallClock(5,0); // initialization
    cpuClock(5,0);

    wallClock(5,1);  // start clock
    cpuClock(5,1);

    if(lancProc->calc_CM == 1) {
      result = 0.0;

/*************************  ikke ferdig  ************      
      result = id_lanczo_matrix_element(sdBlock[Rank],spBas,
                 sdBas,&lancProc->sdStoreCM,initVec, initVec);
****************************************/

      MPI_Allreduce(MPI_IN_PLACE,&result,1, MPI_DOUBLE, MPI_SUM,
                                                MPI_COMM_WORLD);

      lancProc->CM_matrix[( (lancProc->storedLancVec - 1) 
                   *(lancProc->storedLancVec + 2))/2]  = result;
    } // end lancProc->calc_CM == 1


    // initialization of new lanczos vector
  
    finalPtr = finalVec;
    for(k = 0; k < sdBas->tot_dimSD; k++) {
      *finalPtr = D_ZERO; 
      finalPtr++;
    }
         /*
         ** do
         **     H|initVec> = |finalVec>
         ** where the energy operator H specified through
         ** its matrix elements in m-scheme.
	 ** Contribution to |finalVec> is distributed over 
	 ** through all processes
         */


    // time control for H()*vec() calculation

    wallClock(6,0); // initialization
    cpuClock(6,0);

    wallClock(6,1);  // start clock
    cpuClock(6,1);

    id_matrix_vector_calc(spBas, sdBas, &lancProc->sdStoreVeff,
                                              initVec,finalVec);

    wallClock(6,2);  // stop clock
    cpuClock(6,2);
    
         /*
	 ** Add together all parts of |lancVec> and 
	 ** store the total |lancVec> in all processes
	 */

    MPI_Allreduce(MPI_IN_PLACE,finalVec,sdBas->tot_dimSD, MPI_DOUBLE,
                                          MPI_SUM,MPI_COMM_WORLD); 
  
    // all processes contribute  to diagonal energy matrix element

    diagElem = scalar_product_of_id_vectors(sdBlock[Rank],
                         initVec  + sdBlock[Rank].startSDnum, 
                         finalVec + sdBlock[Rank].startSDnum);

    // sum_of(diagElem) --> diagElem and stores in all processes

    MPI_Allreduce(MPI_IN_PLACE,&diagElem,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD); 

    lancProc->h_matrix[( (lancProc->storedLancVec - 1) 
	    *(lancProc->storedLancVec + 2))/2] = diagElem;          

         /*
         ** orthogonalize |finalVec>  against |initVec> and all previous
	 ** Lanczos vectors stored on file. 
	 ** Each process contains sdBas[Rank.tot_numSD amplitudes of 
	 ** |lancVec> stored locally from 
	 **     lancVec + sbBlock[Rank} + sbBlock[Rank].startSDnum
	 ** Finally, |lancVec> is normalized.
	 */
  
    if(!orthogonalizationProc(sdBlock[Rank], diagElem,lancProc,
                              initVec  + sdBlock[Rank].startSDnum,
                              finalVec + sdBlock[Rank].startSDnum,
                              lancProc->storedLancVec - 1))    {
      lancProc->run_code = 2;  // no more lanczo vector available
      break;
    }

    // Using MPI_Allgatherv produce a total |finalVec> in all processes

    MPI_Allgatherv(finalVec + sdBlock[Rank].startSDnum,
                   sdBlock[Rank].tot_numSD, MPI_DOUBLE,
                   finalVec, rcounts, displs,MPI_DOUBLE,
                   MPI_COMM_WORLD);   

    if(lancProc->calc_CM == 1) {

      // off-diagonal elements <lanc_vec(k)|H(CM)|final_vec>,
      //          k = 0,....,storedLancVec - 1 
      
      CM_energy_matr_elem(sdBlock[Rank], lancProc, spBas, sdBas,
                   initVec, finalVec, lancProc->storedLancVec);
    }

  
    // new lanczo vector found - append to file

    storeLancVec(lancProc->storedLancVec, sdBlock[Rank].tot_numSD,
                       lancProc,finalVec + sdBlock[Rank].startSDnum);

    lancProc->storedLancVec++;

    // time control after each Lanczos iteration

    wallClock(5,2);  // stop clock
    cpuClock(5,2);

    // Time controll for lanczo iyeration process

/****************************************
    if(   (iterate == 1) 
	|| ((((iterate + 1)/25)*25)==(iterate + 1)))  {
      writeTimeControl(lancProc,wallClock(5,3),cpuClock(5,3), 
		                wallClock(6,3),cpuClock(6,3));

    }
***********************************/

    if(iterate == 1) {
      writeTimeControl(lancProc,wallClock(5,3),cpuClock(5,3), 
		                wallClock(6,3),cpuClock(6,3));
    }


         /*
         ** When dimension of h_matrix has reached the required
	 ** eigenvalues, the energy spectrum is tested against
	 ** a given criterium.
	 ** If function returns TRUE energy convergence has been
	 ** reached. and the Lanczos process terminates
         */
  
    if((lancProc->storedLancVec - 1) >= lancProc->states) {

     if(Rank == MASTER) {
	control = eigen_value_convergence(lancProc);
      } // end MASTER
	
      MPI_Bcast(&control, 1, MPI_INT, MASTER, MPI_COMM_WORLD);

      if(control == TRUE) {
	lancProc->run_code = 1;
	break; // convergence reached
      }
    } // end test energy convergences


    temp     = initVec; // interchange pointers and continue
    initVec  = finalVec;
    finalVec = temp;

  } // end lanczo iteration loop

  lancProc->num_iterations = iterate;
  
  writeTimeControl(lancProc,wallClock(5,3),cpuClock(5,3), 
		                wallClock(6,3),cpuClock(6,3));

     /*
     ** Possible modification of the number
     ** of calculated energy eigenvectors
     **     lancProc->states
     ** and number of iterations
     **     lancProc->num_iterations
     */

  if(lancProc->run_code == 2) {
    lancProc->states = MIN(lancProc->states, 
                           lancProc->storedLancVec);
  }
  else {
    lancProc->states = MIN(lancProc->states,
                           lancProc->storedLancVec - 1);
  }

  // remove all files containing <SD|VEFF|SD>

  k = 1;  
  for( ; ;) {
    sprintf(sectFile,"%s%s%d%s%s%d",lancProc->scratchFile,
                     "Rank",Rank,SD_STORE_VEFF,NONDIAG,k);
    k++;
    if(remove(sectFile)) break;
  }

  free(displs);
  free(rcounts);

}  // End: function lancIterationProc()

    /*
    ** The function
    **              orthogonalizationProc()
    ** runs in parallel over all processes, orthogonalize
    ** |finalVec> against |initVec> and all previous Lanczos 
    ** vectors stored on file. 
    ** Each process treats sdBas[Rank].tot_numSD amplitudes of 
    ** |lancVec> stored locally from 
    ** lancVec> ---> lancVec + sbBlock[Rank} + sbBlock[Rank].startSDnum
    ** Finally, |lancVec> is normalized.
    ** if(norm > matrix_limit) the function return TRUE. Otherwise FALSE
    ** Finally, MASTER calculates and store <n|H|n+1> 
    */

static int orthogonalizationProc(SD_BLOCK sdBlock, double diagElem, 
       LANC_PROC *lancProc,double *initVec, double *finalVec, int n)
{
  int      k;
  double   factor, sqFactor;

       /*
       ** orthogonalize to |initVec>
       **   (|finalVec[]> - diag_elem * |initVec[]>)  
       **                ---> |finalVec[]>  
       ** Each process modifies its local part of |finalVec>
       */

  add_a_vector(- diagElem, sdBlock.tot_numSD,initVec,finalVec);

      /*
      ** Orthogonalize |finalVec> to all previous Lanczo vectors
      ** Again, each process modifies its part of |lancVec>
      */

  for(k = 0; k < n; k++)  {

    // read local part of lancVec(k) from file into |initVec[]>

    readLancVec(k, sdBlock.tot_numSD, lancProc, initVec);

    // calculate factor = <initVec|finalVec>

    factor = scalar_product_of_id_vectors(sdBlock, initVec,finalVec);

    // produce total normalization factor over all processes

    MPI_Allreduce(MPI_IN_PLACE,&factor, 1, MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD); 

    if(fabs(factor) < ZERO_LIMIT){

      // no contribution from current |lancVec(k)>

      continue;
    }
        /*
	** Orthogonalize |finalVec> to lancVec[k]>  
	**   |finalVec[]> <--- |finalVec> - factor * |initVec>
        ** Again, each modifies its part of |lancVec>
	*/
 
    add_a_vector(-factor, sdBlock.tot_numSD,initVec,finalVec);

  } // end of loop k

  // close last sectFiles for lancVec: k = -1

  if(k > 0) readLancVec(-1, sdBlock.tot_numSD, lancProc, initVec);

        /*
	** Normalize |finalVec[]> over all processes
	** Calculate off-diagonal matrix element
	**         <n|H|n-1>
	** The new lanc_vec must  have  norm > MATRIX_LIMIT.    
	*/

  sqFactor = SquareNorm_of_a_vector(sdBlock, finalVec);

  MPI_Allreduce(MPI_IN_PLACE,&sqFactor, 1, MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD); 

  factor = sqrt(sqFactor);

  if(factor < MATRIX_LIMIT)  return FALSE;

  lancProc->h_matrix[((n + 1) * (n + 4))/2 - 1] = factor;

  // Normalize |final_vec> over all processes

  factor = 1.0 / factor;

  scale_a_vector(sdBlock.tot_numSD, factor, finalVec);

  return TRUE;

} // End: function orthogonalizationProc()

    /*
    ** The function 
    **     readLancVec()
    ** opens a sect file and reads a lanczos iteration
    ** vector,If current sect file has reacded 
    ** MAX_FILE_SIZE current sect file is closed and 
    ** a new is opened
    */

static void readLancVec(int lancVecNum, int dim, 
                         LANC_PROC *lancProc, double *vec)
{
  static char   basFile[ONE_LINE], sectFile[ONE_LINE];
  static int    maxLancVecPerFile, sectFileNo;
  static float  *floatMem;
  static FILE   *filePtr;
  char          *func = {"readLancVec(): "};
  int           k;
  float         *initPtr;
  double        *finalPtr;       

  if(lancVecNum < 0) {
    free(floatMem);   // release local memory
    fclose(filePtr); // close lancVecfile and return
    return;
  }
  if(lancVecNum == 0) { // initialization

    sprintf(basFile,"%s%s%d%s",lancProc->scratchFile,"Rank",Rank,LANCZO_VECTORS);

    // maximum lanczos vector in a single file

    maxLancVecPerFile 
          = (  (MAX_FILE_SIZE * M_BYTES)
	     / (dim * sizeof(float)+ 2*sizeof(int)));
    sectFileNo = 1;  // initialization
    sprintf(sectFile,"%s%d", basFile, sectFileNo);

    // open first lanczo iteration sectFile

    if((filePtr = fopen(sectFile,"rb")) == NULL) {
      printf("\n\nError in function readLancVec():");
      printf("\nNot allowed to open file %s\n", sectFile);
       MPI_Abort(MPI_COMM_WORLD,Rank);
    }
    floatMem = MALLOC(dim, float, func, "tempMem");
  }
  if((lancVecNum/maxLancVecPerFile) >= sectFileNo) {

    // close current and open next lanczo iteration sectFile

    fclose(filePtr);

    sectFileNo++;  // initialization
    sprintf(sectFile,"%s%d", basFile, sectFileNo);
    if((filePtr = fopen(sectFile,"rb")) == NULL)  { 
      printf("\n\nError in function readLancVec():");
      printf("\nNot allowed to open file %s\n", sectFile);
       MPI_Abort(MPI_COMM_WORLD,Rank);
    }
  }
  readLancVecFromOpenFile(sectFile, &filePtr, lancVecNum, dim, floatMem); 

  // convert float amplitudes to double

  initPtr  = floatMem;
  finalPtr = vec;
  for(k = 0; k < dim; k++)  {
    *finalPtr = (float)(*initPtr);
    initPtr++;
    finalPtr++;
  }

} // End: function readLancVec()

     /*
     ** The function                                           
     **            readLancVecFromOpenFile( )                    
     ** reads from open file the data
     **     int      n             - local vector number on current file
     **     int      dim           - dim of vector[]
     **     float   floatMem[]     - list of vector elements
     ** Note first vector on file is numbered n = ZERO
     */

static void readLancVecFromOpenFile(char *filename, FILE **filePtr,
                    int lancVecNum, int dim, float *floatMem) 
{
  int     n;

  // read vector number
  if(fread((void *)&n,(size_t) sizeof(int), 1, *filePtr) != 1) {
    printf("\nError in function readLancVecFromOpenFile()");
    printf("\nIn reading vector number =  %d", n);
    printf("\nfrom file %s\n",filename);
     MPI_Abort(MPI_COMM_WORLD,Rank);
  }   // end of if-test

  // test vector number

  if(n !=  lancVecNum) {
    printf("\nRank(%d): Error in function readLancVecFromOpenFile()",Rank);
    printf("\nWrong vector number read");
    printf("\nread number = %d - should be lancVecNum = %d\n",
                                                n,lancVecNum);
     MPI_Abort(MPI_COMM_WORLD,Rank);
  }

  // read vector dimension
  if(fread((void *)&n,(size_t) sizeof(int), 1, *filePtr) != 1) {
    printf("\nRank(%d): Error in function readLancVecFromOpenFile()", Rank);
    printf("\nIn reading vector dimension =  %d", n);
    printf("\nfrom file %s\n",filename);
     MPI_Abort(MPI_COMM_WORLD,Rank);
  }   // end of if-test

  // test vector dimension
  if(n !=  dim) {
    printf("\nRank(%d): Error in function readLancVecFromOpenFile()",Rank);
    printf("\nWrong vector dimension read");
    printf("\nread number = %d - should be lancVecNum = %d\n",
                                                n,dim);
     MPI_Abort(MPI_COMM_WORLD,Rank);
  }

  // read vector float amplitudes
  if(fread((void *)floatMem,(size_t)sizeof(float), (size_t) dim,
                                       *filePtr) != (size_t) dim)  {
    printf("\n\nError in function readLancVecFromOpenFile():");
    printf("\nIn reading %d vector amplitudes from file %s\n",
                                                 dim, filename);
     MPI_Abort(MPI_COMM_WORLD,Rank);
  }
} // End: function  readLancVecFromOpenFile()

      /*
      ** The function                                     
      **        eigen_value_convergence()                 
      ** calculates and stores the changes in eigenvalues of the model.states lowest
      ** eigenstates. The changes for the last four iterations are saved. If any 
      ** changes are larger than accuracy the function returns FALSE,
      ** otherwise TRUE.
      */

static int eigen_value_convergence(LANC_PROC *lancProc)
{
  char         *func = {"eigen_value_convergence(): "};
  int          dim_h, num_states, loop, num;
  static int   delta_num = -1;
  double       *list_eigenvalue;

  dim_h      = lancProc->storedLancVec - 1;
  num_states  = lancProc->states;
  list_eigenvalue = MALLOC(dim_h, double, func, "list_eigenvalue[]");

  // Only eigenvalues are calculated

  eigenvalues(dim_h, lancProc->h_matrix, list_eigenvalue);

  if(delta_num < 0) { // save eigenvalues and return
    for(loop = 0; loop < num_states; loop++)  { // save eigenvalues and return
      lancProc->delta_eigen[4 *  num_states + loop] =  list_eigenvalue[loop];
    }
    delta_num++;
    free(list_eigenvalue);
    return FALSE;
  }
  // save changes in eigenvalues

  for(loop = 0, num = delta_num * num_states; loop < num_states; loop++, num++)  {
    lancProc->delta_eigen[num] 
         = fabs(lancProc->delta_eigen[  4 * num_states + loop] 
                                      - list_eigenvalue[loop]);
  }
      /*
      ** Increase delta_num as a pointer to a circular
      ** save-buffer containing  maximum eight sets.                  
      */

  delta_num++;
  delta_num = 3 & delta_num;

  for(loop = 0; loop < num_states; loop++)  { // save eigenvalues
    lancProc->delta_eigen[4 * num_states + loop] = list_eigenvalue[loop];
  }
  if((dim_h - num_states) < 4) {
    free(list_eigenvalue);
    return FALSE;
  }
        /*
	** If any of the saved energy differences are larger than 
	** energy accuracy return FALSE, else return TRUE.
	*/

  for(loop = 0;loop < 4 * num_states; loop++) {
    if(lancProc->delta_eigen[loop] > ENERGY_LIMIT)  {
      free(list_eigenvalue);
      return FALSE;
    }  
  }
  free(list_eigenvalue);
  
  return TRUE; // The process has converged

} // End: function eigen_value_convergence()

    /*
    ** The function                                           
    **              saveFinalData()        
    ** diagonalizes the energy matrix in Lanczos basis and
    ** transforms eigenvectors into a |SD> basis
    ** Save Eigenvalue, ang-mom and eigenvector on file
    */

static void saveFinalData(SD_BLOCK sdBlock,SHELL *model, SP_BAS *spBas,
                                    SD_BAS *sdBas, LANC_PROC *lancProc)
{
  char         *func = {"saveFinalData(): "}, sectFile[ONE_LINE];
  int          dim_h, loop, k;
  double       diagElem, result,*vecPtr, *initVec, *finalVec, *j_occ;
  EIGEN        *listEigen;
  EIGEN_DATA   h_eigenData;

     /*
     ** Calculates and stores in lancProc->op_ang[] the complete
     ** set of m-scheme one- and two-particle matrix elements
     ** of angular momentum operator J**2.
     ** Calculate and store on files diagonal and nondiagonal
     ** slater determinant matrix elements <SD'|J**2|SD>
     */ 

  angularMomentumInteraction(sdBlock, model, spBas, sdBas, &lancProc->sdStoreAng);

    /*
    ** Diagonalization of lancProc->h_matrix[], 
    ** dimension = dim_h. Eigenvalues and eigenvectors 
    ** are calculated and sorted
    */

  dim_h =  (lancProc->run_code == 2) 
          ? lancProc->storedLancVec 
          : (lancProc->storedLancVec - 1);

  listEigen = MALLOC(dim_h, EIGEN, func,"listEigen[]");
  for(loop = 0; loop < dim_h; loop++) {
    listEigen[loop].vector = MALLOC(dim_h, double, func, 
                                    "listEigen[].vector[]");
  }

  eigenvectors(dim_h, lancProc->h_matrix, listEigen);


  //****************** TEST *******************


  // memory for two |lancVec>

  initVec  = lancProc->vec1;
  finalVec = lancProc->vec2;

  // memory for single-particle occopation numbers

  j_occ = MALLOC(spBas->numj_orb,double, func, "j_occ[]"); 


  for(loop = 0; loop < lancProc->states; loop++) {

    // calculate eigenvector in sd-basis

    for(vecPtr = finalVec, k = 0; k < sdBas->tot_dimSD; vecPtr++, k++) {
      *vecPtr = D_ZERO;
    }
  
    // transfer eigenvec with number loop to SD basis


    lancToSD_Bas(lancProc, sdBlock.tot_numSD, 
		 initVec,finalVec + sdBlock.startSDnum,
                      dim_h, listEigen[loop].vector);
         /*
	 ** Add together all parts of |finalVec> and 
	 ** store the total |finalVec> in all processes
	 */

    MPI_Allreduce(MPI_IN_PLACE,finalVec,sdBas->tot_dimSD, MPI_DOUBLE,
                                          MPI_SUM,MPI_COMM_WORLD); 

    h_eigenData.vecNum   = loop; 
    h_eigenData.dim      = sdBas->tot_dimSD; 
    h_eigenData.numj_occ = spBas->numj_orb; 
    h_eigenData.eigenVal = listEigen[loop].value;

         /*
         ** do
         **     J**2|finalVec> = |initVecVec>
         ** where the energy operator J**2 specified through
         ** its matrix elements in m-scheme.
	 ** Contribution to |initVec> is distributed over 
	 ** through all processes
         */

    for(vecPtr = initVec, k = 0; k < sdBas->tot_dimSD; vecPtr++, k++) {
      *vecPtr = D_ZERO;
    }

    id_matrix_vector_calc(spBas, sdBas,&lancProc->sdStoreAng,
			                  finalVec, initVec);
         /*
	 ** Add together all parts of |initVec> and 
	 ** store the total |initVec> in all processes
	 */

    MPI_Allreduce(MPI_IN_PLACE,initVec,sdBas->tot_dimSD, MPI_DOUBLE,
                                          MPI_SUM,MPI_COMM_WORLD); 

    diagElem = scalar_product_of_id_vectors(sdBlock,
                         initVec  + sdBlock.startSDnum, 
                         finalVec + sdBlock.startSDnum);

    // sum_of(diagElem) --> diagElem and stores in all processes

    MPI_Allreduce(MPI_IN_PLACE,&diagElem,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD); 

    h_eigenData.angMom  = diagElem;

    // calculate list of single-particle occupation numbers

    id_j_occupation(sdBlock, spBas, sdBas, j_occ, finalVec);

       /*
       ** Add together all parts of h_eigenData.j_occ[]
       ** and store the total values  in all processes
       */

    MPI_Allreduce(MPI_IN_PLACE,j_occ,spBas->numj_orb,
		  MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD); 

    // save energy eigenvalue and  eigenvector on file in process MASTER

    if(Rank == MASTER) {
      appendEigenVecToFile(lancProc,&h_eigenData,j_occ,finalVec);
 
    }   // end loop:  All eigenvector written to file in process MASTER
  }

  // remove all files containing |lancVec>

  k = 1;  
  for( ; ;) {
    sprintf(sectFile,"%s%s%d%s%d",lancProc->scratchFile,"Rank",
                                    Rank,LANCZO_VECTORS, k++);
    if(remove(sectFile) != 0) break;
  }

  // remove all files containing <SD|ANG|SD>

  sprintf(sectFile,"%s%s%d%s%s",lancProc->scratchFile,"Rank",Rank,SD_STORE_ANG,DIAG);
  remove(sectFile);
  k = 1;  
  for(  ; ;) {
    sprintf(sectFile,"%s%s%d%s%s%d",lancProc->scratchFile,
                     "Rank",Rank,SD_STORE_ANG,NONDIAG, k++);
    if(remove(sectFile) != 0) break;
  }
 
  // remove local memory

  free(j_occ);
						   
  for(loop = dim_h - 1; loop >= 0;loop--) {
    free(listEigen[loop].vector); 
  }

  free(listEigen);


} // End: function saveFinaData()

      /*
      ** The function                                        
      **          lanczoToSD_Basis()                       
      ** transform an eigenvector from Lanczos vector basis 
      ** into SD basis and store the result in |finalVec>   
      ** Each MPI process calculates only 
      **        sdBlock.tot_numSD = totDimSD 
      ** part of |finalVec> 
      */

static void lancToSD_Bas(LANC_PROC *lancProc, int totDimSD, double *initVec, 
                               double *finalVec, int dimH, double *lancAmp)
{
  int  loop;

  // init values for eigenvector
 
  for(loop = 0; loop < totDimSD; loop++) finalVec[loop] = D_ZERO;

  // run through all |lancVec> stored on file 

  for(loop = 0; loop < dimH; loop++) {

    // read |lancVec(loop)> into |initVec>

    readLancVec(loop,totDimSD, lancProc, initVec); 

    // calculate the contribution from |initVec(loop)> 

    add_a_vector(*lancAmp, totDimSD, initVec, finalVec);
    lancAmp++;

  } // end of loop through all lanczo_vectors

  // close sectFiles for |lancVec>

  readLancVec(-1, totDimSD, lancProc, initVec);

} // End: function lancToSD_Bas()

       /*
       ** The function                                       
       **        scalar_product_of_id_vectors()                  
       ** calculates and returns the scalar product of 
       ** two vectors stored vec1 and vec2.
       ** Possible time reversal symmetry is included
       */

static double scalar_product_of_id_vectors(SD_BLOCK sdBlock, 
                                  double *vec1, double *vec2)
{
  int      loop;
  double   result;

  result = D_ZERO; // initalization

  if(sdBlock.MJ == 0) {
    for(loop = 0; loop < sdBlock.numSD[0]; loop++) {
      result += (*vec1) * (*vec2);
      vec1++;
      vec2++;
    }
    result *= 2.0; // end asym contributions
    for( ; loop < sdBlock.tot_numSD; loop++) {
      result += (*vec1) * (*vec2);
      vec1++;
      vec2++;
    } // end sym contributions
  } // end time reversal symmetry
  else {
    for(loop = 0; loop < sdBlock.tot_numSD; loop++) {
      result += (*vec1) * (*vec2);
      vec1++;
      vec2++;
    }
  }   //end no time reversal symmetry

  return result;

} /* End: function scalar_product_of_id_vectors() */

       /*
       ** The function                                       
       **          SquareNorm_of_a_vector()                   
       ** calculates the square norm of a vector stored in 
       ** vec[] with iteration->tot_dimSD elements. 
       ** The function returns 
       **     result = sqrt( Sum(x_i * x_i))
       ** Possible time-reversal symmetry is included   
       */

static double SquareNorm_of_a_vector(SD_BLOCK sdBlock, double *vec)
{
  int      loop; 
  double   *ptr_vec, result;

  result  = D_ZERO;
  ptr_vec = vec;
  if(!sdBlock.MJ) {
    for(loop = 0; loop < sdBlock.numSD[0]; loop++, ptr_vec++) {
      result += (*ptr_vec) * (*ptr_vec);
    }
    result *= 2.0;   // contr from time-reversal comp
    for( loop = sdBlock.numSD[0];  loop < sdBlock.tot_numSD;
                                         loop++, ptr_vec++) {
      result += (*ptr_vec) * (*ptr_vec);
    }
  } // end time reversal symmetry
  else {
    for(loop = 0; loop < sdBlock.tot_numSD; loop++, ptr_vec++) {
      result += (*ptr_vec) * (*ptr_vec);
    }
  }
  return result;

}  // End: function SquareNorm_of_a_vector()

 /*
  ** The function 
  **      CM_energy_matr_elem()
  ** calculates nondiagonal center_of_Mass matrix 
  ** elements 
  **      <lanc_vec(k)|CM_INTfinal_vec(n)>
  ** for k = 0, 1,,...., < n
  */

static void CM_energy_matr_elem(SD_BLOCK sdBlock, LANC_PROC *lancProc, 
                                SP_BAS *spBas, SD_BAS *sdBas,
                                 double *initVec, double *finalVec, int n)
{
  int     k;
  double  result;

  for(k = 1; k <= n  ; k++)  {

    // read lancVec(k) from file into |initVec[]>

    readLancVec(k - 1, sdBlock.tot_numSD, lancProc, initVec);
 
    result = 0.0;

/******************  ikke ferdig  ***************

    result = id_lanczo_matrix_element(sdBlock,spBas,sdBas,
                       &lancProc->sdStoreCM,initVec, finalVec);
******************************/

    lancProc->CM_matrix[((n - 1) * (n + 2))/2 + k] = result;

  } // end of loop k
  
  // close sectFiles for lancVec:

  readLancVec(-1, sdBas->tot_dimSD, lancProc, initVec);

}  // end function CM_energy_matr_elem()

       /*
       ** The function 
       **      angularMomentumInteraction()
       ** calculates and stores on file the complete set of m-scheme
       ** one- and two-particle matrix elements  of the angular
       ** momentum operator J**2. Calculate and store diagonal 
       ** and nondiagonal slater determinant matrix elements 
       ** <SD'|J**2|SD> on files 
       */

static void  angularMomentumInteraction(SD_BLOCK sdBlock,SHELL *model, 
                      SP_BAS *spBas, SD_BAS *sdBas, SD_STORE *sdStore)
{

  // two-particle angular momentum matrix elements

  id_ang_mom_interaction(spBas, &sdStore->op);

  // add singel-particle angular momentum matrix element

  add_id_single_particle_ang(sdBas->part, spBas, &sdStore->op);

  // sdStore ANG data initialization

  sprintf(sdStore->result,"%s", model->title); 
  sprintf(sdStore->scratchFile,"%s%s%d%s" , model->scratchFile,
                                      "Rank",Rank,SD_STORE_ANG); 

  sdStore->tot_file_size = model->file_nondiag;
  sdStore->MJ            = model->MJ;

  sdStore->storedSD_elem = 0;

  sdStore->diagMemPtr    = model->diagMemPtr;
  sdStore->tempMemPtr    = model->tempMemPtr;
  sdStore->tempMemSize   = model->tempMemSize;

  // calculate and store <SD'|ANG|SD> 

   sdStore->typeInt = ANG_INT;

   storeSD_MatrElem(spBas, sdBas, sdStore);

} // End: function  angularMomentumInteraction()

     /*
     ** The function
     **       add_id_single_particle_ang()
     ** adds contributions from single-particle terms of J**2 to 
     ** the two-particle matrix elements  <k.l | J**2 | k.l> 
     */

static void add_id_single_particle_ang(int num_part, SP_BAS *spBas, MATR_OP *op_ang)
{
   int         k, l, limit;
   double      num, *id_diag;

   limit  = spBas->numm_orb - 1;  /* max. k - values  */
   id_diag = op_ang->id_diag;

   num    = 1.0 /((double) (num_part - 1));
   for(k = 0; k < limit; k++) { 
      for(l = k+1; l <= limit; l++, id_diag++) {   
	*id_diag += 0.25 * (double) (spBas->mbas[k].j 
                                  * (spBas->mbas[k].j + 2)
                                       + spBas->mbas[l].j * (spBas->mbas[l].j + 2)) * num;
      } /* end l-loop */
   } /* end k-loop */

} /* End: function  add_id_single_particle_ang() */

      /*
     ** The function                                           
     **            appendEigenVecToFile( )                    
     ** opens the file "file_name" and store 
     **    h_eigenData.vecNum   - local vector number
     **                           on current file
     **               .dim      - dim of vector[]
     **               .eigenVal - eigenvalue amplitudes
     **               .angMom   - < |J**2|>
     **               .val_CM   - < |CM| > if val_CM = 1
     **               .numj_occ - num j-orbitals
     **
     **    j_occ[]  - list of single-particle occ prob
     **    vector[] - eigenvector ampl in float
     ** Note that the first vector on file is numbered 
     ** vecNum = ZERO and the amplitudes  is transferred 
     ** to type float and stored on file
     */

 static void appendEigenVecToFile(LANC_PROC *lancProc,
                              EIGEN_DATA *h_eigenData,
                          double *j_occ,double *vector) 
{
  char    *func = {"appendEigenVecToFile(): "},
          filename[ONE_LINE];
  int     k;
  float   *tempMem, *finalPtr;
  double  *initPtr;
  FILE    *filePtr;

  sprintf(filename,"%s%s%d%s%d",lancProc->scratchFile,"Rank",
	                                 Rank,EIGEN_VECTORS,1);

  if((filePtr = fopen(filename,"ab+")) == NULL) {
    printf("\nRank(%d): Error in function appendEigenVecToFile();",Rank);
    printf("\nWrong file = %s to store a vector no %d\n",
                              filename, h_eigenData->vecNum);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }

  // struct h_eigenData

  if(fwrite((const void *)h_eigenData,(size_t) sizeof(EIGEN_DATA),1,filePtr) != 1) {
    printf("\nRank(%d): Error in function appendEigenVecToFile()",Rank);
    printf("\nIn writing structure h_eigenData");
    printf("\nto file %s\n",filename);
     MPI_Abort(MPI_COMM_WORLD,Rank);
  }   // end of if-test

  // list of j_occupation data

  if(fwrite((const void *)j_occ,(size_t) sizeof(double), 
               h_eigenData->numj_occ, filePtr) != h_eigenData->numj_occ) {
    printf("\nRank(%d): Error in function appendEigenVecToFile()",Rank);
    printf("\nIn writing %d singel-particle occupation  data",
                                                  h_eigenData->numj_occ);
    printf("\nto file %s\n",filename);
     MPI_Abort(MPI_COMM_WORLD,Rank);
  }   // end of if-test


  tempMem = MALLOC(h_eigenData->dim, float, func, "tempMem");
  initPtr = vector;
  finalPtr = tempMem;
  for(k = 0; k < h_eigenData->dim;k++)  {
    *finalPtr = (float)(*initPtr);
    initPtr++;
    finalPtr++;
  }
  if(fwrite((const void *)tempMem,(size_t) sizeof(float), 
                (size_t) h_eigenData->dim, filePtr) != h_eigenData->dim) { 
    printf("\n\nRank(%d): Error in function appendEigenVecToFile():",Rank);
    printf("\nIn writing eigenvector no %d diag to file %s\n\n",
                               h_eigenData->dim, filename);
     MPI_Abort(MPI_COMM_WORLD,Rank);
  }
  free(tempMem);
 
  fclose(filePtr);

} // End: function appendEigenVecToFile()


   /*
   ** The function                            
   **        id_j_occupation()                   
   ** calculate the j-occupation numbers  of an identicle particle shell model
   ** eigenstate from the lanczo procedure and returned the result in j_occ[]. 
   */

static void id_j_occupation(SD_BLOCK sdBlock, SP_BAS *spBas, SD_BAS *sdBas, 
			                     double *j_occ, double *eigenVec)
{
  int       orb, loop = 0, part;
  ULL       pos, *bas_ptr;
  double    probability, *ampl_ptr;

  for(loop = 0; loop < spBas->numj_orb; loop++) j_occ[loop] = 0.0;

  if(sdBlock.MJ == 0) {   // time-reversal symmetry
    bas_ptr  = sdBas->SD + sdBlock.startSDnum; // initialization
    ampl_ptr = eigenVec + sdBlock.startSDnum;
    if(sdBlock.numSD[0]) {  // contribution from |asymSD()>
      for(loop = 0; loop < sdBlock.numSD[0]; ampl_ptr++, 
                                            bas_ptr++, loop++)   {
	probability = 2.0 * (*ampl_ptr) * (*ampl_ptr);
	pos    = ULL_ONE;  // orbital initialization
	orb    = 0;
	part   = 0;
	do   {
	  for( ;!((*bas_ptr) & pos); orb++, pos <<= 1); // occ. orbit found
	  j_occ[spBas->mbas[orb].orb] += probability;
	  orb++;
	  pos <<= 1;
	  part++;
	} while(part < sdBas->part);      
      } // loop through all |asymSD()>
    } // end contribution from asymmetric contribution
    if(sdBlock.numSD[1]) {              // contribution from |symSD()>
      bas_ptr  = sdBas->SD + sdBlock.startSDnum + sdBlock.numSD[0];
      ampl_ptr = eigenVec + sdBlock.startSDnum + sdBlock.numSD[0];
      for(loop = sdBlock.numSD[0]; loop < sdBlock.tot_numSD; ampl_ptr++, bas_ptr++, loop++)   {
	probability = (*ampl_ptr) * (*ampl_ptr);
	pos    = ULL_ONE;  // orbital initialization
	orb    = 0;
	part   = 0;
	do   {
	  for( ;!((*bas_ptr) & pos); orb++, pos <<= 1); // occ.orbit found
	  j_occ[spBas->mbas[orb].orb] += probability;
	  orb++;
	  pos <<= 1;
	  part++;
	} while(part < sdBas->part);      
      } // loop through all |symSD()>
    } // end contribution frm symmetric terms
  } // end time-reversal case
 
 else  {  // contribution from |nosymSD()>
    bas_ptr  = sdBas->SD + sdBlock.startSDnum; // initialization
    ampl_ptr = eigenVec + sdBlock.startSDnum;
    for(loop = 0; loop < sdBlock.numSD[0]; ampl_ptr++,
                                            bas_ptr++, loop++)   {
      probability = (*ampl_ptr) * (*ampl_ptr);
      pos    = ULL_ONE; // orbital initialization
      orb    = 0;
      part   = 0;
      do   {
	for( ;!((*bas_ptr) & pos); orb++, pos <<= 1); // occ. orbit found.
	j_occ[spBas->mbas[orb].orb] += probability;
	orb++;
	pos <<= 1;
	part++;
      } while(part < sdBas->part);      
    } // loop through all |nosymSD()>
  } // end contribution frm nonsymmetric case
} // End: function id_j_occupation()

     /*
     ** The function
     **     writeTimeControl()
     ** collects time controll data from
     ** all processes to MASTER an writes
     ** the result to output file 
     */

static void writeTimeControl(LANC_PROC *lancProc,
                        TID wallTime, TID cpuTime,
			     TID wallTimeA, TID cpuTimeA)
{   
  FILE         *filePtr;
  int          rank,data[12];


  if((filePtr = fopen(lancProc->title,"a")) == NULL)   {
    printf("\n\nRank%d: Error in function  lancIterationProc():", Rank);
    printf("\nWrong file = %s to open the output data file\n",lancProc->title);
     MPI_Abort(MPI_COMM_WORLD,Rank);
  }

  fprintf(filePtr,"\n\nRank%d:Stored no %d energy lanczo vectors",
                                 Rank,lancProc->storedLancVec);
  fprintf(filePtr,"\nWall time: %llu hour %llu min %llu sec",
                           wallTime.hour, wallTime.min, wallTime.sec);
  fprintf(filePtr,"\nCPU  time: %llu hour %llu min %llu sec",
                            cpuTime.hour,cpuTime.min, cpuTime.sec);

  // Print time control for H() * Vec()

  fprintf(filePtr,"\nRank %d Test time H()*Vec()",Rank);
  fprintf(filePtr,"\nRank%d: Wall time: %llu hour %llu min %llu sec",
                     Rank, wallTimeA.hour, wallTimeA.min, wallTimeA.sec);
  fprintf(filePtr,"\nRank%d: CPU  time: %llu hour %llu min %llu sec",
                     Rank, cpuTimeA.hour, cpuTimeA.min, cpuTimeA.sec);
  fflush(stdout);
  fclose(filePtr);

} // End: function writeTimeControl()
