
/************************ File: SINGLE-pnLancProcess.c *************/

#include "PAR-pnShellModel.h"

     /*
     ** The entrance function
     **         pnEigenvalueLanczoSearch()
     ** performs a Lanczos iteration for a system of identical particle 
     ** over a set MPI connnected processes to a  convergence criterium 
     ** or a maximum iteration limit is reached. 
     ** The corresponding energy matrix is diagonalized in all nodes.
     ** The result is written to output file and eigenvectors in a |SD()>
     ** basis are stored throughout the MPI processes
     */
          /**** local function declarations ****/

static void localLanczoMem(int type,LANC_PROC *lancProc);
     /* 
     ** reserves dynamic memory if(type == 0)
     ** to be used in the Lanczos iteration 
     ** process to obtain shell model eigenvalues
     ** and the corresponding eigenvectors
     ** if(type != 0) release previous 
     ** reserved memory
     */

static void random_start_vector(int totDimSD, double *initVec);
     /*
     ** calculates and return a normalized start vector in vec[] 
     ** where amplitudes are selected through a random procedure
     */

static void storeLancVec(char *title,int lancVecNum,int dim,
			                              double *vec);
    /*
    ** opens a sect file and stores a lanczos iteration
    ** vector, If current sect file has reacded 
    ** MAX_FILE_SIZE current sect file is closed and 
    ** a new is opened
    */

static void appendLancVecToFile(char *filename, int n, int dim, 
                                                double *vector);
    /*
    ** opens the file "file_name" and store 
    **     int      n             - local vector number on current file
    **     int      dim           - dim of vector[]
    **     double   vector[]      - list of vector elementsa
    ** Note that the first vector on file is numbered n = ZERO and the
    ** vector is transferred to type float and stored on file
    */

static void lancIterationProc(GR_BAS *grbas,GROUP *group,
			      LANC_PROC *lancProc);
    /*
    ** takes a start lanczos vector stored in |vec1> and 
    ** performs a Lanczos iteration process 
    ** to convergence is obtained
    */

static double scalarProductOf_pnVectors(int dim, double *vec1, 
                                                 double *vec2);
    /*
    ** calculates and returns the scalar product of 
    ** two vectors stored vec1 and vec2.
    */

static int orthogonalizationProc(int dim, double diagElem, 
	                         LANC_PROC *lancProc,
                                 double *initVec, double *finalVec, 
                                 int n);
    /*
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

static void readLancVec(int lancVecNum, int dim, 
			LANC_PROC *lancProc, double *vec);
    /*
    ** opens a sect file and reads a lanczos iteration
    ** vector,If current sect file has reacded 
    ** MAX_FILE_SIZE current sect file is closed and 
    ** a new is opened
    */

static void readLancVecFromOpenFile(char *filename, FILE **filePtr,
			   int lancVecNum, int dim, float *floatMem);
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

static double SquareNorm_of_a_vector(int dim, double *vec);
    /*
    ** calculates the square norm of a vector stored in 
    ** vec[] with iteration->tot_dimSD elements. 
    ** The function returns 
    **     result = sqrt( Sum(x_i * x_i))
    ** Possible time-reversal symmetry is included   
    */

static void saveFinalData(GR_BAS *grBas,GROUP *group, LANC_PROC *lancProc);
   /*
   ** diagonalizes the energy matrix in Lanczos basis and
   ** transforms eigenvectors into a |SD> basis
   ** Save Eigenvalue, ang-mom and eigenvector on file
   */

static void lancToSD_Bas(LANC_PROC *lancProc, int totDimSD, double *initVec, 
			 double *finalVec, int dimH, double *lancAmp);
   /*
   ** transform an eigenvector in Lanczos vector basis 
   ** into SD basis and store the result in |finalVec>   
   */

static void proton_j_occupation(GR_BAS *grBas, GROUP *group,
			 double *j_occ, double *eigenVec);
   /*
   ** calculate the j-occupation numbers  of protons for a shell model
   ** eigenstate eigenVec[]
   */
   
static void neutron_j_occupation(GR_BAS *grBas, GROUP *group,
			  double *j_occ, double *eigenVec);
   /*
   ** calculate the j-occupation numbers  of neutrons for a shell model 
   ** eigenstate eigenVec[]
   */

static void storeEigenVec(LANC_PROC *lancProc, EIGEN_DATA h_eigenData,
                          double *j_occZ,double *j_occN,double *vec);
   /*
   ** opens a sect file and stores a lanczos 
   ** eigenvector, If current sect file has 
   ** reacded MAX_FILE_SIZE current sect file 
   ** is closed and a new is opened
   */

static void appendEigenVecToFile(char *filename, EIGEN_DATA h_eigenData,
				 double *j_occZ,double *j_occN,double *vector);
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

static void writeTimeControl(LANC_PROC *lancProc,
                        TID wallTime, TID cpuTime,
                        TID wallTimeA, TID cpuTimeA,
			TID wallTimeB, TID cpuTimeB);
   /*
   ** time control data to output file 
   */

          /**** End local function declarations ****/

    /*
    ** The entrance function
    **         pnEigenvalueLanczoSearch()
    ** performs a Lanczos iteration for a system of identical particle 
    ** over a set MPI connnected processes to a  convergence criterium 
    ** or a maximum iteration limit is reached. 
    ** The corresponding energy matrix is diagonalized in all nodes.
    ** The result is written to output file and eigenvectors in a |SD()>
    ** basis are stored throughout the MPI processes
    */

void pnEigenvalueLanczoSearch(GR_BAS *grBas,GROUP *group,
                                      LANC_PROC *lancProc)
{
  TID    wallTime, cpuTime;

  // reserve local memory for Lanczos process

  localLanczoMem(ALLOCATE,lancProc);

  // Calculate a complete |lancStartVec>

  if(Rank == MASTER) {

    random_start_vector(grBas->totNumSD_ZN_oscLim,lancProc->vec1);

  } // end  Rank == MASTER

  // MASTER sends |startVec> to all processes

  MPI_Bcast(lancProc->vec1,grBas->totNumSD_ZN, MPI_DOUBLE, 
                                   MASTER, MPI_COMM_WORLD);

  lancProc->storedLancVec = 0; // initialization

  storeLancVec(lancProc->title,lancProc->storedLancVec,
	                      grBas->sdBlock[Rank].totNumSD,
	    lancProc->vec1 + grBas->sdBlock[Rank].startSDnum);

  lancProc->storedLancVec++;

  // wall time control for total lanczos iteration 

  wallClock(4, 0); // initialization
  cpuClock(4, 0);

  wallClock(4, 1);  // start clock
  cpuClock(4, 1);

  lancIterationProc(grBas,group,lancProc);

  wallClock(4, 2);  // stop clock
  cpuClock(4, 2);

  if(Rank == MASTER) { // print out all lanczo proc
    {
      char   filename[ONE_LINE];
      FILE   *filePtr;

      sprintf(filename,"%s%s",lancProc->title,RESULT_OUTPUT);
      if((filePtr = fopen(filename,"a")) == NULL)   {
	printf("\n\nRank%d: Error in function id_eigenvalue_Lanczos_search():",
                                                                          Rank);
	printf("\nWrong file = %s to open output data file\n",filename);
	MPI_Abort(MPI_COMM_WORLD,Rank);
      }
      wallTime = wallClock(4,3); // wall time for total lanczo process
      cpuTime  = cpuClock(4,3);  // CPU time  for total lanczo process 
      
      fprintf(filePtr,"\n\nTotal time used for Lanczos iteration proc");
      fprintf(filePtr,"\nWall time: %llu hour %llu min %llu sec",
	      wallTime.hour, wallTime.min, wallTime.sec);
      fprintf(filePtr,"\nCPU  time: %llu hour %llu min %llu sec\n",
                            cpuTime.hour, cpuTime.min, cpuTime.sec);
      fflush(stdout);
  
      fclose(filePtr);

    } // end print
  } // MASTER process

  // calculate and save final data

  saveFinalData(grBas,group,lancProc);

  // print out all lanczo proc

  if(Rank == MASTER) {
    {
      char   filename[ONE_LINE];
      FILE   *filePtr;

      sprintf(filename,"%s%s",lancProc->title,RESULT_OUTPUT);
      if((filePtr = fopen(filename,"a")) == NULL)   {
	printf("\n\nRank%d: Error in function id_eigenvalue_Lanczos_search():",
                                                                          Rank);
	printf("\nWrong file = %s to open the output data file\n",filename);
	MPI_Abort(MPI_COMM_WORLD,Rank);
      }
      wallTime = wallClock(4,3);   // wall time for total lanczo process
      cpuTime  = cpuClock(4,3);   // CPU time  for total lanczo process 
      
      fprintf(filePtr,"\n\nTotal time used for Lanczos iteration proc");
      fprintf(filePtr,"\nWall time: %llu hour %llu min %llu sec",
	      wallTime.hour, wallTime.min, wallTime.sec);
      fprintf(filePtr,"\nCPU  time: %llu hour %llu min %llu sec\n",
                            cpuTime.hour, cpuTime.min, cpuTime.sec);
      fflush(stdout);
  
      fclose(filePtr);

    } // end print
  } // MASTER process

  // release local memory for lanczo process

  localLanczoMem(FREE,lancProc);

} // End: function pnEigenvalueLanczosearch()

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

void localLanczoMem(int type, LANC_PROC *lancProc) 
{
  char   *func = {"localLanczoMem(): "};  

  if(type == 0) {// allocate memory

    lancProc->vec1 = MALLOC(lancProc->totNumSD_ZN_oscLim, double, 
                                         func,"lancProc->vec1[]"); 
    lancProc->vec2 = MALLOC(lancProc->totNumSD_ZN_oscLim, double, 
                                         func,"lancProc->vec2[]"); 

    // local memory for Lanczos energy matrix

    lancProc->h_matrix 
      = CALLOC( ((lancProc->max_iterations + 1) 
               *(lancProc->max_iterations + 2))/2, 
                       double, func, "h_matrix[]");
    
      // local memory to save dynamical changes in eigenvalues

    lancProc->delta_eigen 
         = MALLOC(5*lancProc->states,double,func,"delta_eigen[]");
  } // end reserve memory type = 0

  else {  // type = 1 // release memory

    free(lancProc->delta_eigen);
    lancProc->delta_eigen = NULL;


    free(lancProc->h_matrix);
    lancProc->h_matrix= NULL;
    free(lancProc->vec2);
    lancProc->vec2 = NULL;
    free(lancProc->vec1);
    lancProc->vec1 = NULL;
  }
} // End function localLanczoMem()

      /*
      ** The function                           
      **        random_start_vector()                
      ** calculates and return a normalized start vector in vec[] 
      ** where amplitudes are selected through a random procedure
      */

static void random_start_vector(int dim, double *initVec)
{
  int      phase, loop;
  double   *ptrVec, value, norm;


/***********************
  srand((unsigned long) time(NULL)); 
  phase = (rand() & TRUE) ? +1 : -1;
********************/

  phase = 1;

  ptrVec = initVec;
  norm = D_ZERO;

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
    printf("\n\nRank%d: Error in function  random_start_vector():",Rank);
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

} // End: function random_start_vector()

    /*
    ** The function 
    **     storeLancVec()
    ** opens a sect file and stores a lanczos iteration
    ** vector, If current sect file has reacded 
    ** MAX_FILE_SIZE current sect file is closed and 
    ** a new is opened
    */

static void storeLancVec(char *title,int lancVecNum,int dim, double *vec)
{
  static char   basFile[ONE_LINE], sectFile[ONE_LINE];
  static int    maxLancVecPerFile, sectFileNo;
  FILE          *filePtr;

  if(lancVecNum == 0) { // initialization

    sprintf(basFile,"%s%d%s%s","Rank",Rank,title,LANCZO_VECTORS);

    // maximum lanczos vector in a single file

    maxLancVecPerFile 
          = (  (MAX_FILE_SIZE * M_BYTES)
	     / (dim * sizeof(float)+ 2*sizeof(int)));
    sectFileNo = 1;  // initialization
    sprintf(sectFile,"%s%d", basFile, sectFileNo);

    // create a new section file to store lanczo iteration vectors

    if((filePtr = fopen(sectFile,"wb")) == NULL) {
      printf("\n\nRank%d: Error in function storeLancVec():",Rank);
      printf("\nNot allowed to create the file %s\n", sectFile);
      MPI_Abort(MPI_COMM_WORLD,Rank);
    }
    fclose(filePtr);
  }
  if((lancVecNum/maxLancVecPerFile) >= sectFileNo) {

    // create a new section file

    sectFileNo++;
    sprintf(sectFile,"%s%d", basFile, sectFileNo);
    if((filePtr = fopen(sectFile,"wb")) == NULL)  { 
      printf("\n\nRank%d: Error in function storeLancVec():",Rank);
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
    printf("\nRank%d: Error in function appendLancIterationVecToFile();",Rank);
    printf("\nWrong file = %s to store a vector no %d\n",filename, n);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }

  // vector number

  if(fwrite((const void *)&n,(size_t) sizeof(int), 1, filePtr) != (size_t)1) {
    printf("\n\nRank%d: Error in function appendLancVecToFile()",Rank);
    printf("\nIn writing vector number =  %d", n);
    printf("\nto file %s\n",filename);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }   // end of if-test

  // dimension 

  if(fwrite((const void *)&dim,(size_t)sizeof(int),1,filePtr)!= (size_t)1) { 
    printf("\n\nRank%d: Error in function appendLancVecToFile()",Rank);
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
                              (size_t) dim, filePtr) != (size_t)dim) { 
    printf("\n\nRank%d: Error in function appendLancVecToFile():",Rank);
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

static void lancIterationProc(GR_BAS *grBas,GROUP *group,
                                     LANC_PROC *lancProc)
{
  char      *func = {"lancIterationProc(): "};
  int       iterate, k, control = 0, *rcounts, *displs;
  double    diagElem, *initVec, *finalVec,
            *finalPtr, *temp;

  // local data for the MPI_Allgatherv process

  rcounts = MALLOC(NumProcs, int, func,"rcounts[]");
  displs  = MALLOC(NumProcs, int, func,"displs[]");

  for(k = 0; k < NumProcs; k++) {
    rcounts[k] = grBas->sdBlock[k].totNumSD;
    displs[k]  = grBas->sdBlock[k].startSDnum;
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

    // initialization of new lanczos vector
  
    finalPtr = finalVec;
    for(k = 0; k < grBas->totNumSD_ZN_oscLim; k++) {
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

    // time control for read <SD'|H|SD> from file

    wallClock(7,0); // initialization
    cpuClock(7,0);

    pnMatrixVectorCalc(grBas->typeProc, grBas,group, initVec,finalVec);

    wallClock(6,2);  // stop clock
    cpuClock(6,2);
    
         /*
	 ** Add together all parts of |lancVec> and 
	 ** store the total |lancVec> in all processes
	 */

    MPI_Allreduce(MPI_IN_PLACE,finalVec,grBas->totNumSD_ZN, MPI_DOUBLE,
                                          MPI_SUM,MPI_COMM_WORLD); 

    // all processes contribute  to diagonal energy matrix element

    diagElem = scalarProductOf_pnVectors(grBas->sdBlock[Rank].totNumSD,
			    initVec  + grBas->sdBlock[Rank].startSDnum,
			    finalVec + grBas->sdBlock[Rank].startSDnum);

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
  
    if(!(orthogonalizationProc
       (grBas->sdBlock[Rank].totNumSD, diagElem,lancProc,
        initVec  + grBas->sdBlock[Rank].startSDnum,
        finalVec + grBas->sdBlock[Rank].startSDnum,
                      lancProc->storedLancVec - 1))) {
      lancProc->run_code = 2;  // no more lanczo vector available
      break;
    }
  
    // Using MPI_Allgatherv produce a total |finalVec> in all processes

    MPI_Allgatherv(finalVec + grBas->sdBlock[Rank].startSDnum,
                   grBas->sdBlock[Rank].totNumSD, MPI_DOUBLE,
                   finalVec, rcounts, displs,MPI_DOUBLE,
                   MPI_COMM_WORLD);   

    // new lanczo vector found - append to file

    storeLancVec(lancProc->title,lancProc->storedLancVec,
	                      grBas->sdBlock[Rank].totNumSD,
	         finalVec + grBas->sdBlock[Rank].startSDnum);

    lancProc->storedLancVec++;

    // time control after each Lanczos iteration

    wallClock(5,2);  // stop clock
    cpuClock(5,2);

    // Time controll for lanczo iteration process
  
    if(Rank == MASTER) {
      if(   (iterate == 1) 
	    || ((((iterate + 1)/25)*25)==(iterate + 1)))  {
	writeTimeControl(lancProc,wallClock(5,3),cpuClock(5,3), 
                                  wallClock(6,3),cpuClock(6,3),
                                 wallClock(7,3),cpuClock(7,3));
      }
    } // end MASTER output

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

}  // End: function lancIterationProc()

       /*
       ** The function                                       
       **         scalarProductOf_pnVectors()
       ** calculates and returns the scalar product of 
       ** two vectors stored vec1 and vec2.
       */

static double scalarProductOf_pnVectors(int dim, double *vec1, 
                                                 double *vec2)
{
  int      loop;
  double   result;

  result = D_ZERO; // initalization

  for(loop = 0; loop < dim; loop++) {
    result += (*vec1) * (*vec2);
    vec1++;
    vec2++;
  }
  return result;

} // End: function  scalarProductOf_pnVectors()

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

static int orthogonalizationProc(int dim, double diagElem, 
                                 LANC_PROC *lancProc,
                                 double *initVec, double *finalVec,
                                 int n)
{
  int      k;
  double   factor, sqFactor;

       /*
       ** orthogonalize to |initVec>
       **   (|finalVec[]> - diag_elem * |initVec[]>)  
       **                ---> |finalVec[]>  
       ** Each process modifies its local part of |finalVec>
       */

  add_a_vector(- diagElem,dim,initVec,finalVec);

      /*
      ** Orthogonalize |finalVec> to all previous Lanczo vectors
      ** Again, each process modifies its part of |lancVec>
      */

  for(k = 0; k < n; k++)  {

    // read local part of lancVec(k) from file into |initVec[]>

    readLancVec(k,dim, lancProc, initVec);

    // calculate factor = <initVec|finalVec>

    factor =  scalarProductOf_pnVectors(dim, initVec,finalVec);

    // produce total normalization factor over all processes

    MPI_Allreduce(MPI_IN_PLACE,&factor, 1, MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD); 

    if(fabs(factor) < ZERO_LIMIT){

      // no contribution from current |lancVec(k)>

      continue;
    }
        /*
	** Orthogonalize |finalVec> to lancVec[k]>  
	**   |finalVec[]> <--- |finalVec> - factor * |initVec>
        ** Again, each Rank modifies its part of |lancVec>
	*/
 
    add_a_vector(-factor,dim,initVec,finalVec);


  } // end of loop k

  // close last sectFiles for lancVec: when k = -1

  if(k > 0) readLancVec(-1,dim, lancProc, initVec);

        /*
	** Normalize |finalVec[]> over all processes
	** Calculate off-diagonal matrix element
	**         <n|H|n-1>
	** The new lanc_vec must  have  norm > MATRIX_LIMIT.    
	*/

  sqFactor = SquareNorm_of_a_vector(dim, finalVec);

  MPI_Allreduce(MPI_IN_PLACE,&sqFactor, 1, MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);

  factor = sqrt(sqFactor);

  if(factor < MATRIX_LIMIT)  return FALSE;

  lancProc->h_matrix[((n + 1) * (n + 4))/2 - 1] = factor;

  // Normalize |final_vec> over all processes

  factor = 1.0 / factor;

  scale_a_vector(dim, factor, finalVec);

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

    sprintf(basFile,"%s%d%s%s","Rank",Rank,
            lancProc->title,LANCZO_VECTORS);

    // maximum lanczos vector in a single file

    maxLancVecPerFile 
          = (  (MAX_FILE_SIZE * M_BYTES)
	     / (dim * sizeof(float)+ 2*sizeof(int)));
    sectFileNo = 1;  // initialization
    sprintf(sectFile,"%s%d", basFile, sectFileNo);

    // open first lanczo iteration sectFile

    if((filePtr = fopen(sectFile,"rb")) == NULL) {
      printf("\n\nRank%d: Error in function readLancVec():",Rank);
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
      printf("\n\nRank%d: Error in function readLancVec():",Rank);
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
  if(fread((void *)&n,(size_t) sizeof(int), 1, *filePtr) != (size_t)1) {
    printf("\n\nRank%d: Error in function readLancVecFromOpenFile()",Rank);
    printf("\nIn reading vector number =  %d", n);
    printf("\nfrom file %s\n",filename);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }   // end of if-test

  // test vector number

  if(n !=  lancVecNum) {
    printf("\n\nRank%d: Error in function readLancVecFromOpenFile()",Rank);
    printf("\nWrong vector number read");
    printf("\nread number = %d - should be lancVecNum = %d\n",
                                                n,lancVecNum);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }

  // read vector dimension
  if(fread((void *)&n,(size_t) sizeof(int), 1, *filePtr) != (size_t)1) {
    printf("\n\nRank%d: Error in function readLancVecFromOpenFile()",Rank);
    printf("\nIn reading vector dimension =  %d", n);
    printf("\nfrom file %s\n",filename);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }   // end of if-test

  // test vector dimension
  if(n !=  dim) {
    printf("\n\nRank%d: Error in function readLancVecFromOpenFile()",Rank);
    printf("\nWrong vector dimension read");
    printf("\nread number = %d - should be lancVecNum = %d\n",
                                                n,dim);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }

  // read vector float amplitudes
  if(fread((void *)floatMem,(size_t)sizeof(float), (size_t) dim,
                                       *filePtr) != (size_t) dim)  {
    printf("\n\nRank%d: Error in function readLancVecFromOpenFile():",Rank);
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
       **          SquareNorm_of_a_vector()                   
       ** calculates the square norm of a vector stored in 
       ** vec[] with iteration->tot_dimSD elements. 
       ** The function returns 
       **     result = sqrt( Sum(x_i * x_i))
       ** Possible time-reversal symmetry is included   
       */

static double SquareNorm_of_a_vector(int dim, double *vec)
{
  int      loop; 
  double   *ptr_vec, result;

  result  = D_ZERO;
  ptr_vec = vec;

  for(loop = 0; loop < dim; loop++, ptr_vec++) {
    result += (*ptr_vec) * (*ptr_vec);
  }

  return result;

}  // End: function SquareNorm_of_a_vector()

    /*
    ** The function                                           
    **              saveFinalData()        
    ** diagonalizes the energy matrix in Lanczos basis and
    ** transforms eigenvectors into a |SD> basis
    ** Save Eigenvalue, ang-mom and eigenvector on file
    */

static void saveFinalData(GR_BAS *grBas,GROUP *group, LANC_PROC *lancProc)
{
  char         *func = {"saveFinalData(): "}, sectFile[ONE_LINE];
  int          calcCM,dim_h, loop, k;
  double       diagElem, *vecPtr, *initVec, *finalVec, *j_occ[2];
  EIGEN        *listEigen;
  EIGEN_DATA   h_eigenData;

  // TEST POINTERS  

  double   *START_Ptr,  *STOP_Ptr;
 

     /*
     ** Calculate <SD'|CM|SD> and store in memory
     ** where <SD'|VEFF|SD> are located
     */

  if((grBas->typeProc == VEFF_LANC) &&(grBas->calc_CM == YES)) {
    calcCM = YES;
    pnStoreSD_MatrElem(CM_CALC,calcCM,grBas,group);
  }
     /*
     ** Calculate <SD'|J**2|SD> and store in the free
     ** memory following <SD'|CM|SD>
     */

     /* Step 1:
     ** Calculates and stores in 
     ** grBas->op_pp[], grBas->op_nn[] and grBas->op_pn[]
     ** the complete set of m-scheme two-particle matrix elements
     ** of angular momentum operator J**2.
     */

  if(grBas->spBas[PROTON]->part >= 2) {
    id_two_part_angJ(grBas->spBas[PROTON],&grBas->op_pp);
  }
  if(grBas->spBas[NEUTRON]->part >= 2) {
    id_two_part_angJ(grBas->spBas[NEUTRON],&grBas->op_nn);
  }
  pn_two_part_angJ(grBas->parZ,grBas->spBas[PROTON],
              grBas->spBas[NEUTRON],&grBas->op_pn);  

     /* Step 2:
     ** Calculate and store in memory all 
     ** diagonal proton-neutron matrix elements
     **     <SD(Z)g,p;SD(N)g,n|J**2|SD(Z)g,p;SD(N)g,n>
     ** and as many as possible nondiag proton-neutron
     ** matrix elements 
     **     <SD(Z)g',p';SD(N)g',n'|J**2|SD(Z)g,p;SD(N)g,n>
     ** Next, calculate nondiagonal proton-neutron 
     ** matrix elements
     **     <SD(Z)g',p';SD(N)g',n'|J**2|SD(Z)g,p;SD(N)g,n> 
     */

  pnStoreSD_MatrElem(ANG_CALC,grBas->calc_CM,grBas,group);

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

  // memory for two |lancVec>

  initVec  = lancProc->vec1;
  finalVec = lancProc->vec2;

  // memory for single-particle occopation numbers

  j_occ[PROTON]  = MALLOC(grBas->spBas[PROTON]->numj_orb,double, func,
                                                      "j_occ[PROTON]"); 
  j_occ[NEUTRON] = MALLOC(grBas->spBas[NEUTRON]->numj_orb,double, func,
			                             "j_occ[NEUTRON]"); 
  for(loop = 0; loop < lancProc->states; loop++) {

    // calculate eigenvector in sd-basis

    for(vecPtr = finalVec,k = 0;k < grBas->totNumSD_ZN_oscLim;vecPtr++,k++) {
      *vecPtr = D_ZERO;
    }

    // transfer eigenvec with number loop to SD basis

    lancToSD_Bas(lancProc,grBas->sdBlock[Rank].totNumSD, 
       initVec,finalVec + grBas->sdBlock[Rank].startSDnum,
                      dim_h, listEigen[loop].vector);

    MPI_Allreduce(MPI_IN_PLACE,finalVec,grBas->totNumSD_ZN, MPI_DOUBLE,
                                          MPI_SUM,MPI_COMM_WORLD); 

    h_eigenData.vecNum            = loop; 
    h_eigenData.dim               = grBas->totNumSD_ZN_oscLim; 
    h_eigenData.numj_occ[PROTON]  = grBas->spBas[PROTON]->numj_orb;
    h_eigenData.numj_occ[NEUTRON] = grBas->spBas[NEUTRON]->numj_orb;
    h_eigenData.eigenVal          = listEigen[loop].value;

    if(grBas->calc_CM == 1) {

         /*
         ** do
         **     R**2|finalVec> = |initVecVec>
         ** where the energy operator R**2 specified through
         ** its matrix elements in m-scheme.
	 ** Contribution to |initVec> is distributed over 
	 ** through all processes
         */

      for(vecPtr = initVec,k = 0;k < grBas->totNumSD_ZN_oscLim;vecPtr++,k++) {
	*vecPtr = D_ZERO;
      }

      pnMatrixVectorCalc(CM_CALC,grBas,group,finalVec,initVec);

         /*
	 ** Add together all parts of |initVec> and 
	 ** store the total |initVec> in all processes
	 */

      MPI_Allreduce(MPI_IN_PLACE,initVec,grBas->totNumSD_ZN, MPI_DOUBLE,
                                          MPI_SUM,MPI_COMM_WORLD);

      diagElem 
	= scalarProductOf_pnVectors(grBas->sdBlock[Rank].totNumSD,
				    initVec  + grBas->sdBlock[Rank].startSDnum,
                                    finalVec + grBas->sdBlock[Rank].startSDnum);

      // sum_of(diagElem) --> diagElem and stores in all processes

      MPI_Allreduce(MPI_IN_PLACE,&diagElem,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD); 

      h_eigenData.valCM  = diagElem;
    }
         /*
         ** do
         **     J**2|finalVec> = |initVecVec>
         ** where the energy operator J**2 specified through
         ** its matrix elements in m-scheme.
	 ** Contribution to |initVec> is distributed over 
	 ** through all processes
         */

    for(vecPtr = initVec, k = 0; k < grBas->totNumSD_ZN_oscLim; vecPtr++, k++) {
      *vecPtr = D_ZERO;
    }

    pnMatrixVectorCalc(ANG_CALC,grBas,group,finalVec,initVec);

         /*
	 ** Add together all parts of |initVec> and 
	 ** store the total |initVec> in all processes
	 */

    MPI_Allreduce(MPI_IN_PLACE,initVec,grBas->totNumSD_ZN, MPI_DOUBLE,
                                          MPI_SUM,MPI_COMM_WORLD);
    diagElem 
      = scalarProductOf_pnVectors(grBas->sdBlock[Rank].totNumSD,
	    initVec  + grBas->sdBlock[Rank].startSDnum,
	    finalVec + grBas->sdBlock[Rank].startSDnum);

   // sum_of(diagElem) --> diagElem and stores in all processes

    MPI_Allreduce(MPI_IN_PLACE,&diagElem,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD); 

    h_eigenData.angMom  = diagElem;


    // calculate list of single-particle occupation numbers

    proton_j_occupation(grBas,group,j_occ[PROTON],finalVec);


    // Summed j_occ[PROTON]occupation elements

    MPI_Allreduce(MPI_IN_PLACE,j_occ[PROTON],grBas->spBas[PROTON]->numj_orb,
                                          MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD); 

    neutron_j_occupation(grBas,group,j_occ[NEUTRON],finalVec);

    MPI_Allreduce(MPI_IN_PLACE,j_occ[NEUTRON],grBas->spBas[NEUTRON]->numj_orb,MPI_DOUBLE,
                                                           MPI_SUM,MPI_COMM_WORLD); 

   // save energy eigenvalue and  eigenvector on file


    if(Rank == MASTER) {

      storeEigenVec(lancProc,h_eigenData,j_occ[0],j_occ[1],finalVec);
    }

  }   // end loop:  all eigenvector written to file

  // remove all files containing |lancVec>

  k = 1;  
  for( ; ;) {
    sprintf(sectFile,"%s%d%s%s%d","Rank",Rank,lancProc->title,
                                    LANCZO_VECTORS, k++);
    if(remove(sectFile) != 0) break;
  }

  // remove all files containing <SD|ANG|SD>

  sprintf(sectFile,"%s%s%d%s%s",lancProc->title,"Rank",Rank,SD_STORE_ANG,DIAG);
  remove(sectFile);
  k = 1;  
  for(  ; ;) {
    sprintf(sectFile,"%s%s%d%s%s%d",lancProc->title,
                     "Rank",Rank,SD_STORE_ANG,NONDIAG, k++);
    if(remove(sectFile) != 0) break;
  }
 
  free(j_occ[NEUTRON]);  // remove local memory
  free(j_occ[PROTON]);
  for(loop = dim_h - 1; loop >= 0;loop--) {
    free(listEigen[loop].vector); 
  }
  free(listEigen);

} // End: function saveFinaData()

      /*
      ** The function                                        
      **          lanczoToSD_Basis()                       
      ** transform an eigenvector in Lanczos vector basis 
      ** into SD basis and store the result in |finalVec>   
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
    **        proton_j_occupation()                   
    ** calculate the j-occupation numbers  of protons for a shell model
    ** eigenstate eigenVec[]
    */

void proton_j_occupation(GR_BAS *grBas, GROUP *group,
                                 double *j_occ, double *eigenVec)
{
  int       k,pStep, nStep,orb, part, currGr,locBitNo;
  ULL       pos, sd_val, *protonSD;
  double    *startGrAmplPtr, probability, *amplPtr;

  // initialization
  
  for(k = 0; k < grBas->spBas[PROTON]->numj_orb; k++) j_occ[k] = 0.0;


  // run through all proton-neutron groups

  for(currGr = 0; currGr < grBas->num_gr; currGr++, group++)  {

    if((group->numSD_ZN_oscLim) == 0) continue;  // no basis states

    // points to start eigenvector amplitude fir current group

    startGrAmplPtr = eigenVec + group->startAmp;


    // run through all |SD(PROTON)> in current group 

    locBitNo = 0;

    for(pStep = 0,protonSD = group->SD[PROTON];pStep < group->numSD[PROTON];
	                protonSD++, pStep++, locBitNo += group->numSD[NEUTRON]){

      amplPtr = startGrAmplPtr + Rank;

      // sum ampl**2 for each pStep

      probability = 0.0;
      for(nStep = Rank; nStep < group->numSD[NEUTRON];nStep += NumProcs) { 
	if(BARR_TEST(group->bitArray,locBitNo + nStep)) continue;
	probability += ((*amplPtr) * (*amplPtr));
	amplPtr += NumProcs;
      }

      //      if(probability < MATRIX_LIMIT) continue;
		              
      sd_val = *protonSD; // orbital initialization
      orb    = 0;
      pos    = ULL_ONE;
      part   = 0;
      do   {
	for( ;!(sd_val & pos); orb++, pos <<= 1); // occ. orbit found
	j_occ[grBas->spBas[0]->mbas[orb].orb] += probability;
	orb++;
	pos <<= 1;
      } while(++part < grBas->spBas[PROTON]->part);      


      startGrAmplPtr += group->numSD[NEUTRON];

    } // end pStep
  } // end all groups
}  // End: function proton_j_occupation()

   /*
   ** The function                            
   **        neutron_j_occupation()                   
   ** calculate the j-occupation numbers  of neutrons for a shell model 
   ** eigenstate eigenVec[]
   */

  void neutron_j_occupation(GR_BAS *grBas, GROUP *group,
                                    double *j_occ, double *eigenVec)
{
  int        k, pStep,nStep,orb, part, currGr, locBitNo;
  ULL        pos, sd_val, *neutronSD;
  double     probability, *amplPtr;
  double     *startAmplPtr;

  for(k = 0; k < grBas->spBas[NEUTRON]->numj_orb; k++) j_occ[k] = 0.0;
  
  // run through all proton-neutron groups

  for(currGr = 0; currGr < grBas->num_gr; currGr++, group++)  {

    if((group->numSD_ZN_oscLim) == 0) continue;  // no basis states

    startAmplPtr = eigenVec + group->startAmp;
    locBitNo = 0;
    for(pStep = 0; pStep < group->numSD[PROTON]; pStep++,
                 locBitNo += group->numSD[NEUTRON]) {

      neutronSD = group->SD[NEUTRON] + Rank;
      amplPtr   = startAmplPtr + Rank;
      for(nStep = Rank; nStep < group->numSD[NEUTRON]; nStep += NumProcs,
                                             neutronSD += NumProcs) {
	if(BARR_TEST(group->bitArray,locBitNo + nStep)) continue;

	probability = ((*amplPtr) * (*amplPtr));

	amplPtr += NumProcs;

	sd_val = *neutronSD;
	pos    = ULL_ONE;
	orb    = 0;
	part   = 0;
	do   {
	  for( ;!(sd_val & pos); orb++, pos <<= 1); // occ. orbit found
      
	  j_occ[grBas->spBas[NEUTRON]->mbas[orb].orb] += probability;
	  orb++;
	  pos <<= 1; 
	} while(++part < grBas->spBas[NEUTRON]->part);
      } // end nStep

      startAmplPtr += group->numSD[NEUTRON];

    } // end pStep

  }    // loop through all groups

} // End: function neutron_j_occupation()

    /*
    ** The function 
    **     storeEigenVec()
    ** opens a sect file and stores a lanczos 
    ** eigenvector, If current sect file has 
    ** reacded MAX_FILE_SIZE current sect file 
    ** is closed and a new is opened
    */

static void storeEigenVec(LANC_PROC *lancProc, EIGEN_DATA h_eigenData,
                          double *j_occZ,double *j_occN,double *vec)
{
  static char   basFile[ONE_LINE], sectFile[ONE_LINE];
  static int    maxLancVecPerFile, sectFileNo;
  FILE          *filePtr;

  if(h_eigenData.vecNum == 0) { // initialization

    sprintf(basFile,"%s%s",lancProc->title,EIGEN_VECTORS);

    // maximum eigenvectors in a single file

    maxLancVecPerFile 
          = (  (MAX_FILE_SIZE * M_BYTES)
	       / (h_eigenData.dim * sizeof(float) + sizeof(EIGEN_DATA)));
    sectFileNo = 1;  // initialization
    sprintf(sectFile,"%s%d", basFile, sectFileNo);

    // create a new section file to store lanczo eigenvectors

    if((filePtr = fopen(sectFile,"wb")) == NULL) {
      printf("\n\nRank%d: Error in function storeEigenVec():",Rank);
      printf("\nNot allowed to create the file %s\n", sectFile);
      MPI_Abort(MPI_COMM_WORLD,Rank);
    }
    fclose(filePtr);
  }
  if((h_eigenData.vecNum/maxLancVecPerFile) >= sectFileNo) {

    // create a new section file

    sectFileNo++;
    sprintf(sectFile,"%s%d", basFile, sectFileNo);
    if((filePtr = fopen(sectFile,"wb")) == NULL)  { 
      printf("\n\nRank%d: Error in function storeEigenVec():",Rank);
      printf("\nNot allowed to create the file %s\n", sectFile);
      MPI_Abort(MPI_COMM_WORLD,Rank);
    }
    fclose(filePtr);
  }
  appendEigenVecToFile(sectFile, h_eigenData, j_occZ,j_occN, vec); 

} // End: function storeEigenVec()

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

static void appendEigenVecToFile(char *filename, EIGEN_DATA h_eigenData,
               double *j_occZ,double *j_occN,double *vector) 
{
  char    *func = {" appendEigenVecToFile(): "};
  int     k;
  float   *tempMem, *finalPtr;
  double  *initPtr;
  FILE    *filePtr;

  if((filePtr = fopen(filename,"ab+")) == NULL) {
    printf("\n\nRank%d: Error in function appendEigenVecToFile();",Rank);
    printf("\nWrong file = %s to store a vector no %d\n",
                              filename, h_eigenData.vecNum);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }
  // struct h_eigenData

  if(fwrite((const void *)&h_eigenData,(size_t) sizeof(EIGEN_DATA), 1, 
                                                filePtr) != (size_t)1) {
    printf("\n\nRank%d: Error in function appendEigenVecToFile()",Rank);
    printf("\nIn writing structure h_eigenData");
    printf("\nto file %s\n",filename);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }   // end of if-test

  // list of j_occupation data

  if(fwrite((const void *)j_occZ,(size_t) sizeof(double), 
	    h_eigenData.numj_occ[0],filePtr) != (size_t)h_eigenData.numj_occ[0]) {
    printf("\n\nRank%d: Error in function appendEigenVecToFile()",Rank);
    printf("\nIn writing %d proton singel-particle occupation  data",
                                                  h_eigenData.numj_occ[0]);
    printf("\nto file %s\n",filename);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }   // end of if-test


  if(fwrite((const void *)j_occN,(size_t) sizeof(double), 
               h_eigenData.numj_occ[1],filePtr) != (size_t)h_eigenData.numj_occ[1]) {
    printf("\n\nRank%d: Error in function appendEigenVecToFile()",Rank);
    printf("\nIn writing %d neutron singel-particle occupation  data",
                                                  h_eigenData.numj_occ[1]);
    printf("\nto file %s\n",filename);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }   // end of if-test

  tempMem = MALLOC(h_eigenData.dim, float, func, "tempMem");
  initPtr = vector;
  finalPtr = tempMem;
  for(k = 0; k < h_eigenData.dim; k++)  {
    *finalPtr = (float)(*initPtr);
    initPtr++;
    finalPtr++;
  }
  if(fwrite((const void *)tempMem,(size_t) sizeof(float), 
                (size_t) h_eigenData.dim, filePtr) != (size_t)h_eigenData.dim) { 
    printf("\n\nRank%d: Error in function appendEigenVecToFile():",Rank);
    printf("\nIn writing eigenvector no %d diag to file %s\n\n",
                               h_eigenData.dim, filename);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }
  free(tempMem);
 
  fclose(filePtr);

} // End: function appendEigenVecToFile()

     /*
     ** The function
     **     writeTimeControl()
     ** time control data to output file 
     */

static void writeTimeControl(LANC_PROC *lancProc,
                        TID wallTime, TID cpuTime,
                        TID wallTimeA, TID cpuTimeA,
                        TID wallTimeB, TID cpuTimeB)
{   
  char         filename[ONE_LINE];
  FILE         *filePtr;

  sprintf(filename,"%s%s",lancProc->title, RESULT_OUTPUT);

  if((filePtr = fopen(filename,"a")) == NULL)   {
    printf("\n\nRank%d: Error in function  lancIterationProc():",Rank);
    printf("\nWrong file = %s to open the output data file\n",filename);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }

    fprintf(filePtr,"\n\nStored no %d energy lanczo vectors",
                                 lancProc->storedLancVec);
    fprintf(filePtr,"\nWall time: %llu hour %llu min %llu sec",
                           wallTime.hour, wallTime.min, wallTime.sec);
    fprintf(filePtr,"\nCPU  time: %llu hour %llu min %llu sec",
                            cpuTime.hour,cpuTime.min, cpuTime.sec);

    // Print time control for H() * Vec()

    fprintf(filePtr,"\nTest time H()*Vec()");
    fprintf(filePtr,"\nWall time: %llu hour %llu min %llu sec",
                      wallTimeA.hour, wallTimeA.min, wallTimeA.sec);
    fprintf(filePtr,"\nCPU  time: %llu hour %llu min %llu sec",
                     cpuTimeA.hour, cpuTimeA.min, cpuTimeA.sec);

    // Print time control for read <SD'|HSD> from file 

    fprintf(filePtr,"\nTest time file reading");
    fprintf(filePtr,"\nWall time: %llu hour %llu min %llu sec",
                     wallTimeB.hour, wallTimeB.min, wallTimeB.sec);
    fprintf(filePtr,"\nCPU  time: %llu hour %llu min %llu sec  ",
                     cpuTimeB.hour, cpuTimeB.min, cpuTimeB.sec);
 
    fflush(stdout);
    fclose(filePtr);

} // End: function writeTimeControl()
