/*******************  The module lanc-id-process.c  ******************/

#include "shell-model.h"

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
     ** reserves dynamic memory to bu used
     ** in the Lanczos iteration process to 
     ** obtain shell model eigenvalues and 
     ** the corresponding eigenvectors
     */

static void initialization_process(char *type_calc, SP_BAS *spBas, SD_BAS *sdBas, 
                                                                  double *initVec);
     /*
     ** returns a start vector for a Lanczos iteration 
     ** process in vec[].
     */

static void random_start_vector(SP_BAS *spBas, SD_BAS *sdBas, double *init_vec);
     /*
     ** calculates and return a normalized start vector where
     **  amplitudes are selected through a random procedure
     */

 void storeLancVec(int lancVecNum, int dim, 
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

static void lancIterationProc(SP_BAS *spBas, SD_BAS *sdBas,
			      LANC_PROC *lancProc);
     /*
     ** takes a start lanczos vector stored in |vec1> and 
     ** performs a Lanczos iteration process 
     ** to convergence is obtained
     */

static int orthogonalizationProc(LANC_PROC *lancProc,SD_BAS *sdBas, 
				double *initVec, double *finalVec, int n);
     /*
     ** orthogonalizes the new vector iteration->final_vec to all previous 
     ** Lanczos vectors, calculates a new off diagonal matrix element and 
     ** finally normalize the new vector.
     ** if(norm > matrix_limit) the function return TRUE. Otherwise FALSE
     ** indicating that no more Lanczos vectors are available 
     */

void readLancVec(int lancVecNum, int dim, 
		 LANC_PROC *lancProc, double *vec);
     /*
     ** opens a sect file and reads a lanczos iteration
     ** vector,If current sect file has reacded 
     ** MAX_FILE_SIZE current sect file is closed and 
     ** a new is opened
     */

void readLancVecFromOpenFile(char *fileName, FILE **filePtr, int fileNum, 
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

static void saveFinalData(SHELL *model, SP_BAS *spBas,
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

static double scalar_product_of_id_vectors(SD_BAS *sdBas, double *vec1, double *vec2);
     /*
     ** calculates and returns the scalar product of 
     ** two vectors stored vec1 and vec2.
     ** Possible time reversal symmetry is included
     */

static double norm_of_a_vector(SD_BAS *sdBas, double *vec);
     /*
     ** calculates the norm of a vector stored in 
     ** vec[] with iteration->tot_dimSD elements. 
     ** The function returns 
     **     result = sqrt( Sum(x_i * x_i))
     ** Possible time-reversal symmetry is included   
     */

static void CM_energy_matr_elem(LANC_PROC *lancProc,SP_BAS *spBas,
                                  SD_BAS *sdBas, double *init_vec, 
                           double *final_vec, int stored_lanc_vec);
     /*
     ** calculates nondiagonal center_of_Mass matrix 
     ** elements 
     **      <lanc_vec|CM_INTfinal_vec>
     ** for the stored lanc_vec below |final_vec> and
     ** store the result in one-dimesional CM_matrix[].
     */
static void  angularMomentumInteraction(SHELL *model, SP_BAS *spBas, 
                                 SD_BAS *sdBas, SD_STORE *sdStoreAng);
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

void storeEigenVec(LANC_PROC *lancProc, EIGEN_DATA h_eigenData,
                                    double *j_occ,double *vec);
   /*
   ** opens a sect file and stores a lanczos 
   ** eigenvector,If current sect file has 
   ** reacded MAX_FILE_SIZE current sect file 
   ** is closed and a new is opened
   */

void appendEigenVecToFile(char *filename, EIGEN_DATA h_eigenData,
                                    double *j_occ, double *vector);
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

static void id_j_occupation(SP_BAS *spBas, SD_BAS *sdBas, 
			    double *j_occ, double *eigenVec);
   /*
   ** calculate the j-occupation numbers  of an identicle particle shell model
   ** eigenstate from the lanczo procedure and returned the result in j_occ[]. 
   */

               /**** End: function declarations ****/

               /**** The function definitions  ****/ 
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

void id_eigenvalue_Lanczos_search(SHELL *model, LANC_PROC *lancProc, 
                                  SP_BAS *spBas, SD_BAS *sdBas)
{
  TID    ex_time;

  // reserve local memory for lanczo process

  localLanczoMem(0, sdBas->tot_dimSD, lancProc);

  // normalized start vector for Lanczos iteration process

  initialization_process(lancProc->type_calc, spBas, sdBas, lancProc->vec1);

  // |init_lancVec> = lancProc->vec1  found  - append to file 

  lancProc->storedLancVec = 0; // initialization
  storeLancVec(lancProc->storedLancVec++,sdBas->tot_dimSD,
                                  lancProc,lancProc->vec1);

  time_step(0, 1);   // time control for total lanczos process

  lancIterationProc(spBas,sdBas,lancProc);

  // time control for the lanczos process

  ex_time = time_step(0,2); 
  printf("\n\nTotal time used for Lanczos iteration procedure");
  printf("\nrun time: %lld hour %lld min %lld sec\n",ex_time.hour, 
                                          ex_time.min, ex_time.sec);
  fflush(stdout);

  // calculate and save final data

  saveFinalData(model, spBas, sdBas, lancProc);

  // release local memory for lanczo process

  localLanczoMem(1, sdBas->tot_dimSD, lancProc);
  
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
    lancProc->vec1 = MALLOC(dimSD, double,func,"lanc_vec1[]");
    lancProc->vec2 = MALLOC(dimSD, double,func,"lanc_vec2[]");

    // local memory for the Lanczos energy matrix

    lancProc->h_matrix 
      = CALLOC(((lancProc->max_iterations + 1) *(lancProc->max_iterations + 2))/2, 
                                   double, func, "h_matrix[]");
    
    if(lancProc->calc_CM == 0) {
      lancProc->CM_matrix = NULL;
    }
    else {
      lancProc->CM_matrix 
	= CALLOC(( (lancProc->max_iterations + 1) 
		   *(lancProc->max_iterations + 2))/2,
		         double, func, "CM_matrix[]");
    } 

    // local memory to save dynamical changes in eigenvalues

    lancProc->delta_eigen = MALLOC(5 * lancProc->states, double, 
                                         func, "delta_eigen[]");
  } // end reserve memory

  else {
    free(lancProc->delta_eigen);    // remove local lanczos local memory 
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
      **   initialization_process();
      ** returns a start vector for a Lanczos iteration 
      ** process in vec[].
      */

static void initialization_process(char *type_calc, SP_BAS *spBas, SD_BAS *sdBas,
                                                                    double *initVec)
{
  switch(strlen(type_calc)) {
      case 12 :                                   /* random-start */
               random_start_vector(spBas, sdBas, initVec);
               break;
      case 15:                                   /* random-continue */
               printf("\n\nError in function initialization_process():");
               printf("\ntype_calc = %s", type_calc);
               printf("\n Not implemented!!!!!!!!!\n");
               exit(1);
               break;
      case 16 :                                   /* fixed-J-continue */
               printf("\n\nError in function initialization_process():");
               printf("\ntype_calc = %s\n",type_calc);
               printf("\n Not implemented!!!!!!!!!\n");
               exit(1);
               break;
      default :
               printf("\n\nError in function initialization_process():");
               printf("\nWrong type_calc = %s",type_calc);
               printf("\n Not implemented!!!!!!!!!\n");
               exit(1);
               break;
  } /* end switch() */

} /* End: function initialization_process() */

      /*
      ** The function                           
      **        random_start_vector()                
      ** calculates and return a normalized start vector in vec[] 
      ** where amplitudes are selected through a random procedure
      */

static void random_start_vector(SP_BAS *spBas, SD_BAS *sdBas, double *initVec)
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
      *(ptrVec++) = (value = 1.0+((double)loop)/dim,value*phase);
      norm +=  value * value;

/*********************
      phase = (rand() & TRUE) ? +1 : -1;
***********************/
 
      phase = -phase;

    }
    norm *= 2.0;     /* contr from time-reversal comp */

    for( ; loop < dim; loop++) {
      *(ptrVec++) = (value = 1.0+((double)loop)/dim,value*phase);
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
      exit(1);
    }
    norm = 1.0 / sqrt(norm);
    ptrVec = initVec;
    for(loop = 0; loop < dim; loop++)  {
      *(ptrVec++) *= norm;
    }

  }  // end case: time-reversal symmetry

  else  {                            /* case: no time-reversal symmetry */

    for(loop = 0; loop < dim; loop++) {
      *(ptrVec++) = (value = 1.0+((double)loop)/dim,value*phase);
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
      exit(1);
    }
    norm = 1.0 / sqrt(norm);

    ptrVec  = initVec;
    for(loop = 0; loop < dim; loop++)  {
      *(ptrVec++) *= norm;
    }
  } // end case: no time-reversal symmetry

} // End: function random_start_vector()

    /*
    ** The function 
    **     storeLancVec()
    ** opens a sect file and stores a lanczos iteration
    ** vector,If current sect file has reacded 
    ** MAX_FILE_SIZE current sect file is closed and 
    ** a new is opened
    */

void storeLancVec(int lancVecNum, int dim, 
                         LANC_PROC *lancProc, double *vec)
{
  static char   basFile[ONE_LINE], sectFile[ONE_LINE];
  static int    maxLancVecPerFile, sectFileNo;
  FILE          *filePtr;

  if(lancVecNum == 0) { // initialization

    sprintf(basFile,"%s%s",lancProc->title, LANCZO_VECTORS);

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
      exit(1);
    }
    fclose(filePtr);
  }
  if((lancVecNum/maxLancVecPerFile) >= sectFileNo) {

    // create a new section file

    sectFileNo++;
    sprintf(sectFile,"%s%d", basFile, sectFileNo);
    if((filePtr = fopen(sectFile,"wb")) == NULL)  { 
      printf("\n\nError in function storeLancVec():");
      printf("\nNot allowed to create the file %s\n", sectFile);
      exit(1);
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

void appendLancVecToFile(char *filename, int n, int dim, double *vector) 
{
  char    *func = {" appendLancIterationVecToFile(): "};
  int     k;
  float   *tempMem, *finalPtr;
  double  *initPtr;
  FILE    *filePtr;

  if((filePtr = fopen(filename,"ab+")) == NULL) {
    printf("\nError in function appendLancIterationVecToFile();");
    printf("\nWrong file = %s to store a vector no %d\n",filename, n);
    exit(1);
  }

  // vector number

  if(fwrite((const void *)&n,(size_t) sizeof(int), 1, filePtr) != 1) {
    printf("\nError in function appendLancVecToFile()");
    printf("\nIn writing vector number =  %d", n);
    printf("\nto file %s\n",filename);
    exit(1);
  }   // end of if-test

  // dimension 

  if(fwrite((const void *)&dim,(size_t) sizeof(int),1,filePtr)!=1) { 
    printf("\nError in function appendLancVecToFile()");
    printf("\nIn writing numSD =  %d", dim);
    printf("\nto file %s\n",filename);
    exit(1);
  }  // end of if-test

  tempMem = MALLOC(dim, float, func, "tempMem");
  initPtr = vector;
  finalPtr = tempMem;
  for(k = 0; k < dim; k++)  {
    *(finalPtr++) = (float)(*(initPtr++));
  }
  if(fwrite((const void *)tempMem,(size_t) sizeof(float), 
                              (size_t) dim, filePtr) != dim) { 
    printf("\n\nError in function appendLancVecToFile():");
    printf("\nIn writing %d diag matrix elements to file %s\n\n",
                                                  dim, filename);
    exit(1);
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

static void lancIterationProc(SP_BAS *spBas, SD_BAS *sdBas,
                                             LANC_PROC *lancProc)
{
  int       iterate;
  double    *initVec, *finalVec, *temp;
  TID       ex_time;

         /*  lancProc->run_code:
 	 **    = 0 - reached maximum lanczo iterations
         **    = 1 - Lanczos process has converged
	 **    = 2 - no more Lanczos vectors
	 */ 

  lancProc->run_code = 0; // initialization
  initVec            = lancProc->vec1; 
  finalVec           = lancProc->vec2;

  for(iterate = 1; iterate < lancProc->max_iterations; iterate++) {

    time_step(1, 1);     // time control for each lanczo iterartion

    if(lancProc->calc_CM == 1) {

      // diagonal element: <init_vec|H(CM)|init_vec> 

      lancProc->CM_matrix[( (lancProc->storedLancVec - 1) 
                           *(lancProc->storedLancVec + 2))/2] 
              = id_lanczo_matrix_element(spBas,sdBas,&lancProc->sdStoreCM,
                                                        initVec, initVec);
    }
         /*
         ** do
         **     H|initVec> = |finalVec>
         ** where the energy operator H specified through
         ** its matrix elements in m-scheme.
         */

    id_matrix_vector_calc(spBas, sdBas, &lancProc->sdStoreVeff,
                                           initVec, finalVec);

    // diagonal element in energy matrix

    lancProc->h_matrix[( (lancProc->storedLancVec - 1) 
                        *(lancProc->storedLancVec + 2))/2] 
             = scalar_product_of_id_vectors(sdBas, initVec, finalVec);

         /*
         ** orthogonalize |finalVec>  against all previous Lanczos
	 ** vectors stored on file. Then normalize and return a new 
	 ** normalized Lanzcos vector in |finalVec>
	 */

    if(!orthogonalizationProc(lancProc, sdBas,
             initVec, finalVec, lancProc->storedLancVec - 1))    {
      lancProc->run_code = 2;  // no more lanczo vector available
      break;
    }
    
    if(lancProc->calc_CM == 1) {

      // off-diagonal elements <lanc_vec(k)|H(CM)|final_vec>,
      //          k = 0,....,storedLancVec - 1 
      
      CM_energy_matr_elem(lancProc, spBas, sdBas,
                   initVec, finalVec, lancProc->storedLancVec);
    }

    // new lanczo vector found - append to file

    storeLancVec(lancProc->storedLancVec++, sdBas->tot_dimSD, 
                                            lancProc, finalVec);

    // time control after each Lanczos iteration

    ex_time = time_step(1,2); 
    printf("\nStored %d energy lanczo vectors",lancProc->storedLancVec);
    printf("\nrun time: %lld hour %lld min %lld sec",ex_time.hour,
                                        ex_time.min, ex_time.sec);
    fflush(stdout);

         /*
         ** When dimension of h_matrix has reached the required eigenvalues,
	 ** the energy spectrum is tested against a given criterium.
	 ** If function returns TRUE energy convergence has been reached.
	 ** and the Lanczos process terminates
         */

    if((lancProc->storedLancVec - 1) >= lancProc->states) {
      if(eigen_value_convergence(lancProc)) {
	lancProc->run_code = 1;
	break;
      }
    }
    temp     = initVec; // interchange the pointers
    initVec  = finalVec;
    finalVec = temp;

  } /// end lanczo iteration loop

     /*
     ** Possible modification of the number
     ** of calculated enegry eigenvectors
     **     lancProc->states
     */

  if(lancProc->run_code == 2) {
    lancProc->states = MIN(lancProc->states, lancProc->storedLancVec);
  }
  else {
    lancProc->states = MIN(lancProc->states, lancProc->storedLancVec - 1);
  }

  // remove all files containing <SD|VEFF|SD>

/***************************************
  sprintf(sectFile,"%s%s%s",lancProc->title,SD_STORE_VEFFJ,DIAG);
  remove(sect_file);
  k = 1;  
  for( ; ;) {
    sprintf(sect_file,"%s%s%s%d",lancProc->title,SD_STORE_VEFFJ,NONDIAG, k++);
    if(remove(sect_file)) break;
  }
*****************************/
 
 // remove all files containing <SD|CM_INT|SD>


}  // End: function lancIterationProc()

    /*
    ** The function
    **                  orthogonalizationProc()                        
    ** orthogonalizes the new vector iteration->final_vec to all previous 
    ** Lanczos vectors, calculates a new off diagonal matrix element and 
    ** finally normalize the new vector.
    ** if(norm > matrix_limit) the function return TRUE. Otherwise FALSE
    ** indicating that no more Lanczos vectors are available 
    */

static int orthogonalizationProc(LANC_PROC *lancProc,SD_BAS *sdBas, 
                               double *initVec, double *finalVec, int n)
{
  int      k;
  double   factor;

       /*
       ** orthogonalize to last stored Lanczos vector
       **   (|finalVec[]> - diag_elem * |initVec[]>)  
       **                ---> |finalVec[]>  
       */

  add_a_vector(-lancProc->h_matrix[(n * (n + 3))/2],
                sdBas->tot_dimSD, initVec, finalVec);

  // Orthogonalize to all previous Lanczo vectors

  for(k = 0; k < n ; k++)  {

    // read lancVec(k) from file into |initVec[]>

    readLancVec(k, sdBas->tot_dimSD, lancProc, initVec);

    // factor = <initVec[] |finalVec[]>

    if(fabs(factor = scalar_product_of_id_vectors(sdBas, initVec, finalVec)) 
                                                    < ZERO_LIMIT) continue;
  
        /*
	** Orthogonalize the new vector to lanczo vector k  
	**       |finalVec[]> <--- |finalVec[]> - factor * |initVec[]>
	*/

    add_a_vector(-factor, sdBas->tot_dimSD, initVec, finalVec);

  } // end of loop k

  // close sectFiles for lancVec: k = -1

  if(k > 0) readLancVec(-1, sdBas->tot_dimSD, lancProc, initVec);

        /*
	** Normalize |finalVec[]> and produce a new lanc_vec. 
	** Calculate the off-diagonal matrix element to the 
	** new Lanczos vector.
	** The new lanc_vec must  have  norm > MATRIX_LIMIT.    
	*/

  if((factor = norm_of_a_vector(sdBas, finalVec)) < MATRIX_LIMIT)
        return FALSE;

  lancProc->h_matrix[((n + 1) * (n + 4))/2 - 1] = factor;
  factor = 1.0 / factor;

  // Normalize the new lanczo vector |final_vec>

  scale_a_vector(sdBas->tot_dimSD, factor, finalVec);

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

void readLancVec(int lancVecNum, int dim, 
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

    sprintf(basFile,"%s%s",lancProc->title, LANCZO_VECTORS);

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
      exit(1);
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
      exit(1);
    }
  }
  readLancVecFromOpenFile(sectFile, &filePtr, lancVecNum, dim, floatMem); 

  // convert float amplitudes to double

  initPtr  = floatMem;
  finalPtr = vec;
  for(k = 0; k < dim; k++)  {
    *(finalPtr++) = (float)(*(initPtr++));
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

void readLancVecFromOpenFile(char *filename, FILE **filePtr,
                    int lancVecNum, int dim, float *floatMem) 
{
  int     n;

  // read vector number
  if(fread((void *)&n,(size_t) sizeof(int), 1, *filePtr) != 1) {
    printf("\nError in function readLancVecFromOpenFile()");
    printf("\nIn reading vector number =  %d", n);
    printf("\nfrom file %s\n",filename);
    exit(1);
  }   // end of if-test

  // test vector number
  if(n !=  lancVecNum) {
    printf("\nError in function readLancVecFromOpenFile()");
    printf("\nWrong vector number read");
    printf("\nread number = %d - should be lancVecNum = %d\n",
                                                n,lancVecNum);
    exit(1);
  }

  // read vector dimension
  if(fread((void *)&n,(size_t) sizeof(int), 1, *filePtr) != 1) {
    printf("\nError in function readLancVecFromOpenFile()");
    printf("\nIn reading vector dimension =  %d", n);
    printf("\nfrom file %s\n",filename);
    exit(1);
  }   // end of if-test

  // test vector dimension
  if(n !=  dim) {
    printf("\nError in function readLancVecFromOpenFile()");
    printf("\nWrong vector dimension read");
    printf("\nread number = %d - should be lancVecNum = %d\n",
                                                n,dim);
    exit(1);
  }

  // read vector float amplitudes
  if(fread((void *)floatMem,(size_t)sizeof(float), (size_t) dim,
                                       *filePtr) != (size_t) dim)  {
    printf("\n\nError in function readLancVecFromOpenFile():");
    printf("\nIn reading %d vector amplitudes from file %s\n",
                                                 dim, filename);
    exit(1);
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

static void saveFinalData(SHELL *model, SP_BAS *spBas,
                     SD_BAS *sdBas, LANC_PROC *lancProc)
{
  char         *func = {"saveFinalData(): "};
  int          dim_h, loop;
  double       *initVec, *finalVec, *j_occ;
  EIGEN        *listEigen;
  EIGEN_DATA   h_eigenData;

     /*
     ** Calculates and stores in lancProc->op_ang[] the complete
     **  set of m-scheme one- and two-particle matrix elements
     **  of angular momentum operator J**2. Calculate and store 
     ** on files diagonal and nondiagonal slater determinant
     ** matrix elements <SD'|J**2|SD>
     */ 

  angularMomentumInteraction(model, spBas, sdBas, &lancProc->sdStoreAng);

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

  j_occ = MALLOC(spBas->numj_orb,double,func,"j_occ[]");

  for(loop = 0; loop < lancProc->states; loop++) {
    h_eigenData.vecNum   = loop; 
    h_eigenData.dim      = sdBas->tot_dimSD;
    h_eigenData.numj_occ = spBas->numj_orb; 
    h_eigenData.eigenVal = listEigen[loop].value;

    // transfer eigenvec with number loop to SD basis

    lancToSD_Bas(lancProc, sdBas->tot_dimSD,initVec,
                       finalVec,dim_h, listEigen[loop].vector);

    // calculate < |J**2| >

/**********************************
   h_eigenData.angMom  = id_lanczo_matrix_element(ANG_INT, spBas,
			           sdBas, &lancProc->sdStoreAng,
                                           finalVec, finalVec);
*********************************/

        /*
         ** do
         **     J**2|finalVec> = |intVec>
         ** where the energy operator J**2 specified through
         ** its matrix elements in m-scheme.
         */

    id_matrix_vector_calc(spBas, sdBas, &lancProc->sdStoreAng,
                                           finalVec, initVec);

    h_eigenData.angMom 
        = scalar_product_of_id_vectors(sdBas, finalVec, initVec);

    if(lancProc->calc_CM == 1) {

      // calculate < |CM| >

      h_eigenData.val_CM = id_lanczo_matrix_element(spBas,
			   sdBas, &lancProc->sdStoreCM, finalVec, finalVec)
                           -1.5;
    }

    // calculate list of single-particle occulation numbers

    id_j_occupation(spBas, sdBas, j_occ, finalVec);

    // save energy eigenvalue and  eigenvector on file

    storeEigenVec(lancProc, h_eigenData, j_occ, finalVec);
  
  }   // end loop:  all eigenvector written to file

  // remove local memory

  for(loop = dim_h - 1; loop >= 0;loop--) {
    free(listEigen[loop].vector); 
  }
  free(listEigen);

} // End: function saveFinaData()

      /*
      ** The function                                        
      **          lanczoToSD_Basis()                       
      ** transform an eigenvector from Lanczos vector basis 
      ** into SD basis and store the result in iteration->final_vec.   
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

    add_a_vector(*(lancAmp++), totDimSD, initVec, finalVec);

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

static double scalar_product_of_id_vectors(SD_BAS *sdBas, double *vec1, double *vec2)
{
  int      loop;
  double   *init_vec, *final_vec, result;

  init_vec  = vec1;                             /* initalization */
  final_vec = vec2;
  result    = D_ZERO;

  if(!sdBas->MJ) {
    for(loop = 0; loop < sdBas->numSD[0]; loop++) {
      result += (*(init_vec++)) * (*(final_vec++));
    }
    result *= 2.0; // end asym contributions
    for(loop = sdBas->numSD[0]; loop < sdBas->tot_dimSD; loop++) {
      result += (*(init_vec++)) * (*(final_vec++));
    } // end sym contributions
  } /* end time reversal symmetry */
  else {
    for(loop = 0; loop < sdBas->tot_dimSD; loop++) {
      result += (*(init_vec++)) * (*(final_vec++));
    }
  }   /* end no time reversal symmetry */

  return result;

} /* End: function scalar_product_of_id_vectors() */

       /*
       ** The function                                       
       **          norm_of_a_vector()                   
       ** calculates the norm of a vector stored in 
       ** vec[] with iteration->tot_dimSD elements. 
       ** The function returns 
       **     result = sqrt( Sum(x_i * x_i))
       ** Possible time-reversal symmetry is included   
       */

static double norm_of_a_vector(SD_BAS *sdBas, double *vec)
{
  int
           loop; 
  double
           *ptr_vec, result;

  result  = D_ZERO;
  ptr_vec = vec;
  if(!sdBas->MJ) {
    for(loop = 0; loop < sdBas->numSD[0]; loop++, ptr_vec++) {
      result += (*ptr_vec) * (*ptr_vec);
    }
    result *= 2.0;   /* contr from time-reversal comp */
    for( loop = sdBas->numSD[0];  loop < sdBas->tot_dimSD; loop++, ptr_vec++) {
      result += (*ptr_vec) * (*ptr_vec);
    }
  }   /* end time reversal symmetry */
  else {
    for(loop = 0; loop < sdBas->tot_dimSD; loop++, ptr_vec++) {
      result += (*ptr_vec) * (*ptr_vec);
    }
  }
  return sqrt(result);

}  /* End: function norm_of_a_vector() */

 /*
  ** The function 
  **      CM_energy_matr_elem()
  ** calculates nondiagonal center_of_Mass matrix 
  ** elements 
  **      <lanc_vec|CM_INTfinal_vec>
  ** for the stored lanc_vec below |final_vec> and
  ** store the result in one-dimesional CM_matrix[].
  */

static void CM_energy_matr_elem(LANC_PROC *lancProc, SP_BAS *spBas, 
           SD_BAS *sdBas, double *initVec, double *finalVec, int n)
{
  int       k;
  
  for(k = 1; k <= n  ; k++)  {

    // read lancVec(k) from file into |initVec[]>

    readLancVec(k-1, sdBas->tot_dimSD, lancProc, initVec);

    lancProc->CM_matrix[((n - 1) * (n + 2))/2 + k] 
       = id_lanczo_matrix_element(spBas,sdBas, 
           &lancProc->sdStoreCM, initVec, finalVec);
  } // end of loop k

  // close sectFiles for lancVec:

  readLancVec(-1, sdBas->tot_dimSD, lancProc, initVec);

}    // CM_energy_matr_elem()

       /*
       ** The function 
       **      angularMomentumInteraction()
       ** calculates and stores on file the complete set of m-scheme
       ** one- and two-particle matrix elements  of the angular
       ** momentum operator J**2. Calculate and store diagonal 
       ** and nondiagonal slater determinant matrix elements 
       ** <SD'|J**2|SD> on files 
       */

static void  angularMomentumInteraction(SHELL *model, SP_BAS *spBas,
                                 SD_BAS *sdBas, SD_STORE *sdStore)
{
  // two-particle angular momentum matrix elements

  id_ang_mom_interaction(spBas, &sdStore->op);

  // add singel-particle angular momentum matrix element

  add_id_single_particle_ang(sdBas->part, spBas, &sdStore->op);

  // sdStore_ANG data initialization

  sprintf(sdStore->result,"%s%s", model->title, RESULT_OUTPUT); 
  sprintf(sdStore->title,"%s%s" , model->title, SD_STORE_ANG); 

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

  // calculate and store <SD'|ANG|SD> 

  sdStore->typeInt = ANG_INT;

  storeSD_MatrElemInMem(spBas, sdBas, sdStore);
  
  storeSD_MatrElemOnFile(spBas, sdBas, sdStore);


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
    **     storeEigenVec()
    ** opens a sect file and stores a lanczos 
    ** eigenvector,If current sect file has 
    ** reacded MAX_FILE_SIZE current sect file 
    ** is closed and a new is opened
    */

void storeEigenVec(LANC_PROC *lancProc, EIGEN_DATA h_eigenData,
                                      double *j_occ,double *vec)
{
  static char   basFile[ONE_LINE], sectFile[ONE_LINE];
  static int    maxLancVecPerFile, sectFileNo;
  FILE          *filePtr;

  if(h_eigenData.vecNum == 0) { // initialization

    sprintf(basFile,"%s%s",lancProc->title, EIGEN_VECTORS);

    // maximum eigenvectors in a single file

    maxLancVecPerFile 
          = (  (MAX_FILE_SIZE * M_BYTES)
	       / (h_eigenData.dim * sizeof(float) + sizeof(EIGEN_DATA)));
    sectFileNo = 1;  // initialization
    sprintf(sectFile,"%s%d", basFile, sectFileNo);

    // create a new section file to store lanczo eigenvectors

    if((filePtr = fopen(sectFile,"wb")) == NULL) {
      printf("\n\nError in function storeEigenVec():");
      printf("\nNot allowed to create the file %s\n", sectFile);
      exit(1);
    }
    fclose(filePtr);
  }
  if((h_eigenData.vecNum/maxLancVecPerFile) >= sectFileNo) {

    // create a new section file

    sectFileNo++;
    sprintf(sectFile,"%s%d", basFile, sectFileNo);
    if((filePtr = fopen(sectFile,"wb")) == NULL)  { 
      printf("\n\nError in function storeEigenVec():");
      printf("\nNot allowed to create the file %s\n", sectFile);
      exit(1);
    }
    fclose(filePtr);
  }
  appendEigenVecToFile(sectFile, h_eigenData, j_occ, vec); 

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

void appendEigenVecToFile(char *filename, EIGEN_DATA h_eigenData,
                                     double *j_occ,double *vector) 
{
  char    *func = {" appendEigenVecToFile(): "};
  int     k;
  float   *tempMem, *finalPtr;
  double  *initPtr;
  FILE    *filePtr;

  if((filePtr = fopen(filename,"ab+")) == NULL) {
    printf("\nError in function appendEigenVecToFile();");
    printf("\nWrong file = %s to store a vector no %d\n",
                              filename, h_eigenData.vecNum);
    exit(1);
  }
  // struct h_eigenData
  if(fwrite((const void *)&h_eigenData,(size_t) sizeof(EIGEN_DATA), 1, filePtr) != 1) {
    printf("\nError in function appendEigenVecToFile()");
    printf("\nIn writing structure h_eigenData");
    printf("\nto file %s\n",filename);
    exit(1);
  }   // end of if-test

  // list of single-j-occupation 

  if(fwrite((const void *)j_occ,(size_t) sizeof(double), h_eigenData.numj_occ,
                                             filePtr) != h_eigenData.numj_occ) {
    printf("\nError in function appendEigenVecToFile()");
    printf("\nIn writing single-particle j_occ[]");
    printf("\nto file %s\n",filename);
    exit(1);
  }   // end of if-test

  tempMem = MALLOC(h_eigenData.dim, float, func, "tempMem");
  initPtr = vector;
  finalPtr = tempMem;
  for(k = 0; k < h_eigenData.dim; k++)  {
    *(finalPtr++) = (float)(*(initPtr++));
  }
  if(fwrite((const void *)tempMem,(size_t) sizeof(float), 
                (size_t) h_eigenData.dim, filePtr) != h_eigenData.dim) { 
    printf("\n\nError in function appendEigenVecToFile():");
    printf("\nIn writing eigenvector no %d diag to file %s\n\n",
                               h_eigenData.dim, filename);
    exit(1);
  }
  free(tempMem);
 
  fclose(filePtr);

} // End: function appendEu\igenVecToFile()

   /*
   ** The function                            
   **        id_j_occupation()                   
   ** calculate the j-occupation numbers  of an identicle particle shell model
   ** eigenstate from the lanczo procedure and returned the result in j_occ[]. 
   */

static void id_j_occupation(SP_BAS *spBas, SD_BAS *sdBas, 
			   double *j_occ, double *eigenVec)
{
  int       orb, loop, part;
  ULL       pos, *bas_ptr;
  double    probability, *ampl_ptr;

  for(loop = 0; loop < spBas->numj_orb; loop++) j_occ[loop] = 0.0;

  if(sdBas->MJ == 0) {   // time-reversal symmetry
    bas_ptr  = sdBas->SD;
    ampl_ptr = eigenVec;
    loop     = 0;
    if(sdBas->numSD[0]) {  // contribution from |asymSD()>
      for(loop = 0; loop < sdBas->numSD[0]; ampl_ptr++, 
                                            bas_ptr++, loop++)   {
	probability = 2.0 * (*ampl_ptr) * (*ampl_ptr);
	pos    = ULL_ONE;  // orbital initialization
	orb    = 0;
	part   = 0;
	do   {
	  for( ;!((*bas_ptr) & pos); orb++, pos <<= 1); // occ. orbit found
	  j_occ[spBas->mbas[orb].orb] += probability;	  orb++;
	  pos <<= 1; 
	} while(++part < sdBas->part);      
      } // loop through all |asymSD()>
    } // end contribution from asymmetric contribution
    if(sdBas->numSD[1]) {              // contribution from |symSD()>
      for( ; loop < sdBas->tot_dimSD; ampl_ptr++, bas_ptr++, loop++)   {
	probability = (*ampl_ptr) * (*ampl_ptr);
	pos    = ULL_ONE;  // orbital initialization
	orb    = 0;
	part   = 0;
	do   {
	  for( ;!((*bas_ptr) & pos); orb++, pos <<= 1); // occ.orbit found
	  j_occ[spBas->mbas[orb].orb] += probability;
	  orb++;
	  pos <<= 1; 
	} while(++part < sdBas->part);      
      } // loop through all |symSD()>
    } // end contribution frm symmetric terms
  } // end time-reversal case
  else  {  // contribution from |nosymSD()>
    bas_ptr  = sdBas->SD; // initialization
    ampl_ptr = eigenVec;
    for(loop = 0; loop < sdBas->numSD[0]; ampl_ptr++,
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
      } while(++part < sdBas->part);      
    } // loop through all |nosymSD()>
  } // end contribution frm nonsymmetric case
} // End: function id_j_occupation()
