
//***************  The module id-lanc-result.c ******************  

#include "shell-model.h"

     /*
     ** The entrance function
     **         id_lanczos_result()
     ** takes the calculated result from identical particle
     ** Lanczos process and write it to output file
     */

	    /**** local function declarations ****/

void readEigenData(int eigenVecNum, int dim, LANC_PROC *lancProc,
	      EIGEN_DATA *h_eigenData, int numj_orb, double *j_occ);
   /*
   ** opens a sect file and reads a lanczos 
   ** eigenvector.
   ** If current sect file has reacded 
   ** MAX_FILE_SIZE current sect file is
   ** closed and a new is opened
   */

void readEigenDataFromOpenFile(char *filename, FILE *filePtr, int eigenVecNum,
                                       EIGEN_DATA *h_eigenData, double *j_occ); 
   /*
   ** opens the file "file_name" and read 
   **    h_eigenData.vecNum   - local vector number
   **                           on current file
   **               .dim      - dim of vector[]
   **               .eigenVal - eigenvalue amplitudes
   **               .angMom   - < |J**2|>
   **               .val_CM   - < |CM| > if val_CM = 1
   **               .numj_occ - number of j-orbits
   ** Note that the first vector on file is numbered 
   ** vecNum = ZERO. The amplitudes  is transferred 
   ** to type float and stored on file
   */ 

               /**** End: function declarations ****/

               /**** The function definitions  ****/ 
 
     /*
     ** The entrace function
     **         id_Lanczos_result()
     ** takes the calculated result from identical particle
     ** Lanczos process and write it to output file
     */

void id_lanczos_result(SP_BAS *spBas, SD_BAS *sdBas, LANC_PROC *lancProc)
{
  char          *func = {"id_lanczos_result(): "}, filename[ONE_LINE];
  int           loop, loop_1; 
  FILE          *filePtr;
  double        *j_occ, E_0;
  EIGEN_DATA    h_eigenData;

  // Start the output procedure

  sprintf(filename,"%s%s", lancProc->title,RESULT_OUTPUT);

  if((filePtr = fopen(filename,"a")) == NULL)   {
    printf("\n\nError in function id_Lanczos_result():");
    printf("\nWrong file = %s to open the output data file\n",filename);
    exit(1);
  }
  switch(lancProc->run_code) {
    case 0:
      fprintf(filePtr,"\n\nThe Lanczos iteration process has reached ");
      fprintf(filePtr," maximum iterartions = %d - specified in input data\n",
                                                     lancProc->max_iterations);
       break;
    case 1:  
      fprintf(filePtr,"\n\nThe Lanczos iteration process has converged ");
      fprintf(filePtr," after %d iterations.\n", lancProc->num_iterations);
      break;
    case 2:  
      fprintf(filePtr,"\n\nNo more Lanczos vectors");
      fprintf(filePtr," after %d iterations.", lancProc->num_iterations);
      fprintf(filePtr, "\nNo energy convergence\n");
               break;
    default:
       fprintf(filePtr,"\n\nAbnormal termination of the Lanczos iteration");
       fprintf(filePtr," process - run_code = %d\n", lancProc->run_code);
       exit(1);
  } // end switch()

  fprintf(filePtr,"\nFinal Eigenvalues");

  // The total parity of the eigenstates

  if(sdBas->P == +1)
    fprintf(filePtr,"\nThe total parity is positive\n");
  else
    fprintf(filePtr,"\nThe total parity is negative\n");


  // Energy eigenstates - open eigenvector files

  E_0 = 0;    // initialization

  j_occ = MALLOC(spBas->numj_orb,double, func, "j-occ[]");

  for(loop = 0; loop < lancProc->states; loop++)    {

    readEigenData(loop, sdBas->tot_dimSD, lancProc,
		  &h_eigenData, spBas->numj_orb, j_occ);

    if(loop == 0) {
      E_0 = h_eigenData.eigenVal;    // ground state energy
      fprintf(filePtr,"\n\nThe ground state energy = %9.4f",E_0);
    }

    // print energy and angular momentum 

    fprintf(filePtr,"\n\nE(%d)= %9.4f", loop, h_eigenData.eigenVal- E_0);
    fprintf(filePtr,"  <J**2> = %7.4f", h_eigenData.angMom);

   // print orbital j-occupation
    
    for(loop_1 = 0; loop_1 < spBas->numj_orb; loop_1++) {
      if( ((spBas->jbas[loop_1].l)/2)*2 == spBas->jbas[loop_1].l) {
	fprintf(filePtr,"  %2d/2(+) ",spBas->jbas[loop_1].j);
      }
      else {
	fprintf(filePtr,"  %2d/2(-) ",spBas->jbas[loop_1].j);
      }
    }
    fprintf(filePtr,"\n                    N(j):          ");
    for(loop_1 = 0; loop_1 < spBas->numj_orb; loop_1++) {
      fprintf(filePtr,"%6.3f    ", j_occ[loop_1]);
    }

   // print center_of_Mass value < |CM| >

    if(lancProc->calc_CM == 1) {
      fprintf(filePtr,"\n  Center_of_Mass: < |H(CM)| > = %7.4f",
                                             h_eigenData.val_CM);
    }

 }  // All eigenvalues/eigenvectors  printed 


  // close eigenvector file

  readEigenData(-1, sdBas->tot_dimSD, lancProc, &h_eigenData,
                                     spBas->numj_orb, j_occ);

  fclose(filePtr);  // close output file

  free(j_occ);   // remove local memory 

           // remove all files containing <SD'|J**2|SD>

  sprintf(filename,"%s%s%s",lancProc->title, SD_STORE_ANG,DIAG);
  remove(filename);

  for(loop = 1; ;) {
    sprintf(filename,"%s%s%s%d",lancProc->title, SD_STORE_ANG,NONDIAG, loop++);
    if(remove(filename)) break;
  }

}  // End: function id_Lanczos_result()

    /*
    ** The function 
    **     readEigenData()
    ** opens a sect file and reads a lanczo 
    ** based eigen data.
    ** If current sect file has reacded 
    ** MAX_FILE_SIZE current sect file is
    ** closed and a new is opened
    */

void readEigenData(int eigenVecNum, int dim, LANC_PROC *lancProc,
               EIGEN_DATA *h_eigenData, int numj_orb, double *j_occ)
{
  static char   sectFile[ONE_LINE];
  static int    maxEigenVecPerFile, sectFileNo, sectEigenVecNum;
  static FILE   *filePtrBas;
  int           k;
  FILE          *filePtr;

  if(eigenVecNum < 0) {
    fclose(filePtrBas); // close eigenvector file and return
    return;
  }
  if(eigenVecNum == 0) { // initialization

    // maximum eigenvectors in a single file

    maxEigenVecPerFile 
          = (  (MAX_FILE_SIZE * M_BYTES)
	     / (dim * sizeof(float)+ sizeof(EIGEN_DATA)));
    sectEigenVecNum = 0;
    sectFileNo = 1;  // initialization
    sprintf(sectFile,"%s%s%d",lancProc->title,EIGEN_VECTORS,
                                                sectFileNo);

    // open first eigenvector sectFile

    if((filePtrBas = fopen(sectFile,"rb")) == NULL) {
      printf("\n\nError in function readEigenData():");
      printf("\nNot allowed to open file %s\n", sectFile);
      exit(1);
    }
  }
  if((eigenVecNum/maxEigenVecPerFile) >= sectFileNo) {

    // close current and open next eigenvector sectFile

    fclose(filePtrBas);

    sectEigenVecNum += sectFileNo * maxEigenVecPerFile; 

    sectFileNo++;  // initialization
    sprintf(sectFile,"%s%s%d",lancProc->title,EIGEN_VECTORS,
                                                sectFileNo);
    if((filePtrBas = fopen(sectFile,"rb")) == NULL)  { 
      printf("\n\nError in function readEigenData():");
      printf("\nNot allowed to open file %s\n", sectFile);
      exit(1);
    }
  }
  k       = eigenVecNum - sectEigenVecNum * maxEigenVecPerFile; 
  filePtr = filePtrBas;

  fseek(filePtr,(long)(k*(h_eigenData->dim*sizeof(float)
			  + sizeof(EIGEN_DATA)
		          + numj_orb*sizeof(double))), 
                            SEEK_SET);

  readEigenDataFromOpenFile(sectFile, filePtr, eigenVecNum, 
				  h_eigenData, j_occ); 

} // End: function readEigenData()

   /*
   ** The function                                           
   **            readEigenDataFromOpenFile( )                    
   ** opens the file "file_name" and read 
   **    h_eigenData.vecNum   - local vector number
   **                           on current file
   **               .dim      - dim of vector[]
   **               .eigenVal - eigenvalue amplitudes
   **               .angMom   - < |J**2|>
   **               .val_CM   - < |CM| > if val_CM = 1
   **               .numj_occ - number of j-orbits
   ** Note that the first vector on file is numbered 
   ** vecNum = ZERO. The amplitudes  is transferred 
   ** to type float and stored on file
   */ 

void readEigenDataFromOpenFile(char *filename, FILE *filePtr, int eigenVecNum, 
                                        EIGEN_DATA *h_eigenData, double *j_occ) 
{
  // read structure EIGEN_DATA h_eigenData

  if(fread((void *)h_eigenData,(size_t) sizeof(EIGEN_DATA), 1, 
                                                 filePtr) != 1) {
    printf("\nError in function readEigenDataFromOpenFile()");
    printf("\nIn reading h_eigendata");
    printf("\nfrom file %s\n",filename);
    exit(1);
  }   // end of if-test

  // test vector number

  if(h_eigenData->vecNum !=  eigenVecNum) {
    printf("\nError in function readEigenDataFromOpenFile()");
    printf("\nWrong vector number read");
    printf("\nread number = %d - should be lancVecNum = %d\n",
                             h_eigenData->vecNum,eigenVecNum);
    exit(1);
  }

  // read j_occupation amplitude

  if(fread((void *)j_occ,(size_t)sizeof(double), (size_t) h_eigenData->numj_occ,
                                 filePtr) != (size_t) h_eigenData->numj_occ)  {
    printf("\n\nError in function readEigenDataFromOpenFile():");
    printf("\nIn reading %d j_occ[] amplitudes from file %s\n",
                                 h_eigenData->numj_occ, filename);
    exit(1);
  }
} // End: function  readEigenDataFromOpenFile()
