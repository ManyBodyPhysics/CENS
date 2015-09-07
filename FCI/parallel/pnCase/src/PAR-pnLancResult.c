
//***************  The module PAR-id-lanc-result.c ******************  

#include "PAR-pnShellModel.h"

     /*
     ** The entrance function
     **         PNLancResult()
     ** takes the calculated result from identical particle
     ** Lanczos process and write it to output file
     */

	    /**** local function declarations ****/

static void readEigenData(int eigenVecNum, GR_BAS *grbas,EIGEN_DATA *h_eigenData,
		   double *protonj_orb, double *neutronj_orb);
   /*
   ** opens a sect file and reads a lanczos 
   ** eigenvector.
   ** If current sect file has reacded 
   ** MAX_FILE_SIZE current sect file is
   ** closed and a new is opened
   */

static void readEigenDataFromOpenFile(char *filename, FILE *filePtr,int eigenVecNum, 
	   EIGEN_DATA *h_eigenData, double *protonj_orb,double *neutronj_orb);
   /*
   ** opens the file "file_name" and read 
   **    h_eigenData.vecNum       - local vector number
   **                               on current file
   **               .dim          - dim of vector[]
   **               .protonj_orb  - num j-orbits
   **                neutronj_orb - num j-orbits
   **               .eigenVal  - eigenvalue amplitudes
   **               .angMom    - < |J**2|>
   **               .num_j_occ - num j-orbits
   ** Note that the first vector on file is numbered 
   ** vecNum = ZERO.
   */ 
               /**** End: function declarations ****/

               /**** The function definitions  ****/ 
 
     /*
     ** The entrace function
     **         pnLancResult()
     ** takes the calculated result from pn particle 
     ** Lanczos process and write it to output file
     */

void pnLancResult(GR_BAS *grBas, LANC_PROC *lancProc)
{
  char          *func = {"pnLancResult(): "}, filename[ONE_LINE];
  int           loop, loop_1; 
  FILE          *filePtr;
  double        *protonj_occ, *neutronj_occ,E_0;
  EIGEN_DATA    h_eigenData;

  // Start the output procedure

  sprintf(filename,"%s%s",grBas->title,RESULT_OUTPUT);

  if((filePtr = fopen(filename,"a")) == NULL)   {
    printf("\n\nRank%d: Error in function pnLancResult():",Rank);
    printf("\nWrong file = %s to open the output data file\n",filename);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }

  switch(lancProc->run_code) {
    case 0:
      fprintf(filePtr,"\n\nThe Lanczos iteration process has reached ");
      fprintf(filePtr,"\nmaximum iterations = %d - specified in input data\n",
                                                     lancProc->num_iterations);
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
       MPI_Abort(MPI_COMM_WORLD,Rank);
  } // end switch()

  fprintf(filePtr,"\nFinal Eigenvalues");

  // The total parity of the eigenstates

  if(grBas->P == +1)
    fprintf(filePtr,"\nThe total parity is positive\n");
  else
    fprintf(filePtr,"\nThe total parity is negative\n");

  // Energy eigenstates - open eigenvector files

  E_0 = 0;    // initialization

  protonj_occ  = MALLOC(grBas->spBas[0]->numj_orb,double, func, "proton_j[]");
  neutronj_occ = MALLOC(grBas->spBas[1]->numj_orb,double, func, "neutron_j[]");

  for(loop = 0; loop < lancProc->states; loop++)    {
    readEigenData(loop,grBas,&h_eigenData,protonj_occ, neutronj_occ);
    if(loop == 0) {
      E_0 = h_eigenData.eigenVal;    // ground state energy
      fprintf(filePtr,"\n\nThe ground state energy = %9.4f",E_0);
    }

    // print energy and angular momentum 

    fprintf(filePtr,"\n\nE(%d)= %7.4f", loop, h_eigenData.eigenVal- E_0);
    fprintf(filePtr," <J**2> = %7.4f", h_eigenData.angMom);

    if(grBas->calc_CM == 1) {
      fprintf(filePtr,"   <R**2> = %7.4f", h_eigenData.valCM);
    }

   // print proton orbital j-occupation

    fprintf(filePtr,"\nProton single-part.  distrib");

    for(loop_1 = 0; loop_1 < grBas->spBas[PROTON]->numj_orb; loop_1++) {
      if( ((grBas->spBas[PROTON]->jbas[loop_1].l)/2)*2 
	  == grBas->spBas[PROTON]->jbas[loop_1].l)      {
	fprintf(filePtr," %2d/2+",grBas->spBas[PROTON]->jbas[loop_1].j);
      }
      else {
	fprintf(filePtr," %2d/2-",grBas->spBas[PROTON]->jbas[loop_1].j);
      }
    }
    fprintf(filePtr,"\n                     N(j)    ");

    for(loop_1 = 0; loop_1 < grBas->spBas[PROTON]->numj_orb; loop_1++) {
      fprintf(filePtr,"%5.3f ", protonj_occ[loop_1]);
    }

   // print neutron orbital j-occupation

    fprintf(filePtr,"\nNeutron single-part. distrib");

    for(loop_1 = 0; loop_1 < grBas->spBas[NEUTRON]->numj_orb; loop_1++) {
      if( ((grBas->spBas[NEUTRON]->jbas[loop_1].l)/2)*2 
	  == grBas->spBas[NEUTRON]->jbas[loop_1].l)      {
	fprintf(filePtr," %2d/2+",grBas->spBas[NEUTRON]->jbas[loop_1].j);
      }
      else {
	fprintf(filePtr," %2d/2-",grBas->spBas[NEUTRON]->jbas[loop_1].j);
      }
    }
    fprintf(filePtr,"\n                     N(j)    ");

    for(loop_1 = 0; loop_1 < grBas->spBas[NEUTRON]->numj_orb; loop_1++) {
      fprintf(filePtr,"%5.3f ", neutronj_occ[loop_1]);
    }

   // print center_of_Mass value < |CM| >

    if(grBas->calc_CM == 1) {
      fprintf(filePtr,"\n  Center_of_Mass: < |H(CM)| > = %7.4f",
                                             h_eigenData.valCM);
    }

}  // All eigenvalues/eigenvectors  printed 

  fclose(filePtr);

  free(neutronj_occ);   // remove local memory 
  free(protonj_occ);

  // remove all files containing <SD'|J**2|SD>

/***************************
  sprintf(filename,"%s%s",grBas->title,SD_STORE_ANG,DIAG);
  remove(filename);

  for(loop = 1; ;) {
    sprintf(filename,"%s%s%s%d",lancProc->title, SD_STORE_ANG,NONDIAG, loop++);
    if(remove(filename)) break;
  }
********************************/

}  // End: function pnLancResult()

    /*
    ** The function 
    **     readEigenData()
    ** opens a sect file and reads an 
    ** eigenvector.
    ** If current sect file has reached 
    ** MAX_FILE_SIZE current sect file is
    ** closed and a new is opened
    */

void readEigenData(int eigenVecNum, GR_BAS *grBas,EIGEN_DATA *h_eigenData,
		   double *protonj_occ, double *neutronj_occ)
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
	     / (grBas->totNumSD_ZN_oscLim * sizeof(float)+ sizeof(EIGEN_DATA)));
    sectEigenVecNum = 0;
    sectFileNo = 1;  // initialization
    sprintf(sectFile,"%s%s%d",grBas->title,EIGEN_VECTORS,sectFileNo);

    // open first eigenvector sectFile

    if((filePtrBas = fopen(sectFile,"rb")) == NULL) {
      printf("\n\nRank%d: Error in function readEigenData():",Rank);
      printf("\nNot allowed to open file %s\n", sectFile);
      MPI_Abort(MPI_COMM_WORLD,Rank);
    }
  }
  if((eigenVecNum/maxEigenVecPerFile) >= sectFileNo) {

    // close current and open next eigenvector sectFile

    fclose(filePtrBas);

    sectEigenVecNum += sectFileNo * maxEigenVecPerFile; 

    sectFileNo++;  // initialization
    sprintf(sectFile,"%s%s%d",grBas->title,EIGEN_VECTORS,
                                                sectFileNo);
    if((filePtrBas = fopen(sectFile,"rb")) == NULL)  { 
      printf("\n\nRank%d: Error in function readEigenData()",Rank);
      printf("\nNot allowed to open file %s\n", sectFile);
      MPI_Abort(MPI_COMM_WORLD,Rank);
    }
  }
  k       = eigenVecNum - sectEigenVecNum * maxEigenVecPerFile; 
  filePtr = filePtrBas;
  fseek(filePtr,(long)(k*(  sizeof(EIGEN_DATA)
                          + grBas->spBas[PROTON]->numj_orb  * sizeof(double)
                          + grBas->spBas[NEUTRON]->numj_orb * sizeof(double)
			  + h_eigenData->dim*sizeof(float))),
	                  SEEK_SET);

  readEigenDataFromOpenFile(sectFile,filePtr, eigenVecNum,
                            h_eigenData, protonj_occ, neutronj_occ); 

} // End: function readEigenData()

   /*
   ** The function                                           
   **            readEigenDataFromOpenFile( )                    
   ** opens the file "file_name" and read 
   **    h_eigenData.vecNum       - local vector number
   **                               on current file
   **               .dim          - dim of vector[]
   **               .protonj_orb  - num j-orbits
   **                neutronj_orb - num j-orbits
   **               .eigenVal  - eigenvalue amplitudes
   **               .angMom    - < |J**2|>
   **               .num_j_occ - num j-orbits
   ** Note that the first vector on file is numbered 
   ** vecNum = ZERO.
   */ 

void readEigenDataFromOpenFile(char *filename, FILE *filePtr,int eigenVecNum, 
              EIGEN_DATA *h_eigenData, double *protonj_occ,double *neutronj_occ) 
{
  // read structure EIGEN_DATA h_eigenData

  if(fread((void *)h_eigenData,(size_t) sizeof(EIGEN_DATA), 1, 
	                                filePtr) != (size_t)1) {
    printf("\n\nRank%d: Error in function readEigenDataFromOpenFile()",Rank);
    printf("\nIn reading h_eigendata");
    printf("\nfrom file %s\n",filename);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }   // end of if-test

  // test vector number

  if(h_eigenData->vecNum !=  eigenVecNum) {
    printf("\n\nRank%d: Error in function readEigenDataFromOpenFile()",Rank);
    printf("\nWrong vector number read");
    printf("\nread number = %d - should be lancVecNum = %d\n",
                             h_eigenData->vecNum,eigenVecNum);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }

  // read single-particle occupation numbers

  if(fread((void *)protonj_occ,(size_t) sizeof(double), 
                      h_eigenData->numj_occ[0],filePtr) 
                      != (size_t)h_eigenData->numj_occ[0]) {
    printf("\n\nRank%d: Error in function readEigenDataFromOpenFile()",Rank);
    printf("\nIn reading %d j_occ[]",h_eigenData->numj_occ[PROTON]);
    printf("\nfrom file %s\n",filename);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }
  if(fread((void *)neutronj_occ,(size_t) sizeof(double), 
                      h_eigenData->numj_occ[1], filePtr) 
                   != (size_t)h_eigenData->numj_occ[1]) {
    printf("\n\nRank%d: Error in function readEigenDataFromOpenFile()",Rank);
    printf("\nIn reading %d j_occ[]",h_eigenData->numj_occ[NEUTRON]);
    printf("\nfrom file %s\n",filename);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }   // end of if-test

} // End: function  readEigenDataFromOpenFile()
