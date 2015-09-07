



//***************  The module PAR-id-lanc-result.c ******************  

#include "PAR-shell-model.h"

     /*
     ** The entrance function
     **         id_lanczos_result()
     ** takes the calculated result from identical particle
     ** Lanczos process and write it to output file
     */

	    /**** local function declarations ****/

void readEigenData(int eigenVecNum, int dim, LANC_PROC *lancProc,
		   EIGEN_DATA *h_eigenData, int numj_occ, double *j_occ);
    /*
    ** opens the file "filename" and read 
    **    h_eigenData.vecNum   - local vector number
    **                           on current file
    **               .dim       - dim of vector[]
    **               .num_j_occ - num j-orbits
    **               .eigenVal  - eigenvalue amplitudes
    **               .angMom    - < |J**2|>
    **               .val_CM    - < |CM| > if val_CM = 1
    ** Note that the first vector on file is numbered 
    ** vecNum = ZERO.
    ** The amplitudes  is transferred 
    ** from float to double 
    */ 

               /**** The function definitions  ****/ 
 
     /*
     ** The entrance function
     **         id_Lanczos_result()
     ** takes the calculated result from identical particle
     ** Lanczos process and write it to output file
     */

void id_lanczos_result(SP_BAS *spBas,SD_BAS *sdBas,LANC_PROC *lancProc)
{
  char  *func = {"id_lanczos_result()"};
  int           loop, loop_1; 
  FILE          *filePtr;
  double        *j_occ, E_0;
  EIGEN_DATA    h_eigenData;

  // Start the output procedure

  if((filePtr = fopen(lancProc->title,"a")) == NULL)   {
    printf("\n\nRank(%d): Error in function id_Lanczos_result():", Rank);
    printf("\nWrong file = %s to open the output data file\n",lancProc->title);
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

    fprintf(filePtr,"\n\nE(%d)= %7.4f", loop, h_eigenData.eigenVal- E_0);
    fprintf(filePtr," <J**2> = %7.4f", h_eigenData.angMom);

   // print orbital j-occupation

    for(loop_1 = 0; loop_1 < spBas->numj_orb; loop_1++) {
      if( ((spBas->jbas[loop_1].l)/2)*2 == spBas->jbas[loop_1].l) {
	fprintf(filePtr," %2d/2+",spBas->jbas[loop_1].j);
      }
      else {
	fprintf(filePtr," %2d/2-",spBas->jbas[loop_1].j);
      }
    }
    fprintf(filePtr,"\n                 N(j):          ");
    for(loop_1 = 0; loop_1 < spBas->numj_orb; loop_1++) {
      fprintf(filePtr,"%5.3f ", j_occ[loop_1]);
    }

   // print center_of_Mass value < |CM| >

    if(lancProc->calc_CM == 1) {
      fprintf(filePtr,"\n  Center_of_Mass: < |H(CM)| > = %7.4f",
                                             h_eigenData.val_CM);
    }

 }  // All eigenvalues/eigenvectors  printed 

  fclose(filePtr);

  free(j_occ);   // remove local memory 

  // remove all files containing <SD'|J**2|SD>

  //  /*******  Remove all scratch......   ************

}  // End: function id_Lanczos_result()

    /*
    ** The function 
    **     readEigenData()
    ** opens the file "filename" and read 
    **    h_eigenData.vecNum   - local vector number
    **                           on current file
    **               .dim       - dim of vector[]
    **               .eigenVal  - eigenvalue amplitudes
    **               .angMom    - < |J**2|>
    **               .val_CM    - < |CM| > if val_CM = 1
    **               .num_j_occ - num j-orbits
    ** Note that the first vector on file is numbered 
    ** vecNum = ZERO.
    ** The amplitudes  is transferred 
    ** from float to double 
    */ 

void readEigenData(int eigenVecNum, int dim, LANC_PROC *lancProc,
	      EIGEN_DATA *h_eigenData, int numj_occ, double *j_occ)
{
  char   filename[ONE_LINE];
  int    k;
  FILE   *filePtr;

  sprintf(filename,"%s%s%d%s%d",lancProc->scratchFile,"Rank",
	                                Rank,EIGEN_VECTORS,1);
  if((filePtr = fopen(filename,"rb")) == NULL) {
    printf("\n\nError in function readEigenDat():");
    printf("\nNot allowed to open file %s\n",filename);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }
 
  k = eigenVecNum;

  fseek(filePtr,(long)(k*(  sizeof(EIGEN_DATA)
                          + numj_occ * sizeof(double)
			  + dim*sizeof(float))),
	                  SEEK_SET);

  // read structure EIGEN_DATA h_eigenData

  if(fread((void *)h_eigenData,(size_t) sizeof(EIGEN_DATA),1, 
                                               filePtr) != 1) {
    printf("\nError in function readEigenData()");
    printf("\nIn reading h_eigendata");
    printf("\nfrom file %s\n",filename);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }   // end of if-test

  // test vector number

  if(h_eigenData->vecNum != eigenVecNum) {
    printf("\nError in function readEigenData()");
    printf("\nWrong vector number read");
    printf("\nread number = %d - should be lancVecNum = %d\n",
                              h_eigenData->vecNum,eigenVecNum);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }

  // read single-particle occupation numbers to j_occ[]

  if(fread((void *)j_occ,(size_t) sizeof(double), h_eigenData->numj_occ, 
                                    filePtr) != h_eigenData->numj_occ) {
    printf("\nError in function readEigenData()");
    printf("\nIn reading %d j_occ[]",h_eigenData->numj_occ);
    printf("\nfrom file %s\n",filename);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }   // end of if-test

} // End: function  readEigenData()
