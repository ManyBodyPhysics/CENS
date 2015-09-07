/****************  The module  pnMatrVectorCalc.c  ******************/

#include  "PAR-pnShellModel.h"

   /*
   ** The entrance function
   **         pnMatrixVectorCalc()     
   ** takes as input an orthonormalized Lanczos
   ** vector and performs the process 
   **    H(OP)*|iniVec> --> |finalVec> 
   */

               /****  local function declarations ****/

static void pnDiagMemSD_Contr(int typeCalc, GR_BAS *grBas, GROUP *group,
			 double *initVec, double *finalVec);
     /*
     **       H(OP) * |initVec> ---> |finalVec>
     ** from the precalculated diagonal matrix elements
     **  stored as float values in memory
     */

static void pnNondiagMemSD_Contr(int typeCalc, GR_BAS *grBas, GROUP *group,
		       	   double *initVec, double *finalVec);
     /*
     ** calculates contribution to 
     **       H(op) * |initVec> ---> |initVec>
     ** from the precalculated nondiagonal matrix elements
     **  stored as float values in memory
     */

static ULL countBits(BARR_ELTYPE *bitArray, ULL size); 
     /*
     ** counts and return number
     ** of bits set in a bitArray[]
     */

               /****   end local function declarations ****/



               /**** The function definitions  ****/ 
 
   /*
   ** The entrance function
   **         pnMatrixVectorCalc()     
   ** takes as input an orthonormalized Lanczos
   ** vector and performs the process 
   **    H(OP)*|iniVec> --> |finalVec> 
   */

void  pnMatrixVectorCalc(int typeCalc, GR_BAS *grBas, GROUP *group,
		       double *initVec,double *finalVec)
{
  // contribution to |finalVec> from disgMatrElem

  pnDiagMemSD_Contr(typeCalc,grBas,group,initVec,finalVec);


  // contribution to |finalVec> from nondiagMatrElem


  if(grBas->permMemSize  > 0) {

    pnNondiagMemSD_Contr(typeCalc,grBas,group,initVec,finalVec);
  }

  //  explicit calculation of <SD'|VEFF|SD>

  /************************************
  nondiagExplicitCalc(nextInitSD, sdStore, spBas, sdBas,
                                    initVec, finalVec);
  **************************/

  // end contribution explicit calculation

} // End: function   pnMatrixVectorCalc()

     /*
     **        pnDiagMemSD_Contr()
     ** calculates contribution to 
     **       H(OP) * |initVec> ---> |finalVec>
     ** from the precalculated diagonal matrix elements
     **  stored as float values in memory
     */

static void pnDiagMemSD_Contr(int typeCalc, GR_BAS *grBas, GROUP *group,
                             double *initVec, double *finalVec)
{
  int       grNum, numAmplSteps,pStep, nStep, numSD_Z, numSD_N;
  int      numInit,numFinal;
  float     *diagPtr;
  double    *initPtr, *finalPtr;

  diagPtr = (typeCalc == ANG_CALC) ? grBas->diagMemElemANG
                                           : grBas->diagMemElem;
  for(grNum = 0; grNum < grBas->num_gr; grNum++) {

    if(group[grNum].RankNumSD_ZN_oscLim == 0) {
      continue; // no basis states in current group
    }
    initPtr      = initVec  + group[grNum].startAmp;
    finalPtr     = finalVec + group[grNum].startAmp;
    numAmplSteps = group[grNum].numAmplStep;
    numSD_Z      = group[grNum].numSD[PROTON];
    numSD_N      = group[grNum].numSD[NEUTRON];

    for(pStep = 0; pStep < numSD_Z;initPtr += numSD_N,finalPtr+= numSD_N,
	  pStep++) {
      for(nStep = Rank; nStep < numSD_N; nStep += numAmplSteps) {

	finalPtr[nStep] = initPtr[nStep] * ((double)(*diagPtr));
	diagPtr++;

      }  // nStep

    }  // pStep

  } // end groupNum

} // end function  pnDiagMemSD_Contr()

     /*
     ** The function
     **        pnNondiagMemSD_Contr()
     ** calculates contribution to 
     **       H(op) * |initVec> ---> |initVec>
     ** from the precalculated nondiagonal matrix elements
     **  stored as float values in memory
     */

static void pnNondiagMemSD_Contr(int typeCalc, GR_BAS *grBas, GROUP *group,
                             double *initVec, double *finalVec)
{
  int           grNum, numZ, numN, pStep,nStep, 
                startAmp, amplStep,pnNumSD, initAmpl;
  GROUP         *grPtr;
  ULL           topBit;
  STORE_F       *matrPtr,*tempPtr;
  BARR_ELTYPE   *bitArray; 

  topBit = ULL_ONE << (8 * sizeof(ULL) - 1);
  matrPtr =  (typeCalc == ANG_CALC) ? grBas->permMemANG
                                  : grBas->permMem;
  for(grNum = 0; grNum < grBas->num_gr; grNum++) {
    grPtr = &group[grNum];
     if(group[grNum].RankNumSD_ZN_oscLim == 0) {
      continue; // no basis states in current group
    }
    numZ     = grPtr->numSD[PROTON]; // initialization
    numN     = grPtr->numSD[NEUTRON];
    bitArray = grPtr->bitArray;
    startAmp = grPtr->startAmp;
    amplStep = grPtr->numAmplStep; 

    for(pStep = 0; pStep < numZ;pStep++) {

      pnNumSD = pStep*numN + Rank;

      for(nStep = Rank;nStep < numN; nStep += amplStep,
	                           pnNumSD += amplStep) {
	if(BARR_TEST(bitArray,pnNumSD)) continue;

	initAmpl =  startAmp + pnNumSD 
                  - countBits(bitArray,pnNumSD);

	tempPtr = matrPtr;
	while(tempPtr->final != topBit) {
	  finalVec[(int)tempPtr->final] 
             += initVec[initAmpl]*(double)(tempPtr->value);
	  tempPtr++;
	}
      
	while(matrPtr->final != topBit) {
	  finalVec[initAmpl] 
             += initVec[(int)matrPtr->final]*(double)(matrPtr->value);
	  matrPtr++;
	}
      
	matrPtr++;
      } // end nStep

    }// end pStep

  }// grNum

} // end function  pnNondiagMemSD_Contr()


  /*
  ** The function
  **    countBits()
  ** counts and return number
  ** of bits set in a bitArray[]
  */

static ULL countBits(BARR_ELTYPE *bitArray, ULL size) 
{
  int  bit;
  ULL  count;

  count = 0;
  for(bit = 0; bit < size; bit++) {
    if(BARR_TEST(bitArray,bit)) count++;
  }
  return count;

} // End: function  countBits()
