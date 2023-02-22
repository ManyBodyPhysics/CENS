

/*********************** File: PAR-StoreSD_MatrElem.c *************/


#include "PAR-pnShellModel.h"

	    /**** local function declarations ****/

static void  pnCalcSD_DiagMatrElem(int typeCalc,GR_BAS *grBas,GROUP *group);
   /*
   ** calculates the diagonal SD matrix elements 
   **   <SD(Z)g,p;SD(N)g,n|OP|SD(Z)g,p;SD(N)g,n>
   ** with OP = OP(pp) + OP(nn) + OP(pn).
   ** The result is stored in memory pointed to by 
   ** grBAS->diagMem.
   */


static void occOrbits(int type,GR_BAS *grBas,GROUP *group);
    /* 
    ** calculates and stores orbit numbers of all occupied
    ** orbits  in all |SD(type)> in present group 
    ** -- complete |SD(N)> for all processes
    */

static void idTwoPartDiagSD_MatrElem(int type,int typeInt,GR_BAS *grBas,
                                          GROUP *group, double *diagPtr);
     /*
     ** calculates all diagonal matrix elements 
     **           <SD()|OP(pp)(type)| SD()> 
     ** where the two-particle matrix elements of OP(pp) are stored in
     **  grBas->op_pp.id_diag[] using the formula
     **        op(k,l) = ((num - k) * k)/2 + (l - 1)
     ** where num = 2 * m_orb - 3. 
     ** The result is stores in diagPtr[]. 
     */

static void  pnTwoPartDiagSD_MatrElem(int typeInt,GR_BAS *grBas,
              GROUP *group,double *ppDiagElem,double *nnDiagElem, 
				                  float *diagPtr);   /*
    ** calculates diagonal SD matrix elements 
    **    <SD(N),SD(Z)|OP(pn)| SD(Z),SD(N)> 
    ** The operator matrix elements are stored in pn_diag[] by the formula
    **        op(k,l) = num_orbN * k + l
    ** and the result is stores in diag_SD[]. 
    */

static void pppThreePartDiagSD_MatrElem(GR_BAS *grBas,GROUP *group);
    /*
    ** The function 
    **        pppThreePartDiagSD_MatrElem()
    ** calculates all diagonal matrix elements 
    **           <SD()|op_veff(three_part)()| SD()> 
    ** where the three-particle matrix elements of OP() are stored in
    ** op_veff->id_diag[] using the formula
    **  op(k,l,m) =  k*(num1 - k*(num2 -k))/6 + (l*(num3 - l))/2
    **             - num_orb + m     
    ** where num1 = 3*num_orb^2 -12*num_orb + 11
    **       num2 = 3*num_orb - 6
    **       num3 = 2*num_orb - 3
    ** The result is stored in memory pointed to by group->diagMemElem.
    */

static void ppnThreePartDiagSD_MatrElem(GR_BAS *grBas,GROUP *group);
    /*
    ** calculates diagonal SD matrix elements 
    **    <SD(N),SD(Z)|OP(ppn)| SD(Z),SD(N)> 
    */

static void pnnThreePartDiagSD_MatrElem(GR_BAS *grBas,GROUP *group);
    /*
    ** calculates diagonal SD matrix elements 
    **    <SD(N),SD(Z)|OP(ppn)| SD(Z),SD(N)> 
    */

static void nnnThreePartDiagSD_MatrElem(GR_BAS *grBas, GROUP *group);
    /*
    ** calculates all diagonal matrix elements 
    **           <SD()|op_veff(three_part)()| SD()> 
    ** where the three-particle matrix elements of OP() are stored in
    ** op_veff->id_diag[] using the formula
    **  op(k,l,m) =  k*(num1 - k*(num2 -k))/6 + (l*(num3 - l))/2
    **             - num_orb + m     
    ** where num1 = 3*num_orb^2 -12*num_orb + 11
    **       num2 = 3*num_orb - 6
    **       num3 = 2*num_orb - 3
    ** The result is stored in memory pointed to by sdStore->freePermMemPtr.
    */

static int pnSectFile(GR_BAS *grBas);
    /* 
    ** creates and open a new section file for
    ** storage of pnNondiag matrix elements with
    ** a corresponding section number
    ** If enough file space is available the 
    ** function return TRUE, otherwise FALSE
    */

static void  pnCalcSD_NondiagMatrElem(int typeInt,GR_BAS *grBas,
                                                   GROUP *group);
   /*
   ** calculates the nondiagonal SD matrix elements 
   **   <SD(Z)g',p';SD(N)g',n'|OP|SD(Z)g,p;SD(N)g,n>
   ** with OP = OP(pp) + OP(nn) + OP(pn).
   ** The result is stored in memory pointed to by 
   ** grBAS->permMem.
   */

 static int ppReadBlockMatrElemIntoMem(char *filename,int sectBlockNum,
				                  STORE_F *memBlockPtr);
     /*
     ** opens proton section file for nondiag
     ** matrix elements and reads numElem
     ** into memBlockPtr.
     ** The function returns numElem
     */

 static int nnReadBlockMatrElemIntoMem(char *filename,int sectBlockNum,
				                  STORE_F *memBlockPtr);
     /*
     ** opens neutron section file for nondiag
     ** matrix elements and reads numElem
     ** into memBlockPtr.
     ** The function returns numElem
     */

static int nnReadFirstBlockIntoMem(char *filename,STORE_F *memBlockPtr);
   /*
   ** reads first nondiag nnMatrElem data block 
   ** back into memBlockPtr[]
   */

static int singleGrNondiag(int typeCM,GR_BAS *grBas, GROUP *group,int nRankStep,
                         ULL sdZ,int sdZ_valN,ULL sdN,int sdN_valN,
			   STORE_F *pnPtr,ULL checkMemSize,
                         int numElemZ,STORE_F *memNondiagSD_Z,
			 int numElemN,STORE_F *RankMemNondiagSD_N);
   /* calculates the nondiag pn matrix elements
   ** within a single group from 
   **       V_pp, V_nn and V_pn
   ** the function returns the number of elements 
   ** calculated
   */

static int  multipleGrNondiag(int typeCM, GR_BAS *grBas,GROUP *initGroup,
	  GROUP *finalGroup,int nRankStep,ULL sdZ, ULL sdN,
	       STORE_F *pnPtr,ULL checkMemSize);
   /*
   ** calculates the nondiag pn matrix elements
   ** from initGroup to finalGroup due to V_pn
   ** The function returns the number of elements 
   ** calculated
   */

static int addTwoPart_ppConfig(ULL checkMemSize,GR_BAS *grBas, 
                                    GROUP *group,int sdN_valN,
			     STORE_F *pnPtr, STORE_F *pp_pStepPtr);
   /*
   ** converts the nondiag ppMatrix elements
   ** to a pn-scale and add them to the total
   ** list of pnNondiagMatrElemts
   */

static ULL countBits(BARR_ELTYPE *bitArray, ULL size);
   /*
   **  counts and return number
   ** of bits set in a bitArray[]
   */

static int addTwoPart_nnConfig(ULL checkMemSize,GR_BAS *grBas,
                                    GROUP *group,int sdZ_valN,
                             STORE_F *pnPtr, STORE_F *matrPtr);
   /*
   ** converts the nondiag nnMatrix elements
   ** to a pn-scale and add them to the total
   ** list of pnNondiagMatrElemts
   */

static int addTwoPart_pnConfig
             (int typeCM,ULL checkMemSize,int nRankStep,ULL sdZ, ULL sdN,
                  STORE_F *pnPtr,GR_BAS *grBas,GROUP *initGroup,
                             GROUP *finalGroup, PN_INT **table);
    /*
    ** calculates all nondiag pn matrix elements
    ** H(PN)*|init_pn> ---> value * |final_pn>
    ** within a single proton/neutron group
    */

static int pnSortNondiag(const STORE_F *one, const STORE_F *two);
    /*
    ** is a utility function for the library function qsort() in order to
    ** sort nondiagSD matrix elements of type STORE_F store[]
    ** after increasing store[].final.
    */

static int pnCompressNondiagSD_Elem(int num, STORE_F *store);
    /*
    ** takes num calculated nondiag SD elem which are sorted
    ** after increasing store[].final and add together all
    ** elements  with same store[].final
    ** Then all element store{].value < MATRIX_LIMIT
    ** are removed.  .
    ** Note: The function must have num > 1;
    ** The final number of elements is returned.
    */

static void printMatrElemData(int typeCalc,GR_BAS *sdBas,
                              TID wall_time,TID cpu_time);
    /*
    ** collect info about precalculation of <SD'|OP|SD>
    ** and send data to outputData() 
    */

static void outputData(int typeCalc,GR_BAS *grBas,int *data);
    /*
    ** send  info about precalculation of
    ** <SD'|OP|SD> to output file
    */

static int valOscN(int type,GR_BAS *grBas,ULL new_sd);
    /* 
    ** calculates sumN and return for
    ** |new_SD>.
    */

static void errorMemoryOverflow(GR_BAS *grBas,char *outFile);
    /*
    ** prints an error message to 
    ** program output file if no 
    ** more memory is available
    ** for data storage 
    */

	    /**** end local function declarations ****/


             /**** The function definitions  ****/ 

     /*
     ** The entrance function
     **     pnStoreSD_MatrElem()
     ** calculates and stores in memory all 
     ** diagonal proton-neutron matrix elements
     **     <SD(Z)g,p;SD(N)g,n|OP|SD(Z)g,p;SD(N)g,n>
     ** Next, all nondiagonal  proton-neutron
     ** matrix elements
     **     <SD(Z)g',p';SD(N)g',n'OP|SD(Z)g,p;SD(N)g,n> 
     ** stores in memory or on files. The nondiag 
     ** matrix elements are calculated in  two steps:
     **  1. All nondiagonal matrix elements for a single 
     *'     proton-neutron group
     **    <SD(Z)g,p';SD(N)g,n'OP|SD(Z)g,p;SD(N)g,n> 
     **  2. All nondiagonal matrix elements between two  
     *'     proton-neutron groups
     **    <SD(Z)g',p';SD(N)g',n'OP|SD(Z)g,p;SD(N)g,n> 
     */

void  pnStoreSD_MatrElem(int typeCalc,int calcCM,GR_BAS *grBas,
                                                  GROUP *group)
{
  char    *func = {" pnStoreSD_MatrElem(): "};
  TID     wallTime, cpuTime;


  if((typeCalc == VEFF_CALC) && (calcCM == NO))   {

    // diag memory for all groups:

    grBas->diagMemElem    // diag matrix element
      = MALLOC(grBas->RankTotNumSD_ZN_oscLim,float,func,
	               "grBas->diagMemElem[]");
         /*
	 ** permanent memory for nondiag
	 ** matrix element for all groups
	 ** . if(grBas->permMemSize > 0).
	 **    memory to store
	 **    <SD'|OP|SD>
	 */
  
    grBas->permMem = NULL;
    
    if(grBas->permMemSize > 0) {
      grBas->permMem
	= MALLOC(((size_t)grBas->permMemSize),STORE_F,func,"permMem[]");
    }
    grBas->curr_permMem     = grBas->permMem;
    grBas->curr_permMemSize = grBas->permMemSize;
  } // end memAlloc

  if((typeCalc == CM_CALC) && (calcCM == YES)) {
    grBas->curr_permMem     = grBas->permMem;
    grBas->curr_permMemSize = grBas->permMemSize;
  }
  else {
    if((typeCalc == ANG_CALC)&&(calcCM == NO))  { 
      grBas->diagMemElemANG   = grBas->diagMemElem;
      grBas->permMemANG       = grBas->permMem;
      grBas->curr_permMem     = grBas->permMemANG;
      grBas->curr_permMemSize = grBas->permMemSize;
    }
    else if((typeCalc == ANG_CALC)&&(calcCM == YES)) { 
      grBas->diagMemElemANG
	= MALLOC(grBas->RankTotNumSD_ZN_oscLim,float,func,
	                         "grBas->diagMemElemANG[]");
      grBas->curr_permMemSize = grBas->permMemSizeANG;
      grBas->curr_permMem     = grBas->permMemANG;
    }
  }
    // calculate and store diag matrix elements

  pnCalcSD_DiagMatrElem(typeCalc,grBas, group);

      /*
      ** calculate and store in memory or on files
      ** nondiagonal proton-neutron matrix elements 
      **  <SD(Z)g',p';SD(N)g,n'OP|SD(Z)g,p;SD(N)g,n> 
      */

  grBas->totMemNondiagElem  = 0;
  grBas->totFileNondiagElem = 0;

  // Proton/Neutron section files preparation

  grBas->pnSectFileNo = 0;

  pnSectFile(grBas); // create and open the first section file

  // time controll for calculating and store <SD'|OP|SD>

  wallClock(3,0);   // initialization
  cpuClock(3,0);
  wallClock(3,1);    // start clock
  cpuClock(3,1);

  pnCalcSD_NondiagMatrElem(typeCalc,grBas,group);

  // end time controll

  wallClock(3,2);  // stop clock
  cpuClock(3,2);
  wallTime = wallClock(3,3);  // read clock
  cpuTime  = cpuClock(3,3);

  //************  For Master only ****

  if(Rank == MASTER)   {
    printMatrElemData(typeCalc,grBas, wallTime, cpuTime); // end time control
  }

  //************  end MASTER only  *******

  fclose(grBas->pnSectFilePtr); // closed last section file

  if(grBas->totFileNondiagElem == 0) {
    remove(grBas->pnSectFilename);
  }

} //  End function pnStoreSD_MatrElem()

   /*
   ** The function
   **    pnCalcSD_DiagMatrElem()
   ** calculates the diagonal SD matrix elements 
   **   <SD(Z)g,p;SD(N)g,n|OP|SD(Z)g,p;SD(N)g,n>
   ** with OP = OP(pp) + OP(nn) + OP(pn).
   ** The result is stored in memory pointed to by 
   ** grBAS->diagMem.
   */

static void  pnCalcSD_DiagMatrElem(int typeCalc,GR_BAS *grBas,GROUP *group)
{
  char    *func = {" pnCalcSD_DiagMatrElem(): "};
  int     initGrNum;
  float   *pnDiagPtr;
  double  *ppDiagSD_Elem, *nnDiagSD_Elem;
  GROUP   *initGroup;

  // calculate and store diag matrix elements

  pnDiagPtr = (typeCalc == ANG_CALC) ? grBas->diagMemElemANG
                                    : grBas->diagMemElem;

  for(initGrNum = 0; initGrNum < grBas->num_gr; initGrNum++) {

    initGroup = &group[initGrNum];

    if(initGroup->RankNumSD_ZN_oscLim == 0) {
      continue; // no basis states in current process
    }
    // calculate and store occupied orbits in |SD(Z)> and |SD(N)> 

    initGroup->occ[PROTON] 
      = MALLOC((grBas->spBas[PROTON]->part+1)*initGroup->numSD[PROTON],
	      	                 int,func,"p_occ[PROTON]");
    occOrbits(PROTON,grBas,initGroup);

    initGroup->occ[NEUTRON] = NULL;
    if(initGroup->numSD[2] > 0) {
      initGroup->occ[NEUTRON] 
           = MALLOC((grBas->spBas[NEUTRON]->part+1)*initGroup->numSD[2],
                                       int,func,"n_occ[NEUTRON]");
      occOrbits(NEUTRON,grBas,initGroup);
    }

          /* 
	  ** Diagonal matrix elements:for current
	  ** group stored in memory pointed to 
	  ** group[initGrNum]->diagMemElem 
	  */

    ppDiagSD_Elem = 0;  // initialization
    if(grBas->spBas[PROTON]->part >= 2) {
      ppDiagSD_Elem = MALLOC(initGroup->numSD[0],double,
                              func,"ppDiagSD_Elem[PROTON]"); 
    }

    nnDiagSD_Elem = 0;  // initialization
    if(grBas->spBas[NEUTRON]->part >= 2) {
      nnDiagSD_Elem = MALLOC(initGroup->numSD[2],double,
                          func,"nnDiagSD_Elem[NEUTRON]"); 
    }

 //  contribution from three-particle interaction 

/**********************  Ikke fullført  ***********  
   if(grBas->typeCalc == VEFF3_LANC) {
      if(grBas->spBas[PROTON]->part >= 3) {
      pppThreePartDiagSD_MatrElem(grBas,initGroup);
    }
    if((grBas->spBas[PROTON]->part >= 2)&&(grBas->spBas[NEUTRON]->part >= 1)) {
      ppnThreePartDiagSD_MatrElem(grBas,initGroup);
    }
    if((grBas->spBas[PROTON]->part >= 1)&&(grBas->spBas[NEUTRON]->part >= 2)) {
      pnnThreePartDiagSD_MatrElem(grBas,initGroup);
    }
    if(grBas->spBas[NEUTRON]->part >= 3) {
      nnnThreePartDiagSD_MatrElem(grBas,initGroup);
    }

**********************************************/

    if(grBas->spBas[PROTON]->part >= 2) {
      idTwoPartDiagSD_MatrElem(PROTON,typeCalc,grBas,
                             initGroup,ppDiagSD_Elem);
    }
    if(grBas->spBas[NEUTRON]->part >= 2) {
      idTwoPartDiagSD_MatrElem(NEUTRON,typeCalc,grBas,
                              initGroup,nnDiagSD_Elem);
    }
    pnTwoPartDiagSD_MatrElem(typeCalc,grBas,initGroup,
                ppDiagSD_Elem,nnDiagSD_Elem,pnDiagPtr);

    if(grBas->spBas[NEUTRON]->part >= 2) {
      free(nnDiagSD_Elem);
    }
    if(grBas->spBas[PROTON]->part >= 2) {
      free(ppDiagSD_Elem);
    }

    pnDiagPtr += initGroup->RankNumSD_ZN_oscLim;


  } // end diag matrix elements for all groups

} // End function pnCalcSD_DiagMatrElem()

    /* 
    ** The function
    **       occOrbits_orbits()
    ** calculates and stores orbit numbers of all occupied
    ** orbits  in all |SD(type)> in present group 
    ** -- complete |SD(N)> for all processes
    */

static void occOrbits(int type,GR_BAS *grBas,GROUP *group)
{
  int   num_part,start,RankStep,*occ,loop,par,orb,sumN;
  ULL   *sd, pos,sd_value;
  MBAS  *mBasPtr;

  num_part = grBas->spBas[type]->part; // initialization
  occ      = group->occ[type];
  mBasPtr  = grBas->spBas[type]->mbas;

  if(type == PROTON) {
    start    = 0;
    RankStep = 1;
    sd       = group->SD[type];
  }
  else {
    start    = Rank;
    RankStep = NumProcs;
    sd       = group->SD[type] + Rank;
  }

  for(loop = start; loop < group->numSD[type]; loop += RankStep,
                                                 sd += RankStep) {
    sd_value = *sd;
    sumN     = 0;
    for(par = 0, orb = 0, pos = ULL_ONE; par < num_part;
                                 par++, orb++, pos <<= 1)  {
      for(;!(sd_value & pos); orb++, pos <<= 1);
      *(occ++) = orb; //particle found - save orbit
      sumN    += mBasPtr[orb].N; 
    }
    *(occ++) = sumN;
  } // end loop through all |SD> in sd[]

} // End: function occOrbits()

     /*
     ** The function 
     **         idTwoPartDiagSD_MatrElem()
     ** calculates all diagonal matrix elements 
     **           <SD()|OP(pp)(type)| SD()> 
     ** where the two-particle matrix elements of OP(type) are stored in
     **  grBas->op_pp.id_diag[] using the formula
     **        op(k,l) = ((num - k) * k)/2 + (l - 1)
     ** where num = 2 * m_orb - 3. 
     ** The result is stores in diagPtr[]. 
     */

static void idTwoPartDiagSD_MatrElem(int type,int typeCalc,GR_BAS *grBas,
                                           GROUP *group,double *diagPtr)
{
  int     num_part,num, num_part_1,loop,
          *occ, occSize, k, l, *k_ptr,*l_ptr, dim;
  double  *idDiagMatrElem, *table_ptr, value; 

  num_part   = grBas->spBas[type]->part; // initialization
  num_part_1 = num_part - 1;
  num        = (grBas->spBas[type]->numm_orb << 1) - 3;
  occ        = group->occ[type];
  occSize    = num_part + 1;

  if(type == PROTON) {
    idDiagMatrElem = grBas->op_pp.id_diag;

    if(typeCalc == CM_CALC) {
      idDiagMatrElem =  grBas->op_pp.id_CM_diag;
    }
    dim = group->numSD[0];
  }
  else {
    idDiagMatrElem = grBas->op_nn.id_diag;
    if(typeCalc == CM_CALC) {
      idDiagMatrElem =  grBas->op_nn.id_CM_diag;
    }
    dim = group->numSD[2];
  }

  for(loop = 0;loop < dim;loop++,occ += occSize,diagPtr++) {
    value = D_ZERO;
    k_ptr = occ;
    for(k = 0; k < num_part_1; k++) {
      table_ptr = idDiagMatrElem 
	         + (((num - (*k_ptr)) * (*k_ptr)) >> 1) -1;
      l_ptr = ++k_ptr;
      for(l = k+1; l < num_part; l++) {
	value += *(table_ptr + (*(l_ptr++)));
      } // end l-loop
    } // end k-loop
    *diagPtr = value;
  }

} // End: function idTwoPartDiagSD_MatrElem()

     /*
     ** The function 
     **      pnTwoPartDiagSD_MatrElem()
     ** calculates diagonal SD matrix elements 
     **    <SD(N),SD(Z)|OP(pn)| SD(Z),SD(N)> 
     ** The operator matrix elements are stored in pn_diag[] by the formula
     **        op(k,l) = num_orbN * k + l
     ** and the result is stores in group->diagMemElem. 
     */

static void pnTwoPartDiagSD_MatrElem(int typeCalc,GR_BAS *grBas,
            GROUP *group,double *ppDiagElem,double *nnDiagElem, 
                                                float *diagPtr)
{
  int     pStep, nStep, k, l, *k_ptr, *l_ptr,
          *Z_N, occStepZ_N, *N_N, occStepN_N;  
  double  *table, ppElem, nnElem,*pn_table, value; 
  
  table = grBas->op_pn.pn_diag;
  if(typeCalc == CM_CALC) {
    table = grBas->op_pn.pn_CM_diag;
  }
  occStepZ_N = grBas->spBas[0]->part + 1;
  occStepN_N = grBas->spBas[1]->part + 1;

  Z_N = group->occ[0] + grBas->spBas[0]->part;


 for(pStep = 0; pStep < group->numSD[0]; pStep++, Z_N += occStepZ_N) {

    ppElem = (grBas->spBas[PROTON]->part >= 2) ? ppDiagElem[pStep]: D_ZERO; 
	      
    N_N = group->occ[1] + grBas->spBas[1]->part;

    for(nStep = 0; nStep < group->numSD[2]; nStep++, N_N += occStepN_N) {

      if(((*Z_N) + (*N_N))> grBas->maxOscN) continue;

      nnElem = (grBas->spBas[NEUTRON]->part >= 2) ? nnDiagElem[nStep]: D_ZERO; 

      k = grBas->spBas[0]->part;
      k_ptr = group->occ[0] + occStepZ_N*pStep;
      value = D_ZERO;
      do {  // proton loop
	pn_table = table + (*k_ptr)*grBas->spBas[1]->numm_orb;
	k_ptr++;
	l     = grBas->spBas[1]->part;
	l_ptr =  group->occ[1] + occStepN_N*nStep;
	do {
	  value += *(pn_table + (*l_ptr));
          l_ptr++;

	} while(--l); // neutron particle loop

      } while(--k); // proton particle loop

      *(diagPtr++) = (float)(value + ppElem + nnElem);

    } // end run through all |SD(N)>
	
  } // end run through all |SD(Z),SD(N)>
  
} // End: function pnTwoPartDiagSD_MatrElem()

     /*
     ** The function 
     **        pppThreePartDiagSD_MatrElem()
     ** calculates all diagonal matrix elements 
     **           <SD()|op_veff(three_part)()| SD()> 
     ** where the three-particle matrix elements of OP() are stored in
     ** op_veff->id_diag[] using the formula
     **  op(k,l,m) =  k*(num1 - k*(num2 -k))/6 + (l*(num3 - l))/2
     **             - num_orb + m     
     ** where num1 = 3*num_orb^2 -12*num_orb + 11
     **       num2 = 3*num_orb - 6
     **       num3 = 2*num_orb - 3
     ** The result is stored in memory pointed to by group->diagMemElem.
     */

static void pppThreePartDiagSD_MatrElem(GR_BAS *grBas,GROUP *group)
{
  int     num_part, num_part_2, loop, *occ, 
          num1, num2,num3, num_orb, k, l, m,  *k_ptr,*l_ptr, *m_ptr;
  double  *table_k_ptr, *table_kl_ptr, value; 
  float   *diagPtr;

  num_part   = grBas->spBas[PROTON]->part; // initialization
  num_part_2 = num_part - 2;

  num_orb = grBas->spBas[PROTON]->numm_orb;
  num1    = 3*num_orb*num_orb -12*num_orb + 11;
  num2    = 3*num_orb - 6;
  num3    = 2*num_orb - 3;
  occ     = group->occ[PROTON];

  diagPtr = grBas->diagMemElem;
  for(loop = 0;loop < group->numSD[PROTON];loop++,occ += num_part) {
    value = D_ZERO; // <SD|OP()|SD> = 0
    k     = num_part_2;
    k_ptr = occ;
    do {      // k-particle loop
      table_k_ptr = grBas->op_pp.id_diag 
                     + ((*k_ptr)*(num1 - (*k_ptr)*(num2 - (*k_ptr))))/6;
      l     = k;
      l_ptr = (++k_ptr);
      do {    // l-particle loop
	table_kl_ptr = table_k_ptr + ((*l_ptr)*(num3 - (*l_ptr)))/2;
	m = l;
	m_ptr = (++l_ptr);
	do {  //  m-particle loop
	  value += *(table_kl_ptr + (*(m_ptr++) - num_orb));
	} while(--m); // end m particle loop 
      } while(--l);   // end l particle loop 
    } while(--k);     // end k particle loop

    k = group->numSD[2];
    do {
      *diagPtr += (float)value;
      diagPtr++;
    } while(--k);
  }

} // End: function  pppThreePartDiagSD_MatrElem()

     /*
     ** The function 
     **      ppnThreePartDiagSD_MatrElem()
     ** calculates diagonal SD matrix elements 
     **    <SD(N),SD(Z)|OP(ppn)| SD(Z),SD(N)> 
     */

static void ppnThreePartDiagSD_MatrElem(GR_BAS *grBas,GROUP *group)
{
   
  //      Not yet implemented

  printf("\n\nRank%d: Error from ppnThreePartDiagSD_MatrElem(): ",Rank);
  printf("\n Not yet implemented\n\n");
  MPI_Abort(MPI_COMM_WORLD,Rank);

} // End: function  ppnThreePartDiagSD_MatrElem()

     /*
     ** The function 
     **      pnnThreePartDiagSD_MatrElem()
     ** calculates diagonal SD matrix elements 
     **    <SD(N),SD(Z)|OP(ppn)| SD(Z),SD(N)> 
     */

static void pnnThreePartDiagSD_MatrElem(GR_BAS *grBas,GROUP *group)
{
   
  //      Not implemented

  printf("\n\nRank%d: Error from pnnThreePartDiagSD_MatrElem(): ",Rank);
  printf("\n Not yet implemented\n\n");
  MPI_Abort(MPI_COMM_WORLD,Rank);

} // End: function  pnnThreePartDiagSD_MatrElem()

     /*
     ** The function 
     **        nnnThreePartDiagSD_MatrElem()
     ** calculates all diagonal matrix elements 
     **           <SD()|op_veff(three_part)()| SD()> 
     ** where the three-particle matrix elements of OP() are stored in
     ** op_veff->id_diag[] using the formula
     **  op(k,l,m) =  k*(num1 - k*(num2 -k))/6 + (l*(num3 - l))/2
     **             - num_orb + m     
     ** where num1 = 3*num_orb^2 -12*num_orb + 11
     **       num2 = 3*num_orb - 6
     **       num3 = 2*num_orb - 3
     ** The result is stored in memory pointed to by sdStore->freePermMemPtr.
     */

static void nnnThreePartDiagSD_MatrElem(GR_BAS *grBas, GROUP *group)
{
  int     num_part,loop, *occ,  num_part_2,
          num1, num2,num3, num_orb, k, l, m,  *k_ptr,*l_ptr, *m_ptr;
  double  *table_k_ptr, *table_kl_ptr, value; 
  float   *diag, *diagPtr;

  num_part   = grBas->spBas[NEUTRON]->part; // initialization
  num_part_2 = num_part - 2;

  num_orb = grBas->spBas[NEUTRON]->numm_orb;
  num1    = 3*num_orb*num_orb -12*num_orb + 11;
  num2    = 3*num_orb - 6;
  num3    = 2*num_orb - 3;
  occ     = group->occ[1];

  diag     = grBas->diagMemElem; 
  for(loop = 0; loop < group->numSD[2]; loop++,occ += num_part, diag++) {
    value = D_ZERO;
    k     = num_part_2;
    k_ptr = occ;
    do {      // k-particle loop
      table_k_ptr =  grBas->op_nn.id_diag 
                     + ((*k_ptr)*(num1 - (*k_ptr)*(num2 - (*k_ptr))))/6;
      l     = k;
      l_ptr = (++k_ptr);
      do {    // l-particle loop
	table_kl_ptr = table_k_ptr + ((*l_ptr)*(num3 - (*l_ptr)))/2;
	m = l;
	m_ptr = (++l_ptr);
	do {  //  m-particle loop
	  value += *(table_kl_ptr + (*(m_ptr++) - num_orb));
	} while(--m); // end m particle loop 
      } while(--l);   // end l particle loop 
    } while(--k);     // end k particle loop
  

   diagPtr = diag;
    k =  group->numSD[0];
    do {
      *diagPtr += (float)value;
      diagPtr  += group->numSD[0];
    } while(--k);
 
  } // end NEUTRON loop

} // End: function  nnnThreePartDiagSD_MatrElem()

     /* 
     ** The function
     **      pnSectFile()
     ** creates and open a new section file for
     ** storage of pnNondiag matrix elements with
     ** a corresponding section number
     ** If enough file space is available the 
     ** function return TRUE, otherwise FALSE
     */

static int pnSectFile(GR_BAS *grBas)
{
  if(grBas->pnSectFileNo > 0) {

/********  ikke i bruk  ennå  **********
    grBas->pnNondiagTotFileSize 
      -=(grBas->pnSectFileSize - grBas->curr_pnSectFileSize);
    if(grBas->pnNondiagTotFileSize < grBas->blockMemSize) {
      return FALSE; // no more space for pn section files
    }
**************************/

    fclose(grBas->pnSectFilePtr);
  }

  grBas->pnSectFileNo++;      // new section file number
  sprintf(grBas->pnSectFilename,"%s%s%s%d",grBas->title,"pn",
	                                     NONDIAG,grBas->pnSectFileNo);

  grBas->pnSectFileSize = MIN(grBas->pnNondiagTotFileSize,
                             (MAX_FILE_SIZE*M_BYTES)/sizeof(STORE_F)); 
  grBas->curr_pnSectFileSize = grBas->pnSectFileSize; 

  // create and open new pn section file to store nondiag pnMatrElem

  if((grBas->pnSectFilePtr = fopen(grBas->pnSectFilename,"ab")) == NULL) {
    printf("\n\nRank%d: Error in function pnSectFile():",Rank);
    printf("\nNot allowed to open %s\n", grBas->pnSectFilename);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }
  return  TRUE;

} // End: function  pnSectFile()

   /*
   ** The function
   **    pnCalcSD_NondiagMatrElem()
   ** calculates the nondiagonal SD matrix elements 
   **   <SD(Z)g',p';SD(N)g',n'|OP|SD(Z)g,p;SD(N)g,n>
   ** with OP = OP(pp) + OP(nn) + OP(pn).
   ** The result is stored in memory pointed to by 
   ** grBAS->permMem.
   */

static void  pnCalcSD_NondiagMatrElem(int typeCalc,GR_BAS *grBas,
                                                    GROUP *group)
{
  char      *func = {"pnCalcSD_NondiagMatrElem()"};
  int       initGrNum, finalGrNum, maxDelM_Z,
            ppNumMemElem, nnNumMemElem,checkMemSize,
            ppMemBlockSize,nnMemBlockSize,
            pStep,nStep,nRankStep,pnNumElem, numElem, typeCM = 0, 
            *Z_N, occStepZ_N, *N_N, occStepN_N;
  ULL       *sdZ, *sdN,topBit;
  GROUP     *initGroup, *finalGroup;
  STORE_F   *pnPtr,*ppMemBlockPtr,*pp_pStepPtr = NULL,
                 *nnMemBlockPtr,*nn_nStepPtr = NULL;

  topBit  = ULL_ONE << (8 * sizeof(ULL) - 1);

/****************  output pp-test **************

{
  char     filename[ONE_LINE];
  FILE     *filePtr;

  sprintf(filename,"%s%s",grBas->title,RESULT_OUTPUT);

  if((filePtr = fopen(filename,"a"))== NULL) {
    printf("\n\nRank%d: Error in ????????:",Rank);
    printf("\nWrong file = %s for the output data file\n",filename);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }
  fprintf(filePtr,"\n\n\n\n  START pnCalcSD_NondiagMatrElem()\n\n",Rank);

  fclose(filePtr);
 }

****************   end output test **************/

  if((typeCalc == VEFF_CALC)||(typeCalc == ANG_CALC)) {
    typeCM = NO;
  }
  else if(typeCalc == CM_CALC)  {
    if(grBas->calc_CM != 1) {
      printf("\n\nRank%d: Error in pnCalcSD_NondiagMatrElem():",Rank);
      printf("\nNo Center-Of-Mass matrix elements available\n\n");
      MPI_Abort(MPI_COMM_WORLD,Rank);
    }
    typeCM = YES;
  }

  // PROTON: Temporary memory to store pp-NondiagMatrElem

  ppMemBlockSize =(grBas->maxGroupNumSD_Z*(grBas->maxGroupNumSD_Z + 1))/2;
  ppMemBlockPtr  = NULL_PTR;
  if(ppMemBlockSize > 1) {
    ppMemBlockPtr = MALLOC(ppMemBlockSize,STORE_F,func,"ppMemBlock[]");
  } // end PROTON case

  // NEUTRON: Temporary memory to store nn-NondiagMatrElem

  nnMemBlockSize =(grBas->maxGroupNumSD_N*(grBas->maxGroupNumSD_N +1))/2;
  if(NumProcs > 1) nnMemBlockSize /= (NumProcs - 1);
  nnMemBlockPtr = NULL_PTR;

  if(nnMemBlockSize > 1) {
    nnMemBlockPtr
      = MALLOC(nnMemBlockSize,STORE_F,func,"nnMemBlock[]");
  } // end NEUTRON case
  
// maximum change in M_Z between initial and final groups

  maxDelM_Z = 2*grBas->spBas[0]->mbas->m;


  for(initGrNum = 0; initGrNum < grBas->num_gr; initGrNum++) {
    initGroup = &group[initGrNum];
    if((initGroup->RankNumSD_ZN_oscLim) == 0) continue;// no basis states

    //   NONDIAG PROTON-PROTON CASE

    ppNumMemElem = 0;   // initialization

/****************  output pp-test **************

{
  char     filename[ONE_LINE];
  FILE     *filePtr;

  sprintf(filename,"%s%s",grBas->title,RESULT_OUTPUT);

  if((filePtr = fopen(filename,"a"))== NULL) {
    printf("\n\nRank%d: Error in ????????:",Rank);
    printf("\nWrong file = %s for the output data file\n",filename);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }
  fprintf(filePtr,"\n\n\n\n  Start PROTON: func  idTwoPartNondiagSD_Calc()\n\n",Rank);
  fprintf(filePtr,"\n initGrNum = %d ppMemBlockSize = %d ppMemBlockPtr = %x",
	  initGrNum,ppMemBlockSize,ppMemBlockPtr);

  fclose(filePtr);
 }

****************   end output test **************/


    if(  (grBas->spBas[PROTON]->part >= 2)
       &&(initGroup->numSD[PROTON] > 1))    {

      ppNumMemElem = idTwoPartNondiagSD_Calc(PROTON,typeCM,grBas,
                                        initGroup, ppMemBlockPtr);

/****************  output pp-test **************

{
  char     filename[ONE_LINE];
  FILE     *filePtr;

  sprintf(filename,"%s%s",grBas->title,RESULT_OUTPUT);

  if((filePtr = fopen(filename,"a"))== NULL) {
    printf("\n\nRank%d: Error in ????????:",Rank);
    printf("\nWrong file = %s for the output data file\n",filename);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }
  fprintf(filePtr,"\n  END PROTON: func  idTwoPartNondiagSD_Calc()\n\n",Rank);
  fprintf(filePtr,"\n  ppNumMemElem = %d ",ppNumMemElem);
  fclose(filePtr);
 }

****************   end output test **************/



    }  // end PROTON case

    // NONDIAG NEUTRON_NEUTRON CASE

    nnNumMemElem   = 0;   // initialization

/****************  output pp-test **************

{
  char     filename[ONE_LINE];
  FILE     *filePtr;

  sprintf(filename,"%s%s",grBas->title,RESULT_OUTPUT);

  if((filePtr = fopen(filename,"a"))== NULL) {
    printf("\n\nRank%d: Error in ????????:",Rank);
    printf("\nWrong file = %s for the output data file\n",filename);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }
  fprintf(filePtr,"\n\n  Start NEUTRON: func  idTwoPartNondiagSD_Calc()\n\n",Rank);
  fprintf(filePtr,"\n initGrNum = %d  nnMemBlockSize = %d nnMemBlockPtr = %x",
	  initGrNum,nnMemBlockSize,nnMemBlockPtr);

  fclose(filePtr);
 }

****************   end output test **************/




    if(  (grBas->spBas[NEUTRON]->part >= 2)
       &&(initGroup->numSD[NEUTRON]) > 1)  {

      nnNumMemElem= idTwoPartNondiagSD_Calc(NEUTRON,typeCM,grBas,
                                         initGroup,nnMemBlockPtr);

/****************  output pp-test **************

{
  char     filename[ONE_LINE];
  FILE     *filePtr;

  sprintf(filename,"%s%s",grBas->title,RESULT_OUTPUT);

  if((filePtr = fopen(filename,"a"))== NULL) {
    printf("\n\nRank%d: Error in ????????:",Rank);
    printf("\nWrong file = %s for the output data file\n",filename);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }
  fprintf(filePtr,"\n  END NEUTRON: func  idTwoPartNondiagSD_Calc()\n\n",Rank);
  fprintf(filePtr,"\n  nnNumMemElem = %d   nnData.numTotalElem = %d",nnNumMemElem);
  fclose(filePtr);
 }

****************   end output test **************/

    } // end NEUTRON case 

    initGroup->storedSD_Elem = 0;  // initialization
    occStepZ_N = grBas->spBas[0]->part + 1; 
    occStepN_N = grBas->spBas[1]->part + 1;
    Z_N        = initGroup->occ[0] + grBas->spBas[0]->part;
    sdZ = initGroup->SD[PROTON];

    //   Find ppNondiag matr elem

    pp_pStepPtr = NULL_PTR;
    if(ppNumMemElem > 0) {
      pp_pStepPtr = ppMemBlockPtr;
    }
    for(pStep = 0; pStep < initGroup->numSD[PROTON]; pStep++, sdZ++,
                                                 Z_N += occStepZ_N ) {
      if((pStep > 0)&&(ppNumMemElem > 0)) { // search for next topBit 
	while(pp_pStepPtr->final != topBit) pp_pStepPtr++;
	pp_pStepPtr++;
      } // ppCurrBlockPtr points to next nondiag ppElem

      initGroup->currSD_Z = pStep; 
      N_N = initGroup->occ[1] + grBas->spBas[1]->part;
      sdN =  initGroup->SD[NEUTRON] + Rank;

      // find nnNondiag matr elem

      nn_nStepPtr = NULL_PTR;
      if(nnNumMemElem > 0) {
	nn_nStepPtr = nnMemBlockPtr;
      }

      for(nStep = Rank,nRankStep = 0;nStep < initGroup->numSD[NEUTRON]; 
	                nStep += NumProcs,nRankStep++, sdN += NumProcs,
	                                             N_N += occStepN_N) {
       if((nStep > Rank)&&(nnNumMemElem > 0)) { // search for next topBit 
	  while(nn_nStepPtr->final != topBit) nn_nStepPtr++;
	  nn_nStepPtr++;
	} // ppCurrBlockPtr points to next nondiag ppElem
    
	// limits on oscillator excitation

	if(((*Z_N)+(*N_N)) > grBas->maxOscN) continue;

	initGroup->currSD_N = nStep;
	pnPtr               = grBas->curr_permMem;
	checkMemSize        = grBas->curr_permMemSize;
	pnNumElem           = 0;

	numElem = singleGrNondiag(typeCM,grBas,initGroup,nRankStep,
                                 *sdZ,*Z_N,*sdN,*N_N,
				  pnPtr,checkMemSize,
                                  ppNumMemElem,pp_pStepPtr,
                                  nnNumMemElem,nn_nStepPtr);

	pnNumElem    += numElem;
	pnPtr        += numElem;
	checkMemSize -= numElem;


	for(finalGrNum = initGrNum + 1; finalGrNum < grBas->num_gr;
                                                       finalGrNum++) {
	  finalGroup = &group[finalGrNum];

	  // check for change in m-values

	  if((initGroup->m[0]-finalGroup->m[0]) > maxDelM_Z)
	    continue;
	  numElem = multipleGrNondiag(typeCM,grBas,initGroup,finalGroup,
                       nRankStep,*sdZ,*sdN,pnPtr,checkMemSize);


	  pnNumElem    += numElem;
	  pnPtr        += numElem;
	  checkMemSize -= numElem;

	}  // end finalGroup calc

         /* 
	 ** sort and compress nondiag matrix elements for 
	 ** each |init_SD> after increasing |final_SD>
	 */

	if(pnNumElem > 1) {
	  qsort(grBas->curr_permMem, (size_t)pnNumElem, (size_t)sizeof(STORE_F),
	      (int(*)(const void *, const void *))pnSortNondiag);
	  pnNumElem = pnCompressNondiagSD_Elem(pnNumElem,grBas->curr_permMem);
	}

	grBas->curr_permMem       += pnNumElem;
	grBas->curr_permMem->value = D_ZERO;
	grBas->curr_permMem->final = topBit;
	grBas->curr_permMem++;
	pnNumElem++;

	grBas->totMemNondiagElem += pnNumElem;
	grBas->curr_permMemSize  -= pnNumElem;
	initGroup->storedSD_Elem += pnNumElem;

      } // end nStep-loop

    } // end pStep-loop

    free(initGroup->occ[1]); 
    initGroup->occ[1] = NULL_PTR;
    free(initGroup->occ[0]);
    initGroup->occ[0] = NULL_PTR;


/****************  output initGrNum-test **************

if(Rank == MASTER) {

   {
    char     filename[ONE_LINE];
    FILE     *filePtr;

    sprintf(filename,"%s%s",grBas->title,RESULT_OUTPUT);

    if((filePtr = fopen(filename,"a"))== NULL) {
      printf("\n\nRank%d: Error in ????????:",Rank);
      printf("\nWrong file = %s for the output data file\n",filename);
      MPI_Abort(MPI_COMM_WORLD,Rank);
    }
    fprintf(filePtr,"\n\n\n Number of NondiagMatrix for initGrNum = %d",initGrNum);
    fprintf(filePtr,"\n ppMemBlockSize = %d  ppNumMemElem = %d nnMemBlockSize = %d ppNumMemElem = %d",
	    ppMemBlockSize,ppNumMemElem,nnMemBlockSize,nnNumMemElem);
    fprintf(filePtr,"\n pnNumElem = %d  checkMemsize = %d",pnNumElem,checkMemSize);
    fclose(filePtr);
   }

 }  // end MASTER


****************   end output test **************/


  } // end initGroup calc
 
  if(nnMemBlockSize > 1) {
    free(nnMemBlockPtr);
  }
  if(ppMemBlockSize > 1) {
    free(ppMemBlockPtr);
  }
  
  // reserve remaining permMem to store <SD'|ANG|SD>
 
  if(typeCalc == CM_CALC) {
    grBas->permMemSizeANG = grBas->curr_permMemSize;
    grBas->permMemANG     = grBas->curr_permMem;
  }

} // End: function pnCalcSD_NondiagMatrElem()

     /*
     ** The function
     **     ppReadBlockMatrElemIntoMem()
     ** opens proton section file for nondiag
     ** matrix elements and reads numElem
     ** into memBlockPtr.
     ** The function returns numElem
     */

static int ppReadBlockMatrElemIntoMem(char *filename,int sectBlockNum,
                                                  STORE_F *memBlockPtr)
{
  int      numElem,lastBlock;
  ULL      sectBit;
  STORE_F  *memPtr;
  static   FILE *filePtr;

  sectBit = ULL_ONE << (8 * sizeof(ULL) - 2);
  lastBlock = NO;
  if(sectBlockNum == 0) {
    if((filePtr = fopen(filename,"rb")) == NULL) {
      printf("\n\nRank%d: Error in function ppReadBlockMatrElemIntoMem();",Rank);
      printf("\nWrong file = %s to read nondiag ppMatrElem\n\n",filename);
      MPI_Abort(MPI_COMM_WORLD,Rank);
    }
    rewind(filePtr);

    if(fread((void *)&numElem,(size_t) sizeof(int), 1,filePtr) != (size_t)1) {
      printf("\n\nRank%d: Error in function ppReadBlockMatrElemIntoMem();",Rank);
      printf("\nin reading first ppMatrElem in file %s\n\n", filename);
      MPI_Abort(MPI_COMM_WORLD,Rank);
    }
    fseek(filePtr,(long)(numElem*sizeof(STORE_F)),SEEK_SET);

  } // end sectBlockNum = 0
  
  if(fread((void *)&numElem,(size_t) sizeof(int), 1,filePtr) != (size_t)1) {
    printf("\n\nRank%d: Error in function ppReadBlockMatrElemIntoMem();",Rank);
    printf("\nin readin first ppMatrElemin file %s\n\n", filename);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }
  if(numElem < 0)  {
    lastBlock = YES;
    numElem = -numElem;
  }
  memPtr = memBlockPtr;
  if(fread((void *)memPtr,(size_t) sizeof(STORE_F),numElem ,filePtr)
                                            != (size_t)numElem) {
    printf("\n\nRank%d: Error in function ppReadBlockMatrElemIntoMem();",Rank);
    printf("\nin readin first ppMatrElem in file %s\n\n", filename);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }
  memBlockPtr       += numElem; 
  memBlockPtr->value = D_ZERO; 
  memBlockPtr->final = sectBit; 

  if(lastBlock == YES) fclose(filePtr);

  return numElem;

} // End: function  ppReadBlockMatrElemIntoMem

     /*
     ** The function
     **     nnReadBlockMatrElemIntoMem()
     ** opens neutron section file for nondiag
     ** matrix elements and reads numElem
     ** into memBlockPtr.
     ** The function returns numElem
     */

static int  nnReadBlockMatrElemIntoMem(char *filename,int sectBlockNum,
                                                  STORE_F *memBlockPtr)
{
  int      numElem, lastBlock;
  ULL      sectBit, topBit;
  STORE_F  *memPtr;
  static   FILE *filePtr;

  topBit = ULL_ONE << (8 * sizeof(ULL) - 1);
  sectBit = ULL_ONE << (8 * sizeof(ULL) - 2);
  lastBlock = NO;
  if(sectBlockNum == 0) {
    if((filePtr = fopen(filename,"rb")) == NULL) {
      printf("\n\nRank%d: Error in function nnReadBlockMatrElemIntoMem();",Rank);
      printf(" \nWrong file = %s to read nondiag ppMatrElem\n\n",filename);
      MPI_Abort(MPI_COMM_WORLD,Rank);
    }
    rewind(filePtr);

    if(fread((void *)&numElem,(size_t) sizeof(int), 1,filePtr) != (size_t)1) {
      printf("\n\nRank%d: Error in function nnReadBlockMatrElemIntoMem();",Rank);
      printf("\nin reading first nnMatrElem in file %s\n\n", filename);
      MPI_Abort(MPI_COMM_WORLD,Rank);
    }
    fseek(filePtr,(long)(sizeof(int) + numElem*sizeof(STORE_F)),SEEK_SET);

  } // end sectBlockNum = 0
  
  if(fread((void *)&numElem,(size_t) sizeof(int), 1,filePtr) != (size_t)1) {
    printf("\n\nRank%d: Error in function nnReadBlockMatrElemIntoMem();",Rank);
    printf("\nin readin first nnMatrElem in file %s\n\n", filename);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }
  if(numElem < 0)  {
    lastBlock = YES;
    numElem = -numElem;
  }
  memPtr = memBlockPtr;
  if(fread((void *)memPtr,(size_t) sizeof(STORE_F),numElem ,filePtr)
                                           != (size_t)numElem) {
    printf("\n\nRank%d: Error in function nnReadBlockMatrElemIntoMem();",Rank);
    printf("\nin readin first nnMatrElem in file %s\n\n", filename);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }
  if(lastBlock == NO) {
    memBlockPtr        += numElem; 
    memBlockPtr->value  = D_ZERO;
    memBlockPtr->final  = sectBit;
  }

  if(lastBlock == YES) fclose(filePtr);
  
  return numElem;

} // End: function  nnReadBlockMatrElemIntoMem

   /*
   ** The function
   **  nnReadFirstFileBlockIntoMem() 
   ** reads first nondiag nnMatrElem data block 
   ** back into memBlockPtr[]
   */

static int nnReadFirstBlockIntoMem(char *filename,STORE_F *memBlockPtr)
{
  int   numElem;
  ULL   sectBit;
  FILE  *filePtr;

  sectBit = ULL_ONE << (8 * sizeof(ULL) - 2);

  if((filePtr = fopen(filename,"rb")) == NULL) {
    printf("\n\nRank%d: Error in function nnReadFirstBlockIntoMem():",Rank);
    printf("\nWrong file = %s to read nondiag first nnNobdiagMatrElem\n\n\n",
                                                                    filename);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }
  rewind(filePtr);

  if(fread((void *)&numElem,(size_t) sizeof(int), 1,filePtr) != (size_t)1) {
    printf("\n\nRank%d: Error in function nnReadFirstBlockIntoMem();",Rank);
      printf(" in reading numElem in first nnDataBlock in %s\n\n", filename);
      MPI_Abort(MPI_COMM_WORLD,Rank);
  }

  if(fread((void *)memBlockPtr,(size_t) sizeof(STORE_F), numElem,filePtr)
                                                != (size_t)numElem) {
    printf("\n\nRank%d: Error in function  nnReadFirstBlockIntoMem();",Rank);
    printf(" in reading first nnMatrElem into Mem %s\n\n", filename);
    MPI_Abort(MPI_COMM_WORLD,Rank);
  }

  fclose(filePtr);

  // Add an extra sectBit  to matrElem data in memBlockPtr[]

  memBlockPtr       += numElem;
  memBlockPtr->value = D_ZERO;
  memBlockPtr->final = sectBit;

  return numElem;

} // End: function nnReadFirstBlockIntoMem()

   /*
   ** The function 
   **  singleGrNondiag
   ** calculates the nondiag pn matrix elements
   ** within a single group from 
   **       V_pp, V_nn and V_pn
   ** the function returns the number of elements 
   ** calculated
   */

static int singleGrNondiag(int typeCM, GR_BAS *grBas, GROUP *group,
                         int nRankStep,ULL sdZ,int sdZ_valN,ULL sdN,
                         int sdN_valN,STORE_F *pnPtr,ULL checkMemSize,
                         int numElemZ,STORE_F *pp_pStepPtr,
		         int numElemN,STORE_F *nn_nStepPtr)
{
  int   numElem, pnNumElem, pStep, nStep,
        ppStepLimit, nnStepLimit;

  pStep        = group->currSD_Z;
  nStep        = group->currSD_N;
  ppStepLimit  = group->numSD[PROTON] - 1; 
  nnStepLimit  = group->numSD[NEUTRON] - 1;
  pnNumElem    = 0;   // initialization

  if(  (grBas->spBas[PROTON]->part >= 2)
     &&(numElemZ > 0)
     &&(pStep < ppStepLimit))  {

    numElem = 0;
    if(numElemZ > 0) {
      numElem = addTwoPart_ppConfig(checkMemSize,grBas,
                                        group,sdN_valN,
                                     pnPtr,pp_pStepPtr);
    }
    pnNumElem    += numElem;
    pnPtr        += numElem;
    checkMemSize -= numElem;
  }   // end pp_case

  numElem = 0;
  if(  (grBas->spBas[NEUTRON]->part >= 2)
	   &&(nStep < nnStepLimit)
	   &&(numElemN > 0)) {

    if(numElemN > 0) {
      numElem = addTwoPart_nnConfig(checkMemSize,grBas,
                                        group,sdZ_valN,
				     pnPtr,nn_nStepPtr);
    }
    pnNumElem    += numElem;
    pnPtr        += numElem;
    checkMemSize -= numElem;
  }  // end nn_case 

  numElem = addTwoPart_pnConfig(typeCM,checkMemSize,nRankStep,
			 sdZ,sdN,pnPtr,grBas,group,group,
                     grBas->op_pn.kl_list[0]);

  pnNumElem    += numElem;
  pnPtr        += numElem;
  checkMemSize -= numElem;
  

  return pnNumElem;

} // End: function singleGrNondiag()

   /*
   ** The function 
   **  multipleGrNondiag
   ** calculates the nondiag pn matrix elements
   ** from initGroup to finalGroup due to V_pn
   ** The function returns the number of elements 
   ** calculated
   */

static int  multipleGrNondiag(int typeCM, GR_BAS *grBas,GROUP *initGroup,
	  GROUP *finalGroup,int nRankStep,ULL sdZ, ULL sdN,
               STORE_F *pnPtr,ULL checkMemSize)
{
  int    numElem,groupStep_m_p;  

  // changes in M_Z and PAR_Z from init to final group

  if(grBas->parZ) {
    groupStep_m_p =  (initGroup->m[0] - finalGroup->m[0]) >> 1;
  }
  else {
    groupStep_m_p =  (initGroup->m[0] - finalGroup->m[0]) 
	      + (((initGroup->par[0] * finalGroup->par[0]) == +1) 
		 ? 0 : 1);
  }

  numElem = addTwoPart_pnConfig(typeCM, checkMemSize, nRankStep,
                 sdZ,sdN,pnPtr,grBas,initGroup,finalGroup,
                      grBas->op_pn.kl_list[groupStep_m_p]);
  return  numElem;  

} // End: function multipleGrNondiag()

   /*
   ** The function
   **     addTwoPart_ppConfig()
   ** converts all nondiag ppMatrix elements from
   ** current init |pStep> to a pn-scale and add 
   ** them to the total list of 
   ** pnNondiagMatrElemts.
   ** matrPtr points to start position of current
   ** nondiag ppMatrix elements - topBit terminates 
   */

static int addTwoPart_ppConfig(ULL checkMemSize,GR_BAS *grBas,
                                    GROUP *group,int sdN_valN, 
                             STORE_F *pnPtr, STORE_F *pp_pStepPtr)
{
  int     nStep,numElem, numOneBits;
  ULL     topBit, newZ_sd, newZ_N, newPos;

  topBit = ULL_ONE << (8 * sizeof(ULL) - 1);
  nStep  = group->currSD_N;

  numElem = 0;
  
  for( ;((*pp_pStepPtr).final != topBit);pp_pStepPtr++) {

    // check oscillator exitation for |newZ_sd,newN_sd>

    newZ_sd = *(group->SD[0] + (*pp_pStepPtr).final);
    newZ_N  = valOscN(PROTON,grBas,newZ_sd);

    if((newZ_N + sdN_valN) > (ULL)grBas->maxOscN) continue;

    (*pnPtr).value = (*pp_pStepPtr).value;
    newPos = (*pp_pStepPtr).final * group->numSD[NEUTRON]
	       + nStep;
    numOneBits = countBits(group->bitArray,newPos);
    (*pnPtr).final = group->startAmp + newPos - numOneBits;
    numElem++;
    pnPtr++;
  
    checkMemSize--;

    if(checkMemSize <= 0) { 
      errorMemoryOverflow(grBas, "addTwoPart_ppConfig()");
    }
  
  } // end for-loop

  return numElem;

} // End function  addTwoPart_ppConfig()

  /*
  ** The function
  **    countBits()
  ** counts and return number
  ** of bits set in a bitArray[]
  */

static ULL countBits(BARR_ELTYPE *bitArray, ULL size) 
{
  ULL  bit, count;

  count = 0;
  for(bit = 0; bit < size; bit++) {
    if(BARR_TEST(bitArray,bit)) count++;
  }
  return count;

} // End: function  countBits()

   /*
   ** The function
   **     addTwoPart_nnConfig()
   ** converts all nondiag ppMatrix elements from
   ** current init |nStep> to a pn-scale and add 
   ** them to the total list of 
   ** pnNondiagMatrElemts.
   ** matrPtr points to start position of current
   ** nondiag ppMatrix elements.- topBit terminates 
   */

static int addTwoPart_nnConfig(ULL checkMemSize,GR_BAS *grBas,
                                    GROUP *group,int sdZ_valN, 
                             STORE_F *pnPtr, STORE_F *nn_nStepPtr)
{
  int     numElem,pStep, newN_N, numOneBits;
  ULL     topBit,newN_sd, newPos;

  topBit = ULL_ONE << (8 * sizeof(ULL) - 1);
  pStep = group->currSD_Z; 

  numElem = 0;

  for( ; ((*nn_nStepPtr).final != topBit);nn_nStepPtr++) {

    // check oscillator exitation for |newZ_sd,newN_sd>
 
    newN_sd = *(group->SD[NEUTRON] + (*nn_nStepPtr).final);
    newN_N  = valOscN(NEUTRON,grBas,newN_sd);

    if((sdZ_valN + newN_N) >  grBas->maxOscN) continue;

    (*pnPtr).value = (*nn_nStepPtr).value;
    newPos = pStep*group->numSD[NEUTRON] 
	    + (*nn_nStepPtr).final; 

    numOneBits = countBits(group->bitArray,newPos);
    (*pnPtr).final = group->startAmp + newPos - numOneBits;
    numElem++;
    pnPtr++;
    checkMemSize--;

    if(checkMemSize <= 0)  {
      errorMemoryOverflow(grBas,"addTwoPart_nnConfig()");
    }
  } // end for-loop

  return numElem;
    
} // End function  addTwoPart_nnConfig()

   /*
   ** The function
   **     addTwoPart_pnConfig()
   ** calculates all nondiag pn matrix elements
   ** H(PN)*|init_pn> ---> value * |final_pn>
   ** within a single proton/neutron group
   */

static int addTwoPart_pnConfig
           (int typeCM, ULL checkMemSize,int nRankStep,ULL sdZ, ULL sdN,
                STORE_F *pnPtr,GR_BAS *grBas,GROUP *initGroup,
                             GROUP *finalGroup,PN_INT **table)
{
  int        pStep, numElem, k, *k_ptr, l, *l_ptr,
             phase, count, sumN, numOneBits;
  ULL        sdZ_k, sdN_l, new_sdZ, new_sdN, high, low,
             numZ, numN,newPos;
  MASK       maskZ, maskN;
  PN_INT     **tab_p, *ij_ptr;

  pStep   = initGroup->currSD_Z;
  maskZ   = grBas->spBas[PROTON]->mask;
  maskN   = grBas->spBas[NEUTRON]->mask;
  numElem = 0;

  k     = grBas->spBas[PROTON]->part;
  k_ptr = initGroup->occ[PROTON] + (k + 1)*pStep;
  do {  // proton k loop
    tab_p  = table +(*(k_ptr))*grBas->spBas[NEUTRON]->numm_orb;
    sdZ_k = sdZ ^ (ULL_ONE << (*(k_ptr++)));
    l     = grBas->spBas[NEUTRON]->part;
    l_ptr = initGroup->occ[NEUTRON] + (l + 1)*nRankStep;
    do   {  // neutron l loop
      sdN_l = sdN ^ (ULL_ONE << (*l_ptr));
      if((ij_ptr = *(tab_p + (*(l_ptr++)))) == NULL) continue;
      do { // loop all p-n matrElem 
	    /* 
	    ** Pauli principle: 
	    ** proton orbit i and  neutron orbit j 
	    ** must both be empty. 
	    */ 
	while(   ij_ptr->one[0] 
	      && (  (ij_ptr->one[0] & sdZ_k)
		  ||(ij_ptr->one[1] & sdN_l)))   {
	  ij_ptr++;
	}
	if(!ij_ptr->one[0]) break; // end of ij-pairs 
	new_sdZ = (sdZ_k ^ ij_ptr->one[0]);

	// check oscillator exitation for |new_sdZ>

	sumN = valOscN(PROTON,grBas, new_sdZ);
	if(sumN <= grBas->maxOscN_Z) { 

	     /* 
	     ** mask check that |new_sdZ> has allowed number 
	     ** of particles in all j-orbits - only active if
	     **  0 < allowed number < 2*j + 1
             */
	  if(maskZ.num > 0) { // check masking
	    for(count = 0; count < maskZ.num; count++) {
	      high = new_sdZ & maskZ.list[count];
	      phase = 0;
	      while(high) {
		high &= high - 1;
		phase++;
	      }
	      if(  (phase > maskZ.lim[2*count + 1])
		 ||(phase < maskZ.lim[2*count])) {
		break; // |new_sdZ> not accepted
	      }
	    }
	    if(count < maskZ.num) continue;
	  } // end check PROTON masking
		
	  new_sdN = (sdN_l ^ ij_ptr->one[1]);

	  // check oscillator exitation for |new_sdZ;new_sdN>

	  sumN += valOscN(NEUTRON,grBas, new_sdN);
	  if(sumN <= grBas->maxOscN) {

             /* 
	     ** check that |new_sdN> has allowed number of 
	     ** particles in all j-orbits - only active if
	     **  0 < allowed number < 2*j + 1
             */

	    if(maskN.num > 0) { // check masking
	      for(count = 0; count < maskN.num; count++) {
		high = new_sdN & maskN.list[count];
		phase = 0;
		while(high) {
		  high &= high - 1;
		  phase++;
		}
		if(  (phase > maskN.lim[2*count + 1])
		   ||(phase < maskN.lim[2*count])) {
		  break; // |new_sdN> not accepted
		}
	      }
	      if(count < maskN.num) continue;
	    } // end check NEUTRON masking
	    phase = +1; // proton  permutation phase
	    low = new_sdZ & ij_ptr->two[0];
	    while(low) {
	      low &= low - 1;
	      phase = -phase;
	    }

	    low  = 0;           
	    high = (ULL)finalGroup->numSD[PROTON];
	    while (1) {
	      numZ = (low + high) >> 1;
	      if(     new_sdZ < finalGroup->SD[PROTON][numZ]) high = numZ-1;
	      else if(new_sdZ > finalGroup->SD[PROTON][numZ]) low  = numZ+1;
	      else                                            break;
	    }

	    // add neutron  permutation phase

	    low = new_sdN & ij_ptr->two[1];
	    while(low) {
	      low &= low - 1;
	      phase = -phase;
	    }

	    low  = 0;
	    high = (ULL)finalGroup->numSD[NEUTRON];
	    while (1) {
	      numN = (low + high) >> 1;
	      if(     new_sdN < finalGroup->SD[NEUTRON][numN]) high = numN - 1;
	      else if(new_sdN > finalGroup->SD[NEUTRON][numN]) low  = numN + 1;
	      else                                       break;
	    }
	    pnPtr->value = (float)(phase * ij_ptr->val[typeCM]);
	    newPos = numZ*finalGroup->numSD[NEUTRON] + numN;
	    if(BARR_TEST(finalGroup->bitArray,newPos)) continue;
	    numOneBits = countBits(finalGroup->bitArray,newPos);
	    pnPtr->final = finalGroup->startAmp + newPos - numOneBits;
	    pnPtr++;
	    numElem++;
	    checkMemSize--;

	    if(checkMemSize <= 0)  {
	      errorMemoryOverflow(grBas,"addTwoPart_pnConfig()");
	    }
	  }  // end sumN_N test
	}   //  end sumN_Z test

      } while((++ij_ptr)->one[0]); //end loop matrElem for fixed kl-pairs

    } while(--l);// end l-loop for neutrons
  } while(--k); // end k-loop for protons

  return numElem;

} // End: function   addTwoPart_pnConfig()

    /*
    ** The function                         
    **        int pnSortNondiag()                  
    ** is a utility function for the library function qsort() in order to
    ** sort nondiagSD matrix elements of type STORE_F store[]
    ** after increasing store[].final.
    */

static int pnSortNondiag(const STORE_F *one, const STORE_F *two)
{
  if(one->final > two->final)       return +1;
  else  if(one->final < two->final) return -1;
  else                              return  0;
} // End: function pnSortNondiag()

      /*
      ** The function 
      **      pnCompressNondiagSD_Elem()
      ** takes num calculated nondiag SD elem which are sorted
      ** after increasing store[].final and add together all
      ** elements  with same store[].final
      ** Then all element store{].value < MATRIX_LIMIT
      ** are removed.  .
      ** Note: The function must have num > 1;
      ** The final number of elements is returned.
      */

static int pnCompressNondiagSD_Elem(int num, STORE_F *store)
{
  int  lim, k, l;
   
   lim = num - 1;
   k   = 0;
   do {
     for( ; k < lim; k++) {
       if(store[k].final == store[k + 1].final) {
	 store[k].value += store[k + 1].value;
	 for(l = k + 1; l < lim; l++) {
	   store[l].value = store[l + 1].value;
	   store[l].final = store[l + 1].final;
	 } // end loop l
	 lim--;
	 break;
       }  //end if-test
     }  // end loop k
   } while(k < lim);

   lim++;
   k = 0;
   do {
     for( ; k < lim; k++) {
       if(fabs(store[k].value) > MATRIX_LIMIT) continue;
       store[k].value = store[k + 1].value;
       store[k].final = store[k + 1].final;
       for(l = k + 1; l < lim; l++) {
	 store[l].value = store[l + 1].value;
	 store[l].final = store[l + 1].final;
       } // end loop l
       lim--;
       break;
     }  // end loop k
   } while(k < lim);

   return lim;

} // End: function pnCompress_nondiag_SD_elem()

  /*
  ** The function
  **  printMatrElemData()
  ** collect info about precalculation of <SD'|OP|SD>
  ** and send data to outputData() 
  */

static void printMatrElemData(int typeCalc,GR_BAS *grBas,TID wall_time, TID cpu_time)
{ 
  int   data[9];
  
  data[0]  = (int)grBas->totMemNondiagElem;
  data[1]  = 0;                     //sdStore->currSD;
  data[2]  = grBas->totNumSD_ZN;
  data[3]  = (int) wall_time.hour;
  data[4]  = (int) wall_time.min;
  data[5]  = (int) wall_time.sec;
  data[6]  = (int) cpu_time.hour;
  data[7]  = (int) cpu_time.min;
  data[8]  = (int) cpu_time.sec;
 
  outputData(typeCalc,grBas,data);

} // End: function printMatrElemData()

  /*
  ** The function 
  **      outputData()
  ** send  info about precalculation of
  ** <SD'|OP|SD> to output file
  */

static void outputData(int typeCalc,GR_BAS *grBas,int *data)
{
  char  text[1000], filename[ONE_LINE];;
  FILE  *file_ptr;

  switch(typeCalc) {
       case  VEFF_CALC : sprintf(text,"%s","VEFF_INT");  
                         break;
       case  CM_CALC   : sprintf(text,"%s","CM_INT");  
                         break;
       case  ANG_CALC  : sprintf(text,"%s","ANG_INT");  
                         break;
  }

sprintf(filename,"%s%s",grBas->title,RESULT_OUTPUT);

if((file_ptr = fopen(filename,"a"))== NULL) {
  printf("\n\nRank%d: Error in function outputData():",Rank);
  printf("\nWrong file = %s for the output data file\n",filename);
  MPI_Abort(MPI_COMM_WORLD,Rank);
}
fprintf(file_ptr,"\n\nReserved permanent memory: numElem = %llu",
                                            grBas->permMemSize);
fprintf(file_ptr,"\n\nStored nondiag matrix elements in memory");

fprintf(file_ptr,"\nNumber of stored <SD'|%s|initSDZ,SDN>: %d",
                         text,data[0]); 
  fprintf(file_ptr,"\nUp to and including |initSD> = %d",data[1]); 
  fprintf(file_ptr," - Last |initSD> =  %d",data[2]);
  fprintf(file_ptr, "\nWall time used : %d hour %d min %d sec",
	  data[3],data[4],data[5]);
  fprintf(file_ptr, "\nCPU time used : %d hour %d min %d sec",
            	                      data[6],data[7],data[8]);
  fclose(file_ptr);
  
} // End: function outputData()

    /* 
    ** The function
    **       valOscN()
    ** calculates sumN for and return for
    ** |new_SD>.
    */

static int valOscN(int type,GR_BAS *grBas,ULL new_sd)
{
  int   num_part,par, orb, sumN;
  ULL   pos;
  MBAS  *mBasPtr;

  num_part     = grBas->spBas[type]->part; // initialization
  mBasPtr      = grBas->spBas[type]->mbas;

  sumN     = 0;
  for(par = 0, orb = 0, pos = ULL_ONE; par < num_part;
                                 par++, orb++, pos <<= 1)  {
    for(;!(new_sd & pos); orb++, pos <<= 1);
    sumN += mBasPtr[orb].N; 
  }

  return sumN;

} // End: function valOscN()

   /*
   ** The function
   **  errorMemoryOverflow()
   ** prints an error message to 
   ** program output file if no 
   ** more memory is available
   ** for data storage 
   */

static void errorMemoryOverflow(GR_BAS *grBas, char *func)
  {
    char   filename[ONE_LINE];
    FILE   *filePtr;

    sprintf(filename,"%s%s",grBas->title,RESULT_OUTPUT);
      if((filePtr = fopen(filename,"a"))== NULL) {
	printf("\n\nRank%d:Error in  errorMemoryOverflow()",Rank);
	printf("\nWrong file = %s for the output data file\n\n",filename);
	MPI_Abort(MPI_COMM_WORLD,Rank);
      }
      fprintf(filePtr,"\n\nRank%d:Error in function %s:",Rank,func );
      fprintf(filePtr,"\nNo more memory to store nondiag matr elem");
      fprintf(filePtr,"\nTotal permMemSize = %llu\n\n", grBas->permMemSize);
      MPI_Abort(MPI_COMM_WORLD,Rank);

  }  // End: function errorMemoryOverflow()
