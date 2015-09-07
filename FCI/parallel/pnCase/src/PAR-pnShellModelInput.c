/*******************  The module PAR-pnShellModelInput.c  ******************/

#include "PAR-pnShellModel.h"

     /*
     ** The module entrance function 
     **   void input_process(char *inputFile,SHELL *model,)
     ** reads input data from file, store the data in SHELL model
     ** and checks that necessary data are available.
     ** All data are copied to output file
     */ 

        // *** local function declarations ***  

static int inputCalcType(char *inputFile, SHELL *model);
     /*
     ** reads and stores in model 
     ** the calculation title.
     */

static int input_j_orbits(char *input_file, int typaCalc,SHELL *model);
    /*
    ** reads  all proton-neutron data for the shell model
    ** single-particle data included possible additional
    ** harmonic oscillator exitation..
    */

static int extended_j_orbits(char *inputFile, int typeCalc,SHELL *model);
  /*
  ** reads input data and calculates extended
  ** harmonic oscillator j_orbits
  */

static int primary_j_orbits(char *inputFile, SHELL *model);
  /*
  ** reads input data and calculates the totan
  ** number of extended and primary j-orbits
  */

static void sp_osc_extension(SHELL *model, int *num_primary_orb);
       /*
       ** adds single-particle j_orbits from 
       ** additional harmonic oscillator shells 
       */
static int input_rest_data(char *inputFile, int typeCalc, SHELL *model);
    /*
    ** reads rest proton-neutron data for the shell model problem.
    */

static void write_input_data(char *in_data, char *title);
     /*
     ** write to output file all essential data in input data file
     */
static int read_orbit_number(FILE *in_data, int orb_num, JBAS *jbas);
      /*
      ** reads and stores in JBAS jbas[]  all input numbers defining
      ** a spherical  j-orbit. If all data are valid the function
      ** returns TRUE, otherwise FALSE.
      */
static int j_comp(const JBAS *one, const JBAS *two);
      /*
      ** is a utility function for the library function qsort() in order to
      ** sort the single-particle orbits in jbas[] in function read_orbit() 
      ** after decreasing j-values.
      */

                // *** End: function declarations ***  

                // *** The function definitions  ***   

     /*
     ** The module entrance function 
     **   void input_process(char *inputFile,SHELL *model,)
     ** reads input data from file, store the data in SHELL model
     ** and checks that necessary data are available.
     ** All data are copied to output file
     */ 

int input_process(char *inputFile, SHELL *model)
 {    
   int typeCalc, result;

   // read Shell model calculation of proton/neutron system? -- title 

   typeCalc = inputCalcType(inputFile, model);


  // write input data from input file to output file

  if(Rank == MASTER) {
    write_input_data(inputFile, model->title);
  
  }

  // read input data for single-particle orbits

  result = input_j_orbits(inputFile,typeCalc,model);

  // read remaining input data

  result = input_rest_data(inputFile,typeCalc, model);

  return result;

} // End: function input_process()   

     /*
     ** The function
     **    inputCalcType()
     ** reads and stores in model 
     ** the calculation title.
     */

static int inputCalcType(char *inputFile, SHELL *model)
{
  int  const  number_of_data = 2;
    int       typeCalc = -1, i;

  // list of text strings from input data file

  char const *identity[] = // list of text strings from input_file   
              {"Shell model calculation of proton/neutron systemA -- title",
	       "Shell model calculation of proton/neutron systemB -- title"};
  char    *func = {" inputCalcType(): "}, 
          data_string[ONE_LINE],text_string[ONE_LINE];
  FILE    *dataPtr;

  if( (dataPtr = fopen(inputFile,"r")) == NULL_PTR) {
    printf("\n\nRank%d:Rank: Error in function %s",Rank,func);
    printf("\nWrong file = %s for the input data\n",inputFile);
    MPI_Abort(MPI_COMM_WORLD,Rank);

  }
  for( ; ; ) {  // read data from inputFile
    if(!read_text_string(dataPtr, text_string)) break;
  
    for(i = 0; i < number_of_data; i++) {
      if(!(strcmp(text_string, identity[i]))) break;
    }
    switch (i)  {
    case 0://Shell model calculation of proton/neutron systemA -- title
      read_data_string(dataPtr, data_string);
      if(strlen(data_string) == 0) {
	sprintf(model->title,"%s","lanczo");
      }
      else   {
	sprintf(model->title,"%s%s%d",data_string,"RANK",Rank);
      }
      model->calc_CM;
      typeCalc = 0;

      break;
    case 1://Shell model calculation of proton/neutron systemB -- title
      read_data_string(dataPtr, data_string);
      if(strlen(data_string) == 0) {
	sprintf(model->title,"%s","lanczo");
      }
      else   {
	//	sprintf(model->title,"%s%s%d",data_string,"RANK",Rank);
	sprintf(model->title,"%s%s",data_string,"XXX");
      }
      typeCalc = 1;
      break;
    default:
      break;
    } // end of switch
  } // end read input data 

  fclose(dataPtr);

  return typeCalc;
  
  } // End function  inputCalcType()

    /*
    ** The function                              
    **     input_j_orbits()                
    ** reads  all proton-neutron data for the shell model
    ** single-particle data included possible additional
    ** harmonic oscillator exitation..
    */

static int input_j_orbits(char *inputFile, int typeCalc, SHELL *model)
{
  int   result;

   /*
   ** Two sets of possible spherical j-orbits:
   ** 1. extended harmonic oscillator j-orbits
   ** 2. Primary j-orbits
   */

  result = extended_j_orbits(inputFile,typeCalc,model);
  result = primary_j_orbits(inputFile,model);

  return result;

} // End:function input_j_orbits()

  /*
  ** The function
  **      extended_j_orbits()
  ** reads input data and calculates extended
  ** harmonic oscillator j_orbits
  */

static int extended_j_orbits(char *inputFile, int typeCalc, SHELL *model)
{
  int  const number_of_data = 5;

  // list of text strings from input data file

  char const *identity[] = // list of text strings from input_file   
              {"Proton N_Z_low",
               "PROTON N_Z_high",
               "NEUTRON N_N_low",
               "NEUTRON N_N_high",
               "The harmonic oscillator parameter hBar_omega(MeV)"};

  char    *func = {"extended_j_orbit(): "}, 
          text_string[ONE_LINE];
  int     i,type, N, n, l;
  ULL     testData, control = 0;
  FILE    *dataPtr;

  if( (dataPtr = fopen(inputFile,"r")) == NULL_PTR) {
    printf("\n\nRank%d:Error in function %s",Rank,func);
    printf("\nWrong file = %s for the input data\n",inputFile);
   MPI_Abort(MPI_COMM_WORLD,Rank);

  }

  model->N_low[0]    = 0;     // default values
  model->N_low[1]    = 0;
  model->N_high[0]   = 0;
  model->N_high[1]   = 0;
  model->numj_orb[0] = 0;
  model->numj_orb[1] = 0;

  for( ; ; ) {  // read data from inputFile
    if(!read_text_string(dataPtr, text_string)) break;
    for(i = 0; i < number_of_data; i++) {
      if(!(strcmp(text_string, identity[i]))) break;
    }
    switch (i)  {
     case 0: // Proton N_Z_low 
      if(!read_int_number(dataPtr, 1, &model->N_low[0])) break;
      control |= ULL_ONE;
      break;
    case 1: // Proton N_Z_high 
      if(!read_int_number(dataPtr, 1, &model->N_high[0])) break;
      control |= (ULL_ONE << 1);
      break;
    case 2: // NEUTRON N_Z_low 
      if(!read_int_number(dataPtr, 1, &model->N_low[1])) break;
      control |= (ULL_ONE << 2);
      break;
    case 3: // NEUTRON N_Z_high 
      if(!read_int_number(dataPtr, 1, &model->N_high[1])) break;
      control |= (ULL_ONE << 3);
      break;
    case 4: // The harmonic oscillator parameter hBar_omega(MeV)
      if(!read_float_number(dataPtr, 1, &model->oscEnergy)) break;
      control |= (ULL_ONE << 4);
      break;
    default:
      break;
    } // end of switch
  }

  fclose(dataPtr);

  if(typeCalc == 0) return  YES; // calculation systemA 

  // extended calculation systemB:

  // check for correct input data

  if(model->N_low[0] == 0) control |= (ULL_ONE << 1);
  if(model->N_low[1] == 0) control |= (ULL_ONE << 3);
  if((model->N_low[0] + model->N_low[1]) == 0) control |= (ULL_ONE << 4);

  testData = (ULL_ONE << number_of_data) - 1;

  if(testData != control) {
    char   filename[ONE_LINE];
    int    k;
    ULL    pos;
    FILE  *filePtr;

    if(Rank == MASTER) {

      sprintf(filename,"%s%s",model->title,RESULT_OUTPUT);

      if((filePtr = fopen(filename,"a")) == NULL_PTR) {
	printf("\n\nRank%d:Error in function %s",Rank,func);
	printf("\nWrong file = %s for the input data\n",filename);
	printf("from error in input data\n\n");
	MPI_Abort(MPI_COMM_WORLD,Rank);
      }

      fprintf(filePtr,"\n\n\n Missing input data");

      for(k = 0, pos = ULL_ONE; k < number_of_data; k++, pos <<= 1) {
	if(!(testData & pos)) {
	  fprintf(filePtr,"\n%s",identity[k]);
	}
      }
      fclose(filePtr);

    } // end MASTER

    return NO;  // error in input data

  } // end check input data

      /* 
      ** Calculate number of extended j-orbits from 
      ** addtional harmonic oscillator shells
      */

  for(type = 0; type <= 1; type++) {
    if(model->N_low[type] == 0) continue;
    for(N = model->N_low[type];N <= model->N_high[type];N++) {
      for(n = 0; n <= N/2; n++) {
	l = N - 2*n;
	model->numj_orb[type] += ((l == 0) ? 1 : 2);
      } // end n loop
    } // end N loop
  } // P/N loop

  return YES;

} // End: function extended_j_orbits()

  /*
  ** The function
  **      primary_j_orbits()
  ** reads input data and calculates the totan
  ** number of extended and primary j-orbits
  */

static int primary_j_orbits(char *inputFile, SHELL *model)
{
  int  const number_of_data = 4;

  // list of text strings from input data file

  char const *identity[] = // list of text strings from input_file   
              {"The number of proton particle j-orbits",
               "Orbit_Z",
               "The number of neutron particle j-orbits",
               "Orbit_N"};
 
  char    *func = {"primary_j_orbits(): "}, 
          text_string[ONE_LINE];
  int     i, count_orb[2],num_primary_orb[2],num;
  ULL     testData, control = 0;
  FILE    *dataPtr;

  if((dataPtr = fopen(inputFile,"r")) == NULL_PTR) {
    printf("\n\nRank%d:Error in function %s",Rank,func);
    printf("\nWrong file = %s for the input data\n",inputFile);
   MPI_Abort(MPI_COMM_WORLD,Rank);
  }

  num_primary_orb[0] = 0;
  num_primary_orb[1] = 0;
  count_orb[0]       = 0;
  count_orb[1]       = 0;

  for( ; ; ) {  // read data from input file
    if(!read_text_string(dataPtr, text_string)) break;
    for(i = 0; i < number_of_data; i++) {
      if(!(strcmp(text_string, identity[i]))) break;
    }
    switch (i)  {
    case 0:   // The number of primary proton particle j-orbits
      if(!read_int_number(dataPtr, 1, &num_primary_orb[0])) break;

      // memory for total proton orbits in a JBAS structure

      num = model->numj_orb[0] + num_primary_orb[0];
      if(num == 0)  {
	printf("\n\nRank%d:Error in function %s",Rank,func);
	printf("\n Number of proton single-particle orbit = ZERO\n");
	MPI_Abort(MPI_COMM_WORLD,Rank);
      }
      model->jbas[0] = MALLOC(num,JBAS,func,"model->jbas[0]");
      control |= ULL_ONE;
      break;
    case 1:  // Orbit_Z
      if(count_orb[0] > num_primary_orb[0])   {
	printf("\n\nRank%d: Error in function %s:",Rank,func);
	printf("\nNumber of proton orbits = %d", num_primary_orb[0]);
	printf("\nTried to read j_orbit number %d\n",count_orb[0]);
	MPI_Abort(MPI_COMM_WORLD,Rank);
      }
      if(!read_orbit_number(dataPtr,count_orb[0]++, model->jbas[0])) break;
      control |= (ULL_ONE << 1);
      break;
     case 2:   // The number of neutron particle j-orbits
      if(!read_int_number(dataPtr, 1, &num_primary_orb[1])) break;

      // memory for neutron orbits in a JBAS structure

      num = model->numj_orb[1] + num_primary_orb[1];
      if(num == 0)   {
	printf("\n\nRank%d: Error in function %s:",Rank,func);
	printf("\n Number of neutron orbit = ZERO\n");
	MPI_Abort(MPI_COMM_WORLD,Rank);
      }
      model->jbas[1] = MALLOC(num,JBAS,func,"model->jbas[1]");
      control |= (ULL_ONE<< 2);
      break;
    case 3:  // Orbit_N
      if(count_orb[1] > num_primary_orb[1]) {
	printf("\n\nRank%d: Error in function %s",Rank,func);
	printf("\nNumber of neutron orbits = %d",num_primary_orb[1]);
	printf("\nTried to read j_orbit number %d\n",count_orb[1]);
	MPI_Abort(MPI_COMM_WORLD,Rank);
      }
      if(!read_orbit_number(dataPtr,count_orb[1]++,model->jbas[1])) break;
      control |= (ULL_ONE << 3);
      break;
    default:
      break;
   } // end of switch
 
  } // end  input data

  fclose(dataPtr);  // close input data file

  // Input data control test 

  testData = (ULL_ONE << number_of_data) - 1; 

  if(testData != control) {
    char  filename[ONE_LINE];
    int   k;
    ULL   pos;
    FILE  *filePtr;

    if(Rank == MASTER) {
      sprintf(filename,"%s%s",model->title,RESULT_OUTPUT);

      if((filePtr = fopen(filename,"a")) == NULL) { 
	fprintf(filePtr,"\n\nRank%d: Error in function %s",Rank,func);
	fprintf(filePtr,"\nWrong file = %s for the output data\n", filename);
	printf("from error in input data\n\n");
      }

      fprintf(filePtr,"\n\n\n Missing input data");

      for(k = 0, pos = ULL_ONE; k < number_of_data; k++, pos <<= 1) {
	if(!(testData & pos)) {
	  fprintf(filePtr,"\n%s",identity[k]);
	}
      }
      fclose(filePtr);

    } // end MASTER

    return NO;  // error in input data

  }  // end control test 

  /*
  ** Add to model->jbas[] the single-particle
  ** j_orb extensions due to additional 
  ** harmonic oscillator shells
  */

  model->numj_orb[0] += num_primary_orb[0]; // total number of j_orbits
  model->numj_orb[1] += num_primary_orb[1];

  sp_osc_extension(model,num_primary_orb);

  // Re-order proton j_obits in data.j_bas[] after decreasing j-values.

  qsort(model->jbas[0],(size_t) model->numj_orb[0],(size_t)sizeof(JBAS),
                  (int(*)(const void *, const void *))j_comp);

  // Re-order neutron j_obits in data.j_bas[] after decreasing j-values.

  qsort(model->jbas[1],(size_t) model->numj_orb[1],(size_t)sizeof(JBAS),
                  (int(*)(const void *, const void *))j_comp);

  return YES;

} // End: function primary_j_orbits()

       /*
       ** The function
       **      sp_osc_extension()
       ** adds single-particle j_orbits from 
       ** additional harmonic oscillator shells 
       */

static void sp_osc_extension(SHELL *model, int *num_primary_orb)
  {
    int    type, N, n, l, orb_num;
    JBAS   *jBas;
  
    for(type = 0; type <= 1; type++) {
      if(model->N_low[type] == 0) continue;
      jBas = model->jbas[type];
      orb_num = num_primary_orb[type];
      for(N = model->N_low[type];N <= model->N_high[type];N++) {
	for(n = 0; n <= N/2; n++) {
	  l = N - 2*n;
	  jBas[orb_num].N        = N;
	  jBas[orb_num].osc      = n;
	  jBas[orb_num].l        = l;
	  jBas[orb_num].j        = 2*l + 1;
	  jBas[orb_num].min_part = 0;
	  jBas[orb_num].max_part = jBas[orb_num].j  + 1;
	  jBas[orb_num].e        = N * model->oscEnergy;
	  orb_num++;
		
	  if(l > 0) {
	    jBas[orb_num].N        = N;
	    jBas[orb_num].osc      = n;
	    jBas[orb_num].l        = l;
	    jBas[orb_num].j        = 2*l - 1;
	    jBas[orb_num].min_part = 0;
	    jBas[orb_num].max_part = jBas[orb_num].j - 1;
	    jBas[orb_num].e        = N * model->oscEnergy;
	    orb_num++;
	  }
	} // end n loop
      } // end N loop
      model->numj_orb[type] = orb_num;
    } // end P/N loop

  } // End: function  sp_osc_extension()

    /*
    ** The function                              
    **             input_rest_data()                
    ** reads rest proton-neutron data for the shell model problem.
    */

static int input_rest_data(char *inputFile, int typeCalc,SHELL *model)
{

  char    *func = {"primary_j_orbits(): "}, 
          text_string[ONE_LINE],data_string[ONE_LINE];
  int     i, value;
  ULL     testData, control = 0;  
  FILE    *dataPtr;

  static  int  const number_of_data = 17;  // num data lines in input file    
  static  char const *identity[] =         // list text strings in input file   
              {"The proton number",
               "The neutron number",
	       "Twice total projection of angular momentum",
               "Total parity (+, -)", 
               "Type of effective interaction",
               "Input proton-proton   V-effective",
               "Input neutron-neutron V-effective",
               "Input proton-neutron  V-effective",
	       "Center-of-Mass matrix effects",
               "Type of calculation process",
               "Memory size to store nondiag <pn_SD'|OP|pn_SD> in MB (default 500 MB)",
               "Max file size for lanczos vectors    in Gb (default 10 GB)",
               "Max oscillator energy excitation for total P/N delN",
               "Maximum Lanczos iterations (<1000)",
               "Wanted number of converged eigenstates",
               "Max oscillator energy excitation for protons delN_Z",
               "Max oscillator energy excitation for neutrons delN_N"};


  model->file_lanczo  = LANC_FILES * 1000; // default file size in MByte

  if( (dataPtr = fopen(inputFile,"r")) == NULL_PTR) {
    printf("\n\nRank%d:Error in function %s",Rank,func);
    printf("\nWrong file = %s for the input data\n", inputFile);
   MPI_Abort(MPI_COMM_WORLD,Rank);
  }
  for( ; ; )   {                   // read data from input file
    if(!read_text_string(dataPtr, text_string)) break;
    i = 0;
    for(i = 0; i < number_of_data; i++) {
      if(!(strcmp(text_string, identity[i]))) break;
    }
    switch (i)  {
    case 0: // The proton number
      if(!read_int_number(dataPtr, 1, &model->part[0])) break;
      control |= ULL_ONE ;
      break;
    case 1: // The neutron number
      if(!read_int_number(dataPtr, 1, &model->part[1])) break;
      control |= (ULL_ONE << 1);
      break;
    case 2:  // Twice total projection of angular momentum
      if(!read_int_number(dataPtr, 1, &model->MJ)) break;
      
      //  test relation between particle number and total m-value
 
      if(MOD(model->part[0] + model->part[1],2) != MOD(model->MJ,2))  {
	printf("\n\nRank%d: Error in function read_input_data()",Rank);
	printf("\n\n Wrong value of total MJ value !!!!!!");
	printf("\nProton number = %d  Neutron number = %d Total MJ = %d\n\n",
                       model->part[0],  model->part[1], model->MJ);
	MPI_Abort(MPI_COMM_WORLD,Rank);
      }
      control |= (ULL_ONE << 2);
      break;
    case 3: // Total parity (+, -)
      if(!read_data_string(dataPtr, data_string)) break;
      if((data_string[0] != '+') && (data_string[0] !='-'))   {
	printf("\n\nRank%d: Error in function input_data():",Rank);
	printf("\nThe total parity is not + or -\n");
	MPI_Abort(MPI_COMM_WORLD,Rank);
      }
      model->P = (data_string[0] == '+') ? +1 : -1;
      control |= (ULL_ONE << 3);
      break;
    case  4:  // Type of effective interaction
      if(!read_int_number(dataPtr, 1, &model->typeVeff)) break;
      control |= (ULL_ONE << 4);
      break;
    case 5: // Input proton-proton   V-effective
      if(!read_data_string(dataPtr, model->veffMatrElem_pp)) break;
      control |= (ULL_ONE << 5);
      break;
    case 6: // Input neutron-neutron   V-effective
      if(!read_data_string(dataPtr, model->veffMatrElem_nn)) break;
      control |= (ULL_ONE << 6);
      break;
    case 7: // Input proton-neutron   V-effective
      if(!read_data_string(dataPtr, model->veffMatrElem_pn)) break;
      control |= (ULL_ONE << 7);
      break;
    case 8: // Center-of-MassM matrix elements(yes=1,no=0)
      if(!read_int_number(dataPtr, 1, &model->calc_CM)) break;
      control |= (ULL_ONE << 8);
      break;
    case 9: //Type of calculation process
      if(!read_int_number(dataPtr, 1, &model->typeProc)) break;

      control |= (ULL_ONE << 9);
      break;
    case 10://Memory size to store number of nondiag <pn_'SD|OP|pn_SD> in MB (default 500 MB)",
      model->permMemSize = (ULL)(500*((ULL)M_BYTES))
                          /sizeof(STORE_F); //default 500 MB
      if(read_int_number(dataPtr, 1, &value)) {
	if(value > 1) {
	  model->permMemSize = (ULL)(value/sizeof(STORE_F))*M_BYTES;
	}
      }
      control |= (ULL_ONE << 10);
      break;
    case 11: // Max. file size for lanczos vectors 
      model->file_lanczo = (ULL)(10*1000*((ULL)M_BYTES))
                          /sizeof(float);// 10 Gb 
      if(read_int_number(dataPtr, 1, &value)) {
	if(value > 10) {
	  model->file_lanczo = (ULL)(value*1000*((ULL)M_BYTES))
                              *sizeof(float);
	}
      }
      control |= (ULL_ONE << 11);
      break;
    case 12:  // Max oscillator energy excitation delN:
      if(!read_int_number(dataPtr, 1, &model->delN)) break;
      control |= (ULL_ONE << 12);
      break;
    case 13:  //Maximum Lanczos iterations (<1000)
      if(!read_int_number(dataPtr, 1, &model->max_iterations)) break;
      control |= (ULL_ONE << 13);
      break;
    case 14: // Wanted number of converged eigenstates
      if(!read_int_number(dataPtr, 1, &model->states)) break;
      control |= (ULL_ONE << 14);
      break;
    case 15: // Max oscillator energy excitation delN_Z:
      if(!read_int_number(dataPtr, 1, &model->delN_Z)) break;
      control |= (ULL_ONE << 15);
      break;
    case 16: // Max oscillator energy excitation delN_N:
      if(!read_int_number(dataPtr, 1, &model->delN_N)) break;
      control |= (ULL_ONE << 16);
      break;
    default:
      break;
   } // end of switch
 
  } // end  input data

  fclose(dataPtr);  // close input data file

  if(typeCalc == 0) { // calculation system B
    control |= (ULL_ONE << 8);
    control |= (ULL_ONE << 12);
    control |= (ULL_ONE << 15);
    control |= (ULL_ONE << 16);

    model->calc_CM = 0; 
    model->delN    = 100;   // just large values
    model->delN_Z  = 100;
    model->delN_N  = 100;
  }

  testData = (ULL_ONE << number_of_data) - 1; 

  if(testData != control) {
    char  filename[ONE_LINE];
    int   k;
    ULL   pos;
    FILE  *filePtr;

    if(Rank == MASTER) {
      sprintf(filename,"%s%s",model->title,RESULT_OUTPUT);

      if((filePtr = fopen(filename,"a")) == NULL) { 
	fprintf(filePtr,"\n\nRank%d: Error in function %s",Rank,func);
	fprintf(filePtr,"\nWrong file = %s for the output data\n\n", filename);
      }

      fprintf(filePtr,"\n\n\n Missing input data from function input_rest_data()");
      fprintf(filePtr,"\n testData = %lx  control = %lx",testData,control);

      for(k = 0, pos = ULL_ONE; k < number_of_data; k++, pos <<= 1) {
	if(!(testData & pos)) {
	  fprintf(filePtr,"\n%s",identity[k]);
	}
      }
      fclose(filePtr);

    } // end MASTER

    return NO; // error in input data

  }

  return YES;  // input data ok

} // End: function input_rest_data()   

     /*
     ** The function
     **              write_input_data()
     ** write to output file all essential data in input data file
     */

static void write_input_data(char *in_data, char *title)
{
  char    text_string[ONE_LINE],filename[ONE_LINE];
  FILE    *in_data_ptr, *out_data_ptr;

  if( (in_data_ptr = fopen(in_data,"r")) == NULL) { // open in_data file
    printf("\n\nRank%d: Error in function write_input_data()",Rank);
    printf("\nWrong file = %s for the input data\n\n",in_data);
   MPI_Abort(MPI_COMM_WORLD,Rank);
  }
  rewind(in_data_ptr);

  sprintf(filename,"%s%s",title,RESULT_OUTPUT);
  if((out_data_ptr = fopen(filename,"w")) == NULL) { // Open output data file
    printf("\n\nRank%d: Error in function write_input_data()",Rank);
    printf("\nWrong file = %s for the output data\n\n", filename);
   MPI_Abort(MPI_COMM_WORLD,Rank);
  }

  // copy input file to output file

  while(!feof(in_data_ptr))  { 
    fgets(text_string,100, in_data_ptr); // read a line
    fputs(text_string, out_data_ptr);    // print a line
  } // All data are copied   
   
  fprintf(out_data_ptr,"\n                            -----  Results  ----");

  fclose(in_data_ptr);                   // close the input  data file
  fclose(out_data_ptr);                  //  close the output data file

}  // End: function write_input_data()

   /*
   ** The function
   **            read_orbit_number()
   ** reads and stores in JBAS jbas[]  all input
   ** numbers defining a spherical  j-orbit.
   ** If all data are valid the function returns TRUE,
   ** otherwise FALSE.
   */

static int read_orbit_number(FILE *in_data, int orb_num, JBAS *jbas)
{
  int    int_data[6];
  double dbl_val;

       // read five integer numbers   

  if(read_int_number(in_data, 6, int_data) == FALSE) {
    printf("\n\nRank%d: Error in function read_orbit_number();",Rank);
    printf("\nWrong integer data for orbit no %d\n", orb_num);
    return FALSE;
  }
  jbas[orb_num].N        = int_data[0];
  jbas[orb_num].osc      = int_data[1];
  jbas[orb_num].l        = int_data[2];
  jbas[orb_num].j        = int_data[3];
  jbas[orb_num].min_part = int_data[4];
  jbas[orb_num].max_part = int_data[5];

         // read single-particle energy   

  if(read_float_number(in_data, 1, &dbl_val) == FALSE) {
    printf("\n\nRank%d: Error in function read_float_number();",Rank);
    printf("\nWrong floating point data for orbit no %d\n", orb_num);
    return FALSE;
  }
  jbas[orb_num].e = dbl_val;   // store the single-particle energy   

  return TRUE;

} // End: function  read_orbit_number()   

    /*
    ** The function                         
    **        int j_comp()                  
    ** is a utility function for the library function qsort() in order to
    ** sort the single-particle orbits in jbas[] in function read_orbit() 
    ** after decreasing j-values.
    */

static int j_comp(const JBAS *one, const JBAS *two)
{
  if(one->j > two->j)       return -1;
  else  if(one->j < two->j) return +1;
  else if(one->j == two->j)  {
    if(PARITY(one->l)) return +1;
    else               return -1;
  }
  else                       return  0;

} // End: function j_comp()   
