/*******************  The module shell-model-input.c  ******************/

#include "shell-model.h"

     /*
     ** The module entrance function 
     **   void input_process(char *indata, SHELL *model)
     ** reads input data from file, store the data in SHELL model
     ** and checks that necessary data are available.
     ** All data are copied to output file
     */ 
                // **** input data definitions ***  

   static  int  const number_of_data = 14;  // number of data lines in input data file    
   static  char const *identity[] =         // list of text strings in input data file   
               {"The number of particle j-orbits",
                "Orbit",
                "The particle number",
                "Total angular momentum J is (even, odd)",
                "Twice total projection of angular momentum",
                "Total parity (+, -)",
                "Input V-effective matrix element",
                "Maximum Lanczos iterations (<1000)",
                "Wanted number of converged eigenstates",
                "Type of calculation process",
                "Shell model calculation of identical particle system -- title",
                "Center_of_Mass matrix elements included(yes = 1, no = 0)",
                "Type of effective interaction",
	        "If(yes): Filename Center_of_Mass matrix elements"};

static ULL input_data(char *in_data, SHELL *model);
      /*
      ** reads all proton-neutron data for the shell model problem.
      ** The function returns the data in SHELL *model
      */

static void write_input_data(char *in_data, char *out_data);
      /*
      ** write to output file all essential data in input data file
      */

static void check_input_data(ULL control, SHELL *model);
      /*
      ** analyzes all input data and checks that 
      ** necessary input data are submitted
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
     ** The function 
     **         input_process()
     ** reads input data from file, store the data in SHELL model
     ** and checks that necessary data are available
     */

void input_process(char *input_file, SHELL *model)
{ 
  ULL
         control;

        // read input data from file

  control = input_data(input_file, model);

        // write input data from input file to output file

  write_input_data(input_file, model->title);
 
  check_input_data(control, model);           // check and analyze input data

} // End: function input_process()   

    /*
    ** The function                              
    **             input_data()                
    ** reads  all proton-neutron data for the shell model problem.
    ** The function returns the data in SHELL *model..
    */

static ULL input_data(char *input_file, SHELL *model)
{
  char    *func = {"input_data(): "}, 
          text_string[ONE_LINE], data_string[ONE_LINE];
  int     i, count_orb, countOverlap;
  ULL     control = 0;  
  FILE    *data_ptr;

  model->part     = 0;  // initialization    
  model->numj_orb = 0;
  count_orb       = 0;
  countOverlap    = 0;

  if( (data_ptr = fopen(input_file,"r")) == NULL_PTR) {
    printf("\n\nError in function input_data()");
    printf("\nWrong file = %s for the input data\n", input_file);
    exit(1);
  }

  for( ; ; )   {                   //read data from input file
    if(!read_text_string(data_ptr, text_string)) break;
    for(i = 0; i < number_of_data; i++)  {
      if(!(strcmp(text_string, identity[i]))) break;
    }
    switch (i)  {
    case 0:                // The number of particle j-orbits
      if(!read_int_number(data_ptr, 1, &model->numj_orb)) break;
      control |= ULL_ONE;

               // memory for single-particle orbits in a JBAS structure

      model->jbas = MALLOC(model->numj_orb + 1, JBAS, func, "data.jbas");
      break;
    case 1:                                        // j_orbit
      if(!model->numj_orb)   {
	printf("\n\nError in function input_data():");
	printf("\n Number of single-particle orbit = ZERO\n");
	exit(1);
      }
      if(count_orb > model->numj_orb)   {
	printf("\n\nError in function input_data():");
	printf("\nNumber of single-particle orbits = %d", model->numj_orb);
	printf("\nTried to read j_orbit number %d\n",count_orb);
	exit(1);
      }
      if(!read_orbit_number(data_ptr, count_orb++, model->jbas)) break;
      control |= (ULL_ONE << 1);
      break;
    case 2:                // The particle number
      if(!read_int_number(data_ptr, 1, &model->part)) break;
      control |= (ULL_ONE << 2);
      break;
    case 3:                // Total angular momentum J is (even, odd)
      if(!read_data_string(data_ptr, data_string)) break;
      if(!strcmp(data_string, "even")) model->J_type = EVEN;
      else if(!strcmp(data_string, "odd")) model->J_type = ODD; 
      else  {
	printf("\n\nError in function input_data()");
	printf("\nError in init j_type (even or odd) = %s\n", data_string);
	exit(1);
      }
      control |= (ULL_ONE << 3);         
      break;
    case 4:                // Twice total projection of angular momentum
      if(!read_int_number(data_ptr, 1, &model->MJ)) break;
      
             //  test relation between particle number and total m-value
 
      if(MOD(model->part,2) != MOD(model->MJ,2))  {
	printf("\n\nError in function read_input_data()");
	printf("\n\n Wrong value of total MJ value !!!!!!");
	printf("\nParticle number = %d  Total MJ = %d\n\n", model->part, model->MJ);
	exit(1);
      }
      control |= (ULL_ONE << 4);
      break;
    case 5:                // Total parity (+, -)
      if(!read_data_string(data_ptr, data_string)) break;
      if((data_string[0] != '+') && (data_string[0] !='-'))   {
	printf("\n\nError in function input_data():");
	printf("\nThe total parity is not + or -\n");
	exit(1);
      }
      model->P = (data_string[0] == '+') ? +1 : -1;
      control |= (ULL_ONE << 5);
      break;
    case 6:                 // Input two-particle v_effective in J-scheme
      if(!read_data_string(data_ptr, model->veffMatrElem)) break;
      control |= (ULL_ONE << 6);
      break;
    case 7:                // Maximum dimension of the energy matrix
      if(!read_int_number(data_ptr, 1, &model->max_iterations)) break;
      control |= (ULL_ONE << 7);
      break;
    case 8:                // Wanted number of converged eigenstates
      if(!read_int_number(data_ptr, 1, &model->states)) break;
      control |= (ULL_ONE << 8);
      break;
    case 9:                 // Type of calculation process
      if(!read_data_string(data_ptr, model->type_calc)) break;
      control |= (ULL_ONE << 9);
      break;
    case 10:                         // calculation title
      if(!read_data_string(data_ptr, data_string)) break;
      if(strlen(data_string) == 0) {
	strcpy(model->title,"lanczo");    // no calculation title
      }
      else   {
	strcpy(model->title, data_string);    // calculation title
      }
      control |= (ULL_ONE << 10);
      break;
    case 11:  // center_of_Mass included ?  
      if(!read_int_number(data_ptr, 1, &model->calc_CM)) break;
      control |= (ULL_ONE  << 11);
      break;
    case 12:  // Type of effective interaction
      if(!read_int_number(data_ptr, 1, &model->typeVeff)) break;
      control |= (ULL_ONE << 12);
      break;
    case 13:                 // Input two-particle v_effective in J-scheme
      if(!read_data_string(data_ptr, model-> centerOfMass)) break;
      control |= (ULL_ONE << 13);
      break;
    default:
      break;
    } // end of switch
  } // end  input data

  fclose(data_ptr);  // close input data file

  if(count_orb < model->numj_orb)   {
    printf("\n\nError in function input_data():");
    printf("\nNumber of single-particle orbits = %d", model->numj_orb);
    printf("\nOnly %d orbits found in input data file\n",count_orb);
    exit(1);
  }
        // Re-order j_obits in data.j_bas[] after decreasing j-values.

  qsort(model->jbas,(size_t) model->numj_orb,(size_t)sizeof(JBAS),
                  (int(*)(const void *, const void *))j_comp);

      // determine the type of interaction

  return  control;   // return  control parameter

} // End: function input_data()   

     /*
     ** The function
     **              write_input_data()
     ** write to output file all essential data in input data file
     */

static void write_input_data(char *in_data, char *title)
{
  char
          text_string[ONE_LINE], *ptr, filename[ONE_LINE];
  int
          k;
  FILE
          *in_data_ptr, *out_data_ptr;

  if( (in_data_ptr = fopen(in_data,"r")) == NULL) { // open in_data file
    printf("\n\nError in function write_input_data()");
    printf("\nWrong file = %s for the input data\n",in_data);
    exit(1);
  }
  rewind(in_data_ptr);

  sprintf(filename,"%s%s",title,RESULT_OUTPUT);
  if((out_data_ptr = fopen(filename,"w")) == NULL) { // Open output data file
    printf("\n\nError in function write_input_data()");
    printf("\nWrong file = %s for the output data\n", filename);
    exit(1);
  }
  while(!feof(in_data_ptr))  {             // copy input file to output file
    fgets(text_string,100, in_data_ptr);                    // read a line  
    for(k = 0, ptr = text_string;   (*ptr != '*') && (*ptr != '\n')
                                   && ( k < ONE_LINE); k++, ptr++);
    if(*(++ptr) != '*') fputs(text_string, out_data_ptr);  // print a line
  } // All data are copied   
   
  fprintf(out_data_ptr,"\n                            -----  Results  ----");

  fclose(in_data_ptr);                   // close the input  data file
  fclose(out_data_ptr);                  //  close the output data file

}  // End: function write_input_data()

     /*
     ** The function 
     **            check_input_data()
     ** analyzes all input data and checks that 
     ** necessary input data are submitted
     */

static void check_input_data(ULL control, SHELL *model)
{
  char
          filename[ONE_LINE];
  int
          k;
  UL
          test, remove = 0, pos, missing_data;
  FILE
          *file_ptr;

         // initialize number_of_data one-bits - one for each input data

  test = (ULL_ONE << number_of_data) - 1;

  if(   strcmp(model->type_calc,"random-start") 
     || strcmp(model->type_calc,"random-continue")) {
    remove +=  (ULL_ONE << 10) + (ULL_ONE << 11)
             + (ULL_ONE << 12) + (ULL_ONE << 13);
  }
       // Test for type of effective interaction

  if((control & (ULL_ONE << 6)) | (control & (ULL_ONE << 16))) {
    remove +=  (ULL_ONE << 6) + (ULL_ONE << 16);
  }
  test ^= remove;
  pos = ULL_ONE;
  missing_data = 0;
  for(k = 0; k < number_of_data; k++, pos <<= 1) {
    if(!(test & pos)) continue;
    if(!(pos & control)) missing_data += pos;
  }

  if(missing_data) {

            // missing input data   

    printf("\n\nError in function check_input_data()");
    printf("\nmissing input data - see output file for information\n");

    strcpy(filename, model->title);
    strcat(filename,RESULT_OUTPUT);

    if((file_ptr = fopen(filename,"w")) == NULL) { 
      printf("\n\nError in function check_input_data()");
      printf("\nWrong file = %s for the output data\n", filename);
      exit(1);
    }
    fprintf(file_ptr,"\n\nError in function check_input_data():");
    fprintf(file_ptr,"\nThe following data are missing:");

    for(k = 0, pos = ULL_ONE; k < number_of_data; k++, pos <<= 1) {
      if(missing_data & pos) {
            fprintf(file_ptr,"\n%s", identity[k]);
      }
    }
    fclose(file_ptr);
    exit(1);
  } // test: data are missing

} // End: function check_input_data()

   /*
   ** The function
   **            read_orbit_number()
   ** reads and stores in JBAS jbas[]  all input
   ** numbers defining a spherical  j-orbit.
   ** If all data are valid the function returns TRUE,
   ** otherwise FALSE.
   */

static int read_orbit_number(FILE  *in_data, int orb_num, JBAS *jbas)
{
  int
            int_data[5];
  double
            dbl_val;

       // read five integer numbers   

  if(read_int_number(in_data, 5, int_data) == FALSE) {
    printf("\n\nError in function read_orbit_number();");
    printf("\nWrong integer data for orbit no %d\n", orb_num);
    return FALSE;
  }
  jbas[orb_num].osc      = int_data[0];
  jbas[orb_num].l        = int_data[1];
  jbas[orb_num].j        = int_data[2];
  jbas[orb_num].min_part = int_data[3];
  jbas[orb_num].max_part = int_data[4];

         // read single-particle energy   

  if(read_float_number(in_data, 1, &dbl_val) == FALSE) {
    printf("\n\nError in function read_float_number();");
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
