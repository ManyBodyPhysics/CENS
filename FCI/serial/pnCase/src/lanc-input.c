/*******************  The module lanc-input.c  ******************/

#include "shell.h"

     /*
     ** The module entrance function 
     **   void input_process(INPUT *data)
     ** reads input data from file and decides what
     ** type of shell model calculation to perform
     **    1. Identical particle shell model calculation
     **    2. Proton-neutron shell model calculation
     ** and sets up names for output- and storage-files. 
     ** All data are returned in struct INPUT data
     */

             /**** data definitions ****/

   static  int  const  number_of_data = 22;     /* number of data lines in input data file */ 
   static  char const  *identity[] =            /* list of text strings in input data file */
               {"The number of proton particle j-orbits",
                "Orbit_Z",
                "The number of neutron particle j-orbits",
                "Orbit_N",
                "The proton number",
                "The neutron number",
                "Total angular momentum J is (even, odd)",
                "Twice total projection of angular momentum",
                "Total parity (+, -)",
                "Input proton-proton v_effective in J-scheme",
                "Input neutron-neutron v_effective in J-scheme",
                "Input proton-neutron v_effective in J-scheme",
                "Maximum Lanczos iterations",
                "Wanted number of converged eigenstates",
                "Type of calculation process",
                "Total 2J value",
                "Init vector file",
                "Number of start vectors",
                "List of vector no",
                "Shell model calculation of proton/neutron system -- title",
                "Memory to store nondiag <SD(Z)SD(N)|OP|SD(Z)SD(N)>(in Kb)",
                "File size to store nondiag <SD(Z)SD(N)|OP|SD(Z)SD(N)>(in Mb)"};
 
		    /**** local function declarations ****/

static void input_data(INPUT *data);
      /*
      ** reads  all data for the shell model problem.
      ** The function returns the data in INPUT *data
      */

static int read_text_string(FILE *in_data, char *text_string);
      /*
      ** reads a character string which identify a specific line of input data.
      ** The string starts following the symbol "\*" and ends with colon ":"
      ** or semicolon ";". If the string  contains the character \, EOF
      ** or reads text_string "END_OF_INPUT_FILE", the function returns FALSE,
      ** otherwise TRUE.
      */

static int read_int_number(FILE *in_data, int number, int *vector);
      /*
      ** reads (number) integers from FILE *in_data stream and
      ** store them in vector[]. If too few integers are read
      ** the function returns  FALSE, otherwise TRUE.
      */

static int read_orbit_number(FILE *in_data, int orb_num, JBAS *jbas);
      /*
      ** reads and stores in JBAS jbas[]  all input numbers defining
      ** a spherical  j-orbit. If all data are valid the function
      ** returns TRUE, otherwise FALSE.
      */

static int read_float_number(FILE *in_data, int number, double *vector);
      /*
      ** reads (number) float from FILE *in_data stream,
      ** convert them to double and store them in vector[].
      ** If too few doubles are read the function returns 
      ** FALSE, otherwise TRUE.
      */

static int read_data_string(FILE *in_data, char *data_string);
      /*
      ** reads one piece of input data as a text_string containing no blanks. 
      ** Whitespaces in front of the data_string is removed. Data_string must
      ** not contain the character '\'. The function returns FALSE if no 
      ** proper data_string is found, otherwise TRUE.
      */

static int j_comp(const JBAS *one, const JBAS *two);
      /*
      ** is a utility function for the library function qsort() in order to
      ** sort the single-particle orbits in jbas[] in function read_orbit() 
      ** after decreasing j-values.
      */

static int particle_structure(int numZ, int numN);
      /*
      ** calculates and return type of shell model calculation
      ** with the following coding 
      **      type =  0 - protons only
      **              1 - neutrons only
      **              2 - one proton and and one neutron
      **              3 - one proton and several neutrons
      **              4 - several protons and one neutron
      **              5 - several protons and several neutrons
      */

static void write_input_data(INPUT *data);
      /*
      ** write to output file all essential data in input data file
      */

static void check_data(INPUT *data);
      /*
      ** analyzes all input data and checks that 
      ** necessary input data are submitted
      */

                /**** End: function declarations ****/

               /**** The function definitions  ****/ 
     /*
     ** The function 
     **         input_process()
     ** reads input data from file and decides what
     ** type of shell model calculation to perform
     **    1. Identical particle shell model calculation
     **    2. Proton-neutron shell model calculation
     */

void input_process(INPUT *data)
{
  input_data(data);                      /* read input data from file */

         /* identify proton/neutron structure */

  data->part_type = particle_structure(data->part[0], data->part[1]);
   
                         /* name output file for the result */

  strcpy(data->out_data,"result-");
  strcat(data->out_data, data->title);
  strcat(data->out_data, ".dat");

                 /* name eigenvector file in SD basis */

  strcpy(data->eigen_vec,"result-");
  strcat(data->eigen_vec,data->title);
  strcat(data->eigen_vec,"-eigenvectors.dat");

                 /* name data file for a restarting Lanczos process */

  strcpy(data->h_final, "result-");
  strcat(data->h_final,data->title);
  strcat(data->h_final,"_restarting.dat");
 
  write_input_data(data);     /* write input data from input file to output file */

  check_data(data);           /* check and analyze input data */

} /* End: function input_process() */

    /*
    ** The function                              
    **         input_data()                
    ** reads  all data for the shell model problem.
    ** The function returns the data in INPUT *data..
    */

static void input_data(INPUT *data)
{
  char
          *func = {"input_data(): "}, 
          text_string[ONE_LINE], data_string[ONE_LINE];
  int    
          i, number, count_orbZ = 0, count_orbN = 0;
  UL
          control= 0;                                    /* data check */
  FILE
          *data_ptr;

  data->part[0]  = 0;                               /* initialization */
  data->part[1]  = 0;
  data->num_j[0] = 0;
  data->num_j[1] = 0;

  if( (data_ptr = fopen(data->in_data,"r")) == NULL) { /* open in_data file */
    printf("\n\nError in function basic_input_data()");
    printf("\nWrong file = %s for the input data\n", data->in_data);
    exit(1);
  }
  for( ; ; )   { /* read data from input file */
    if(!read_text_string(data_ptr, text_string)) break;
    for(i = 0; i < number_of_data; i++)  {
      if(!(strcmp(text_string,identity[i]))) break;
    }
    switch (i)  {
    case 0:                /* The number of proton particle j-orbits */
      if(!read_int_number(data_ptr, 1, &data->num_j[0])) break;
      control |= UL_ONE;

               /* memory for proton single-particle orbits in a JBAS structure. */

      data->jbas[0] = MALLOC(data->num_j[0] + 1, JBAS, func, "data.jbas[0]");
      break;
    case 1:                /* Orbit_Z */
      if(!data->num_j[0])   {
	printf("\n\nError in function basic_input_data():");
	printf("\n Number of proton single-particle orbit = ZERO\n");
	exit(1);
      }
      if(count_orbZ > data->num_j[0])   {
	printf("\n\nError in function basic_input_data():");
	printf("\nNumber of proton single-particle orbits = %d", data->num_j[0]);
	printf("\nTry to read orbit number %d\n",count_orbZ);
	exit(1);
      }
      if(!read_orbit_number(data_ptr, count_orbZ++, data->jbas[0])) break;
      control |= (UL_ONE << 1);
      break;
    case 2:                /* The number of neutron particle j-orbits */
      if(!read_int_number(data_ptr, 1, &data->num_j[1])) break;
      control |= (UL_ONE << 2);
      
               /* memory for neutron single-particle orbits in a JBAS structure. */

      data->jbas[1] = MALLOC(data->num_j[1] + 1, JBAS, func, "data.jbas[1]");
      break;
    case 3:                /* Orbit_N */
      if(!data->num_j[1])   {
	printf("\n\nError in function basic_input_data():");
	printf("\n Number of neutron single-particle orbit = ZERO\n");
	exit(1);
      }
      if(count_orbN > data->num_j[1])   {
	printf("\n\nError in function basic_input_data():");
	printf("\nNumber of neutron single-particle orbits = %d", data->num_j[1]);
	printf("\nTry to read orbit number %d\n",count_orbN);
	exit(1);
      }
      if(!read_orbit_number(data_ptr, count_orbN++, data->jbas[1])) break;
      control |= (UL_ONE << 3);
      break;
    case 4:                /* The proton number */
      if(!read_int_number(data_ptr, 1, &data->part[0])) break;
      control |= (UL_ONE << 4);
      break;
    case 5:                /* The neutron number */
      if(!read_int_number(data_ptr, 1, &data->part[1])) break;
      control |= (UL_ONE << 5);
      break;
    case 6:                /* Total angular momentum J is (even, odd) */
      if(!read_data_string(data_ptr, data_string)) break;
      if(!strcmp(data_string, "even")) data->j_value = EVEN;
      else if(!strcmp(data_string, "odd")) data->j_value = ODD; 
      else  {
	printf("\n\nError in function basic_input_data()");
	printf("\nError in init j_type (even or odd) = %s\n", data_string);
	exit(1);
      }
      control |= (UL_ONE << 6);         
      break;
    case 7:                /* Twice total projection of angular momentum */
      if(!read_int_number(data_ptr, 1, &data->MJ)) break;
      
             /*  test relation between particle number and total m-value */
 
      if(MOD(data->part[0] + data->part[1],2) != MOD(data->MJ,2))  {
	printf("\n\nError in function read_input_data()");
	printf("\n\n Wrong value of total MJ value !!!!!!");
	printf("\nProton number = %d  Neutron number = %d  Total MJ = %d\n\n",
	                              data->part[0], data->part[1], data->MJ);
	exit(1);
      }
      control |= (UL_ONE << 7);
      break;
    case 8:                /* Total parity (+, -) */
      if(!read_data_string(data_ptr, data_string)) break;
      if((data_string[0] != '+') && (data_string[0] !='-'))   {
	printf("\n\nError in function basic_input_data():");
	printf("\nThe total parity is not + or -\n");
	exit(1);
      }
      data->P = (data_string[0] == '+') ? +1 : -1;;
      control |= (UL_ONE << 8);
      break;
    case 9:                 /* Input proton-proton v_effective in J-scheme */
      if(!read_data_string(data_ptr, data-> veff_pp)) break;
      control |= (UL_ONE << 9);
      break;
    case 10:                /* Input neutron-neutron v_effective in J-scheme" */
      if(!read_data_string(data_ptr, data-> veff_nn)) break;
      control |= (UL_ONE << 10);
      break;
    case 11:                /* Input proton-neutron v_effective in J-scheme */
      if(!read_data_string(data_ptr, data-> veff_pn)) break;
      control |= (UL_ONE << 11);
      break;
    case 12:                /* Maximum dimension of the energy matrix */
      if(!read_int_number(data_ptr, 1, &data->max_iterations)) break;
      control |= (UL_ONE << 12);
      break;
    case 13:                /* Wanted number of converged eigenstates */
      if(!read_int_number(data_ptr, 1, &data->states)) break;
      control |= (UL_ONE << 13);
      break;
    case 14:                 /* Type of calculation process */
      if(!read_data_string(data_ptr, data->type_calc)) break;
      control |= (UL_ONE << 14);
      break;
    case 15:                 /* Total 2J value */
      if(!read_int_number(data_ptr, 1, &(data->tot_2J))) break;
      control |= (UL_ONE  << 15);
      break;
    case 16:                /* Init vector file */
      if(!read_data_string(data_ptr, data->start_vec)) break;
      control |= (UL_ONE << 16);
      break;
    case 17:                /* Number of start vectors */
      if(!read_int_number(data_ptr, 1, &data->num_start_vec)) break;
      control |= (UL_ONE  << 17);
      break;
    case 18:                /* List of vector no */
      if(   !strcmp(data->type_calc,"random-start") 
	 || !strcmp(data->type_calc,"random-continue")) break;
      data->vec_no = MALLOC(data->num_start_vec, int, func, "data.vec_no");
      if(!read_int_number(data_ptr, data->num_start_vec, data->vec_no)) break; 
      control |= (UL_ONE  << 18);
      break;
    case 19:                         /* Calculation title" */
      if(!read_data_string(data_ptr, data_string)) break;
      if(strlen(data_string) == 0) {
	strcpy(data->title,"lanczo");    /* no calculation title */
      }
      else   {
	strcpy(data->title, data_string);    /* no calculation title */
      }
      control |= (UL_ONE << 19);
      break;
    case 20:                            /* memory size for |SD> matrix elements */
      if(!read_int_number(data_ptr, 1, &number)) break;

                  /* transfer to max number of bytes in memory */
      
      if((1000 * number) > INT_MAX) {
	printf("\n\nError in function basic_input_data():");
        printf("\nMemory to store nondiag <SD(Z)SD(N)|OP|SD(Z)SD(N)>(in bytes)");
	printf("\n too large - must not exceed sizeof(int) = %d\n",INT_MAX);
        exit(1);
      }
      data->mem_nondiag  = (1000 * number) / sizeof(STORE);
      control |= (UL_ONE << 20);
    case 21:                            /* file size for |SD> matrix elements */
      if(!read_int_number(data_ptr, 1, &number)) break;

                  /* transfer to max number of bytes on file */
 
      if((1000000 * number) > INT_MAX) {
	printf("\n\nError in function basic_input_data():");
        printf("\nFile size to store nondiag <SD(Z)SD(N)|OP|SD(Z)SD(N)>(in bytes)");
	printf("\ntoo large - must not exceed sizeof(int) = %d\n", INT_MAX);
        exit(1);
      }
      data->file_nondiag = (1000000 * number) / sizeof(STORE);
      control |= (UL_ONE << 21);
    default:
      break;
    } /* end of switch */
  } /* end  input data */

  fclose(data_ptr);                                 /* close input data file */

  if(count_orbZ < data->num_j[0])   {
    printf("\n\nError in function basic_input_data():");
    printf("\nNumber of proton single-particle orbits = %d", data->num_j[0]);
    printf("\nOnly %d orbits found in input data file\n",count_orbZ);
    exit(1);
  }
        /* Re-order proton j_obits in data.j_bas[] after decreasing j-values.*/

  qsort(data->jbas[0],(UL) data->num_j[0],sizeof(JBAS),
                  (int(*)(const void *, const void *))j_comp);
  if(count_orbN < data->num_j[1])   {
    printf("\n\nError in function basic_input_data():");
    printf("\nNumber of neutron single-particle orbits = %d", data->num_j[1]);
    printf("\nOnly %d orbits found in input data file\n",count_orbN);
    exit(1);
  }
      /* Re-order neutron j_obits in data.j_bas[] after decreasing j-values.*/

  qsort(data->jbas[1],(UL) data->num_j[1],sizeof(JBAS),
	        (int(*)(const void *, const void *))j_comp);

  data->control = control;       /* store identification of data read */

} /* End: function input_data() */

    /*
    ** The function                            
    **               read_text_string()                   
    * reads a character string which identify some specific input data.
    ** The string starts following the symbol "\*" and ends with colon ":"
    ** or semicolon ";". If the string  contains the character \, EOF
    ** or reads text_string "END_OF_INPUT_FILE", the function returns FALSE,
    ** otherwise TRUE.
    */

static int read_text_string(FILE *in_data, char *text_string)
{
  char
         *ptr;
  int
         val;

  *text_string =  0x0;                      /* initialization */

  for( ; ; )   {       /* find next '\*' */
    while((val = getc(in_data)) && ((val != EOF) && (((char) val) != '\\')));
    if(val == EOF) return FALSE;
    if((val = getc(in_data)) && (((char) val) == '*')) break;
  }
  while((val = getc(in_data)) && isspace(val)) {  /* remove possible whitespace */
    if(val == EOF) return FALSE;
  }
  if(((char) val) == '\\')  {
    ungetc(val, in_data);     /* return it to the stream */
    return FALSE;
  }
  ptr = text_string;          /* point to start of text string */
  *(ptr++) = (char) val;
  while((val = getc(in_data)) && (val != EOF) && (((char) val) != '\\') 
                      && (((char) val) != ':') && (((char) val) != ';'))
                                                       *(ptr++) = (char) val;
  if(val == EOF)   return FALSE;
  if(((char) val) == '\\')  {
    ungetc(val, in_data);     /* return it to the stream */
    return FALSE;
  }
  *ptr = 0x0;    /* terminate the text_string */

  if(!strcmp(text_string,"END_OF_INPUT_FILE")) return FALSE;  

  return TRUE;    /* correct text_string found */

} /* End: function read_text_string() */

    /*
    ** The function                            
    **             read_int_number()                   
    ** reads (number) integers from FILE *in_data stream and
    ** store them in vector[]. If too few integers are read
    ** the function returns  FALSE, otherwise TRUE.
    */

static int read_int_number(FILE *in_data, int number, int *vector)
{
  char
          data_string[ONE_LINE], *ptr;
  int
          loop, phase, val;

  for(loop = 0; loop < number; loop++) vector[loop] = 0;       /* initialization */

  for(loop = 0; loop < number; loop++)  {

    if(!read_data_string(in_data, data_string))  return FALSE;

        /* convert data_string to an integer number */ 
 
    ptr = data_string;                 /* initialization */
    val = 0;
    phase = +1;
    do  {
      switch(*ptr) {
          case '+':  phase = +1;
	             break;
          case '-':  phase = -1;
	             break;
          default : if(!isdigit((int)(*ptr)))  return FALSE;
                    val *= 10;
                    val += (int) (*ptr - '0');
      }
    } while(*(++ptr));
    vector[loop] = phase * val;
  }  /* end loop through all integer data */
  if(loop < number)  return FALSE;
  else               return TRUE;
}  /* End: function read_int_number()  */

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

       /* read five integer numbers */

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

         /* read single-particle energy */

  if(read_float_number(in_data, 1, &dbl_val) == FALSE) {
    printf("\n\nError in function read_float_number();");
    printf("\nWrong floating point data for orbit no %d\n", orb_num);
    return FALSE;
  }
  jbas[orb_num].e = dbl_val;   /* store the single-particle energy */

  return TRUE;

} /* End: function  read_orbit_number() */

    /*
    ** The function                            
    **             read_float_number()                   
    ** reads (number) float from FILE *in_data stream,
    ** convert them to double and store them in vector[].
    ** If too few doubles are read the function returns 
    ** FALSE, otherwise TRUE.
    */

static int read_float_number(FILE  *in_data, int number, double *vector)
{
  char
          data_string[ONE_LINE], *ptr;
  int
          loop;

  for(loop = 0; loop < number; loop++) vector[loop] = 0.0;       /* initialization */

  for(loop = 0; loop < number; loop++)  {

    if(!read_data_string(in_data, data_string))  return FALSE;

         /* convert data_string to an integer number */ 
  
    ptr = data_string;              /* initialization */
    while(*ptr) {
      if(!(   (*ptr == '+') || (*ptr == '-') 
           || (*ptr == '.') || isdigit(((int)(*ptr))))) return FALSE;
      ptr++;
    }
    vector[loop] = atof(data_string);     /* convert and save float number */

  }  /* end loop through all floating poin numbers */

  if(loop < number)   return FALSE;
  else                 return TRUE;

} /* End: function read_float_number() */

    /*
     * The function                            
     *           read_data_string()                   
     * reads one piece of input data as a text_string containing no blanks. 
     * Whitespaces in front of the data_string is removed. Data_string must
     * not contain the character '\'. The function returns FALSE if no 
     * proper data_string is found, otherwise TRUE.
     */

static int read_data_string(FILE  *in_data, char *data_string)
{
  int
        val;

  *data_string =  0x0;                      /* initialization */

  while((val = getc(in_data)) && isspace(val)); /* remove whitespace before data string */

  if((val == EOF) || (((char) val) == '<')) return FALSE;

  if( (char) val == '\\')  {   /* data_string must not contain the character '\' */
    ungetc(val, in_data);     /* return it to the stream */
    return FALSE;
  }
  *(data_string++) = (char) val;     /* first character in data_string */

           /* read remaining data string */

  while((val = getc(in_data)) && (!isspace(val) && (val != EOF))
                               && ((char) val != '\\') && ((char) val != '<')) 
                                               *(data_string++) = (char) val;

  if( (char) val == '\\')  {   /* data_string must not contain the character '\' */
    ungetc(val, in_data);     /* return it to the stream */
    return FALSE;
  }
  *data_string = 0x0;         /* terminate the string */

  return TRUE;

}  /* End: function read_data_string() */

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

} /* End: function j_comp() */

   /*
   ** The function 
   **        particle_structure()
   ** calculates and return type of shell model calculation
   ** with the following coding 
   **      type =  0 - protons only
   **              1 - neutrons only
   **              2 - one proton and and one neutron
   **              3 - one proton and several neutrons
   **              4 - several protons and one neutron
   **              5 - several protons and several neutrons
   */

static int particle_structure(int numZ, int numN)
{
  if((numZ != 0) && (numN == 0))    return 0;
  if((numZ == 0) && (numN != 0))    return 1;
  if((numZ == 1) && (numN == 1))    return 2;
  if((numZ == 1) && (numN > 1))     return 3;
  if((numZ > 1) && (numN == 1))     return 4;
  if((numZ > 1) && (numN > 1))      return 5;

             /* should not reach here */

  printf("\n\nError in function particle_structure():");
  printf("\nWrong proton - neutron numbers: ");
  printf(" proton number = %d  neutron number = %d \n",numZ, numN);
  exit(1);

} /* End: function particle_structure() */  

     /*
     ** The function
     **              write_input_data()
     ** write to output file all essential data in input data file
     */

static void write_input_data(INPUT *data)
{
  char
          text_string[ONE_LINE], *ptr;
  int
          k;
  FILE
          *in_data_ptr, *out_data_ptr;

  if( (in_data_ptr = fopen(data->in_data,"r")) == NULL) { /* open in_data file */
    printf("\n\nError in function write_input_data()");
    printf("\nWrong file = %s for the input data\n", data->in_data);
    exit(1);
  }
  rewind(in_data_ptr);
  if((out_data_ptr = fopen(data->out_data,"w")) == NULL) { /* Open output data file */
    printf("\n\nError in function write_input_data()");
    printf("\nWrong file = %s for the output data\n", data->out_data);
    exit(1);
  }
  while(!feof(in_data_ptr))  {             /* Copy input file to output file */
    fgets(text_string,100, in_data_ptr);                    /* read a line*/
    for(k = 0, ptr = text_string;   (*ptr != '*') && (*ptr != '\n')
                                   && ( k < ONE_LINE); k++, ptr++);
    if(*(++ptr) != '*') fputs(text_string, out_data_ptr);  /* print a line */
  } /* All data are copied */
   
  fprintf(out_data_ptr,"\n                            -----  Results  ----");

  fclose(in_data_ptr);                   /* close the input  data file */
  fclose(out_data_ptr);                  /*  close the output data file */

}  /* End: function write_input_data() */

     /*
     ** The function 
     **            check_data()
     ** analyzes all input data and checks that 
     ** necessary input data are submitted
     */

static void check_data(INPUT *data)
{
  int
          k;
  UL
          test, remove, pos;
  FILE
          *file_ptr;

  test = 0X3FFFFF;                          /* initialization of 22 bits */

  switch(data->part_type) { 
     case 0:                                    /* protons only */
       remove =  (UL_ONE << 2)  + (UL_ONE << 3) + (UL_ONE << 5) /* protons only */
	       + (UL_ONE << 10) + (UL_ONE << 11);
       break;
     case 1:                                     /* neutrons only */
       remove =  UL_ONE        + (UL_ONE << 1)  + (UL_ONE << 4)
	       + (UL_ONE << 9) + (UL_ONE << 11);
       break;
     case 2:                                    /* one proton/one neutrons */
       remove =  (UL_ONE << 9)  + (UL_ONE << 10); 
       break;
     case 3:                                    /* one proton/several neutrons*/
       remove =  (UL_ONE << 9); 
       break;
     case 4:                                    /* several protons/one neutron*/
       remove =  (UL_ONE << 10); 
       break;
     case 5:                                   /* several protons and neutrons */ 
       remove = 0;
       break;
     default:
       printf("\n\nError in function check_data()");
       printf("\nproton number = %d neutron number = %d are wrong",
	                              data->part[0], data->part[1]);
       printf("\nproducing particle_structure() = %d", data->part_type);
       printf("\nonly 0,....,5 is allowed\n");
       exit(1);
  } /* end switch() */

  if(   strcmp(data->type_calc,"random-start") 
     || strcmp(data->type_calc,"random-continue")) {
    remove +=  (UL_ONE << 15) + (UL_ONE << 16)
             + (UL_ONE << 17) + (UL_ONE << 18);
  }
  data->control |= remove;
  if(test != data->control) {

            /* missing input data */

    printf("\n\nError in function check_data()");
    printf("\nmissing input data - see output file for information\n");

    if((file_ptr = fopen(data->out_data,"w")) == NULL) { 
      printf("\n\nError in function check_data()");
      printf("\nWrong file = %s for the output data\n", data->out_data);
      exit(1);
    }
    fprintf(file_ptr,"\n\nError in function check_data():");
    fprintf(file_ptr,"\nThe following data are missing:");

    for(k = 0, pos = UL_ONE; k < number_of_data; k++, pos <<= 1) {
      if((test & pos) & !(data->control & pos)) {
            fprintf(file_ptr,"\n%s",identity[k]);
      }
    }
    fclose(file_ptr);
    exit(1);
  } /* test: data are missing */

} /* End: function check_data() */
