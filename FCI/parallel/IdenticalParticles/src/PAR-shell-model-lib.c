/*******************  The module PAR-shell-model-lib.c  ******************/

#include "PAR-shell-model.h"

         /* a macro used in function pythag() */

static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

           /**** local function declarations ****/

void  *mem_alloc(size_t size, char *func, char *var);
       /*
       ** reserves dynamic memory in  heap using malloc(). No initialization
       ** of the elements.   
       **         size  - number of bytes reserved       
       **         func  - name of calling function
       **         var   - name of variable  
       ** Returns a void pointer to the reserved memory location
       ** or if(size == 0) returns NULL pointer. 
       */

void  *mem_calloc(size_t num_elem, size_t size, char *func, char *var);
       /*
       ** reserves dynamic memory in  heap using calloc(). All elements are
       ** initialized to ZERO
       **         num_elem   - number of elements
       **         size       - number of bytes for each element       
       **         func       - name of calling function
       **         var        - name of variable  
       ** Returns a void pointer to the reserved memory location
       ** or if(num_elem = 0) returns a NULL pointer.
       */

  void  *mem_realloc(void *ptr, size_t num_elem, size_t size, char *func, char *var);
      /*
      ** reallocates previously reserved dynamic memory in heap using realloc().
      ** In case of increased memory no initialization of the new elements.   
      **         *ptr  - pointer to the previosly reserved memory 
      **         size  - number of bytes reserved       
      **         func  - name of calling function
      **         var   - name of variable  
      ** Returns a void pointer to the new reserved memory location.
      */

int read_text_string(FILE *in_data, char *text_string);
      /*
      ** reads a character string which identify a specific line of input data.
      ** The string starts following the symbol "\*" and ends with colon ":"
      ** or semicolon ";". If the string  contains the character \, EOF
      ** or reads text_string "END_OF_INPUT_FILE", the function returns FALSE,
      ** otherwise TRUE.
      */

int read_int_number(FILE *in_data, int number, int *vector);
      /*
      ** reads (number) integers from FILE *in_data stream and
      ** store them in vector[]. If too few integers are read
      ** the function returns  FALSE, otherwise TRUE.
      */

int read_float_number(FILE *in_data, int number, double *vector);
      /*
      ** reads (number) float from FILE *in_data stream,
      ** convert them to double and store them in vector[].
      ** If too few doubles are read the function returns 
      ** FALSE, otherwise TRUE.
      */

int read_data_string(FILE *in_data, char *model_string);
      /*
      ** reads one piece of input data as a text_string containing no blanks. 
      ** Whitespaces in front of the data_string is removed. Data_string must
      ** not contain the character '\'. The function returns FALSE if no 
      ** proper data_string is found, otherwise TRUE.
      */

double clebsch_gordan(int j1, int j2, int j3, int m1, int m2);
       /*
       ** returns the value of the Clebsch-Gordan coefficient                                  
       ** <j1/2, m1/2, j2/2, m2/2 | j3/2, (m1 + m2)/2> 
       */

static double fac_ratio(int m, int n);
       /*
       ** calculates and returns the ratio (n! / m!).           
       */

TID run_clock(double *reference_time);
       /*
       ** calculates the current wall clock time
       ** and return the time used from a reference point
       ** in format: sec, min, hour
       ** Input data:                            
       **    double   reference_time
       ** Return data 
       **    TID  run_time
       */

void hpcwall(double *retval);
       /*
       ** returns the current time of day 
       */

void run_time(int number, TID *run);
       /*
       ** calculates start and stop time and returns the difference.                
       ** Input data:                            
       **    int number      - = 1 for start     
       **                        2 for stop
       **    TID run  -  returns the time difference. 
       */

TID time_step(int clock_num, int num);
       /*
       ** calculates start and stop time and returns the difference.                
       ** Input data:                            
       **    int clock_num    may include up to 100 different clocks
       **    int type  = 0    total initialization
       **              = 1    for start - zero return values     
       **              = 2    for stop  - return time from last start
       **              = 3    return total cummulative time
       **    TID run        - returns the time difference. 
       */

TID cum_time(int clock_num, int type);
       /*
       ** calculates and sum up time used between e start and stop signal
       ** and returns the sum on request.                
       ** Input data:                            
       **    int clock_num    may include up to 100 different clocks
       **    int type  = 1    for start time     
       **              = 2    for stop  time from last start
       **    int       = 3    retuned total used time
       **    TID run        - returns ex_time. 
       */

TID cumClock(int clockNum, int type);
       /*
       ** calculates and sum up cpu time used between a start and
       ** stop signal and returns the summed CPU time on request.                
       ** Input data:                            
       **    int clock_num    may include up to 100 different clocks
       **    int type  = 0    Total initialization     
       **              = 1    initialization for each start time     
       **              = 2    for stop  time from last start
       **    int       = 3    retuned total used time
       **    TID run        - returns ex_time. 
       */

void add_a_vector(double factor, int dim, double *init_vec, double *final_vec);
     /*
     ** multiply |init_vec> by a constant, add the result to 
     ** a vector |final_vec> and store the result in final_vec>.   
     **        final_vec[] <--- final_vec[] + factor * init_vec[]              
     */

void  scale_a_vector(int dim, double factor, double *ptr);
     /*
     ** multiply each component of a vector by  factor.
     **        ptr[] <--- factor * ptr[]                 
     */

void append_basis_vector_to_file(char *file_namen, int n, int dim, double *vector);
     /*
     ** opens the file "lanc-store.dat" and store a vector with 
     **          int n             -   as local vector number on current file and
     ** If n > 0 the the wave functions are appended to the previous ones.
     ** If n = 0 a new file is created if "lanc-store.dat" does  not exists.
     ** If n = 0 and "lanc-store.dat" does exist the vector overwrites the previous ones.
     ** The vector is stored in vec[].
     ** Note that the first vector on file is numbered n = ZERO.
     */

void writeLanczosEigenVec(char *fileName, EIGEN_DATA eigenData, double *vector);
     /*
     ** opens the file "fileName" and store 
     **     EIGEN_DATA eigenData
     **     double     vector[] converted to float
     ** If n > 0 the vector[] is appended to the previous ones.
     ** If n = 0 a new file is created if "filenane" does  not exists.
     ** If n = 0 and "filename" does exist vector[] overwrites previous ones.
     ** Note that the first vector on file is numbered n = ZERO.
     */

void read_basis_vector_from_open_file(char *file_name, FILE *file_ptr,
                                               int n, int dim, double *vec_ptr, float *temp_mem);
     /*
     ** reads from open file sequencially  the Lanczos vectors in the format
     **          int n             -   as local vector number on current file and
     **          int dim           -   as number of components
     **          double vec_ptr[]  -   dim amplitudes 
     */

void eigenvalues(int dim, double *h_matrix, double *list_eigenvalue);
     /*
     ** diagonalyzes a general real symmetrix matrix of dimension dim 
     ** stored in h_matrix[] in the format specified in function 
     ** pn_eigenvalue_lanczo_process() in file lanc-pn-process.c
     ** The function returns all eigenvalues in increasing order
     ** in list_eigenvalue[].
     */

int eigenvalue_comp(const double *one, const double *two);
     /*
     ** uses by the library function qsort() to put
     ** eigenvalues in an increasing sequence.
     */

void eigenvectors(int dim, double *h_matrix, EIGEN *list_eigen);
     /*
     ** diagonalyzes a general real symmetrix matrix of dimension dim 
     ** stored in h_matrix[] in the format specified in function 
     ** pn_eigenvalue_lanczo_process() in file lanc-pn-process.c
     ** The function returns all eigenvalues and corresponding
     ** eigenvectors  in increasing order in structure 
     ** EIGEN list_eigen[]
     */

int eigenvector_comp(const EIGEN *one, const EIGEN *two);
     /*
     ** uses by the library function qsort() to put eigenvalues and 
     ** corresponding eigenvectors in an increasing sequence.
     */

void eigenvalue_tred2(double **a, int n, double d[], double e[]);
       /*
       ** (C) Copr. 1986-92 Numerical Recipes Software )%.
       */

void eigenvalue_tqli(double d[], double e[], int n);
       /*
       ** (C) Copr. 1986-92 Numerical Recipes Software )%.
       */

static double pythag(double a, double b);
       /*
       ** (C) Copr. 1986-92 Numerical Recipes Software )%.
       */

void eigenvector_tred2(double **a, int n, double d[], double e[]);
       /*
       ** (C) Copr. 1986-92 Numerical Recipes Software )%.
       */

void eigenvector_tqli(double d[], double e[], int n, double **z);
       /*
       ** (C) Copr. 1986-92 Numerical Recipes Software )%.
       */
               /**** End: local function declarations ****/

               /**** The function defin|<itions  ****/ 

       /*
       ** The function 
       **        void  *mem_alloc()
       ** reserves dynamic memory in  heap using malloc(). No initialization
       ** of the elements.   
       **         size  - number of bytes reserved       
       **         func  - name of calling function
       **         var   - name of variable  
       ** Returns a void pointer to the reserved memory location
       ** or if(size == 0) returns NULL pointer. 
       **/

void  *mem_alloc(size_t size, char *func, char *var)
{
  void  *new_memory;

    /* allow a call to library function malloc() during memory allocation */

#undef     malloc   

   if(size == 0) new_memory = NULL_PTR;
   else {
      if((new_memory = (void *) malloc(size)) == NULL)  {
         printf("\n\nError from function mem_alloc():");
         printf("\nmemory allocation from function %s: %s ", func, var); 
         printf("\nMemory request: %d bytes \n", (int)size);
       MPI_Abort(MPI_COMM_WORLD,Rank);
      }
   }

#define malloc             "don't call malloc directly anymore!"

   return new_memory;

} /* End: function  mem_alloc() */

       /*
       ** The function 
       **        void  *mem_calloc()
       ** reserves dynamic memory in  heap using calloc(). All elements are
       ** initialized to ZERO
       **         num_elem   - number of elements
       **         size       - number of bytes for each element       
       **         func       - name of calling function
       **         var        - name of variable  
       ** Returns a void pointer to the reserved memory location
       ** or if(num_elem = 0) returns a NULL pointer.
       */

void  *mem_calloc(size_t num_elem, size_t size, char *func, char *var)
{
   void  *new_memory;

     /* allow a call to library function calloc() during memory allocation */

#undef     calloc   

   if(num_elem == 0) new_memory = NULL_PTR;
   else {
      if((new_memory = (void *) calloc(num_elem, size)) == NULL)  {
         printf("\n\nError from function mem_calloc():");
         printf("\nmemory allocation from function %s: %s ", func, var);
         printf("\nMemory request: num_elem = %d size = %d bytes \n", 
                                           (int)num_elem, (int)size);
       MPI_Abort(MPI_COMM_WORLD,Rank);
      }
   }

#define calloc             "don't call calloc directly anymore!"

   return new_memory;

} /* End: function  mem_calloc() */

      /*
      ** The function 
      **        void  *mem_realloc()
      ** reallocates previously reserved dynamic memory in heap using realloc().
      ** In case of increased memory no initialization of the new elements.   
      **         *ptr  - pointer to the previosly reserved memory 
      **         size  - number of bytes reserved       
      **         func  - name of calling function
      **         var   - name of variable  
      ** Returns a void pointer to the new reserved memory location.
      **/

void  *mem_realloc(void *ptr, size_t num_elem, size_t size, char *func, char *var)
{
   void  *new_memory;

    /* allow a call to library function realloc() during memory allocation */

#undef     realloc   

   if((new_memory = (void *) realloc(ptr, num_elem * size)) == NULL)  {
      printf("\n\nError from function mem_realloc():");
      printf("\nmemory allocation from function %s: %s ", func, var); 
      printf("\nMemory request: %d bytes \n", (int)size);
    MPI_Abort(MPI_COMM_WORLD,Rank);
   }

#define realloc             "don't call realloc directly anymore!"

   return new_memory;

} /* End: function  mem_realloc() */

   /*
    ** The function                            
    **               read_text_string()                   
    * reads a character string which identify some specific input data.
    ** The string starts following the symbol "\*" and ends with colon ":"
    ** or semicolon ";". If the string  contains the character \, EOF
    ** or reads text_string "END_OF_INPUT_FILE", the function returns FALSE,
    ** otherwise TRUE.
    */

int read_text_string(FILE *in_data, char *text_string)
{
  char
         *ptr;
  int
         val;

  *text_string =  0x0;                      // initialization   

  for( ; ; )   {       // find next '\*'   
    while((val = getc(in_data)) && ((val != EOF) && (((char) val) != '\\')));
    if(val == EOF) return FALSE;
    if((val = getc(in_data)) && (((char) val) == '*')) break;
  }
  while((val = getc(in_data)) && isspace(val)) {  // remove possible whitespace   
    if(val == EOF) return FALSE;
  }
  if(((char) val) == '\\')  {
    ungetc(val, in_data);     // return it to the stream   
    return FALSE;
  }
  ptr = text_string;          // point to start of text string   
  *ptr = (char) val;
  ptr++;
  while((val = getc(in_data)) && (val != EOF) && (((char) val) != '\\') 
	&& (((char) val) != ':') && (((char) val) != ';')) {
    *ptr = (char) val;
    ptr++;
  }
  if(val == EOF)   return FALSE;
  if(((char) val) == '\\')  {
    ungetc(val, in_data);     // return it to the stream   
    return FALSE;
  }
  *ptr = 0x0;    // terminate the text_string   

  if(!strcmp(text_string,"END_OF_INPUT_FILE")) return FALSE;  

  return TRUE;    // correct text_string found   

} // End: function read_text_string()   

    /*
    ** The function                            
    **             read_int_number()                   
    ** reads (number) integers from FILE *in_data stream and
    ** store them in vector[]. If too few integers are read
    ** the function returns  FALSE, otherwise TRUE.
    */

int read_int_number(FILE *in_data, int number, int *vector)
{
  char
          data_string[ONE_LINE], *ptr;
  int
          loop, phase, val;

  for(loop = 0; loop < number; loop++) vector[loop] = 0;       // initialization   

  for(loop = 0; loop < number; loop++)  {

    if(!read_data_string(in_data, data_string)) return FALSE;

        // convert data_string to an integer number    
 
    ptr = data_string;                 // initialization   
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
      ptr++;
    } while(*ptr);
    vector[loop] = phase * val;
  }  // end loop through all integer data   
  if(loop < number)  return FALSE;
  else               return TRUE;
}  // End: function read_int_number()    
    /*
    ** The function                            
    **             read_float_number()                   
    ** reads (number) float from FILE *in_data stream,
    ** convert them to double and store them in vector[].
    ** If too few doubles are read the function returns 
    ** FALSE, otherwise TRUE.
    */

int read_float_number(FILE *in_data, int number, double *vector)
{
   char   data_string[ONE_LINE];
   int    loop;

   for(loop = 0; loop < number; loop++) vector[loop] = 0.0;       /* initialization */

   for(loop = 0; loop < number; loop++)  {

      if(!read_data_string(in_data, data_string))  return FALSE;

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

int read_data_string(FILE  *in_data, char *data_string)
{
  int
        val;

  *data_string =  0x0;                      // initialization   

  while((val = getc(in_data)) && isspace(val)); // remove whitespace before data string   

  if((val == EOF) || (((char) val) == '<')) return FALSE;

  if( (char) val == '\\')  {   // data_string must not contain the character '\'   
    ungetc(val, in_data);     // return it to the stream   
    return FALSE;
  }
  *data_string = (char) val;     // first character in data_string   
  data_string++;

           // read remaining data string   

  while((val = getc(in_data)) && (!isspace(val) && (val != EOF))
	&& ((char) val != '\\') && ((char) val != '<')) { 
    *data_string = (char) val;
    data_string++;
  }

  if( (char) val == '\\')  {   // data_string must not contain the character '\'   
    ungetc(val, in_data);     // return it to the stream   
    return FALSE;
  }
  *data_string = 0x0;         // terminate the string   

  return TRUE;

}  // End: function read_data_string()   

       /*
       ** The function                                 
       **           clebsch_gordan()                   
       ** returns the value of the Clebsch-Gordan coefficient                                  
       ** <j1/2, m1/2, j2/2, m2/2 | j3/2, (m1 + m2)/2> 
       */

double clebsch_gordan(int j1, int j2, int j3, int m1, int m2)
{
         /*   fac[n] = n! / sqrt{10^n}  */
   
  static double fac[150] = {
    .1000000000000000E+01,.3162277660168379E-01,.2000000000000000E-02,
    .1897366596101028E-03,.2400000000000000E-04,.3794733192202056E-05,
    .7200000000000002E-06,.1593787940724864E-06,.4032000000000001E-07,
    .1147527317321902E-07,.3628800000000002E-08,.1262280049054092E-08,
    .4790016000000003E-09,.1969156876524384E-09,.8717829120000006E-10,
    .4135229440701207E-10,.2092278988800002E-10,.1124782407870728E-10,
    .6402373705728006E-11,.3846755834917892E-11,.2432902008176642E-11,
    .1615637450665515E-11,.1124000727777609E-11,.8175125500367506E-12,
    .6204484017332402E-12,.4905075300220504E-12,.4032914611266062E-12,
    .3443362860754794E-12,.3048883446117143E-12,.2796010642932893E-12,
    .2652528598121915E-12,.2600289897927591E-12,.2631308369336940E-12,
    .2745906132211536E-12,.2952327990396046E-12,.3267628297331728E-12,
    .3719933267899019E-12,.4352480892045862E-12,.5230226174666021E-12,
    .6450376682011969E-12,.8159152832478993E-12,.1057861775849963E-11,
    .1405006117752883E-11,.1910498367185033E-11,.2658271574788455E-11,
    .3782786767026367E-11,.5502622159812102E-11,.8178384990311005E-11,
    .1241391559253610E-10,.1923556149721149E-10,.3041409320171346E-10,
    .4905068181788930E-10,.8065817517094409E-10,.1351836790901029E-09,
    .2308436973392420E-09,.4014955268976057E-09,.7109985878048654E-09,
    .1281573721857157E-08,.2350561331282885E-08,.4385545276195193E-08,
    .8320987112741415E-08,.1605109571087441E-07,.3146997326038803E-07,
    .6269557984667544E-07,.1268869321858846E-06,.2608136121621699E-06,
    .5443449390774448E-06,.1153317792981115E-05,.2480035542436839E-05,
    .5411367084667394E-05,.1197857166996993E-04,.2689449441079695E-04,
    .6123445837688631E-04,.1413574626231488E-03,.3307885441519399E-03,
    .7845339175584759E-03,.1885494701666058E-02,.4591092485552201E-02,
    .1132428117820634E-01,.2829031189597267E-01,.7156945704626409E-01,
    .1833212210859029E+00,.4753643337012861E+00,.1247684230710655E+01,
    .3314240134565367E+01,.8908465407274079E+01,.2422709538367284E+02,
    .6665313817722467E+02,.1854826422573992E+03,.5220273782040236E+03,
    .1485715964481768E+04,.4275404227490954E+04,.1243841405464136E+05,
    .3658035857041260E+05,.1087366156656748E+06,.3266626020337846E+06,
    .9916779348709544E+06,.3041882150138602E+07,.9426890448883294E+07,
    .2951234062064472E+08,.9332621544394462E+08,.2980746402685117E+09,
    .9614466715035176E+09,.3131572170660985E+10,.1029901674514568E+11,
    .3419676810361796E+11,.1146280563734714E+12,.3878597438312349E+12,
    .1324641819451836E+13,.4565884904381298E+13,.1588245541522752E+14,
    .5574945468249565E+14,.1974506857221085E+15,.7055650984616650E+15,
    .2543559733472202E+16,.9249958440832429E+16,.3393108684451918E+17,
    .1255404359589777E+18,.4684525849754318E+18,.1762838801735966E+19,
    .6689502913449167E+19,.2559641940120622E+20,.9875044200833661E+20,
    .3840998695345006E+21,.1506141741511150E+22,.5953547977784760E+22,
    .2372173242880062E+23,.9526867474051174E+23,.3856204823625829E+24,
    .1573076357315330E+25,.6466855489220516E+25,.2678949036508007E+26,
    .1118248651196012E+27,.4703162928493458E+27,.1992942746161532E+28,
    .8508021737644667E+28,.3659042881952573E+29,.1585214610157954E+30,
    .6917786472619537E+30,.3040758665204989E+31,.1346201247571762E+32,
    .6002457605114648E+32,.2695364137888182E+33,.1218859041294581E+34,
    .5550293832739345E+34,.2544977678223085E+35,.1174997204390920E+36,
    .5462031093002385E+36,.2556323917872885E+37,.1204487096628886E+38};

  double value;       /* The return value */
  int temp, change;   /* interchange parameters */
  int jz1, jz2, jz3,
     jt1, jt2, jt3,
     jt4, jt5, j4, j5; /* ang. mom. relations */

    /* Local variables */

  int     loop_min, loop_max, phase, loop;
  double                           factor;

  value = 0.0;              /* initialization */

       /*
       ** Interchange the angular momenta such that the smallest j-value 
       ** is placed in position 2. A parameter is set to remember the interchange 
       **    change  = -1 - j1 is placed in position 2.  
       **            =  0 - no change                       
       **            =  1 - j3 is placed in position 2.  
       **/

  if(   (j1 < j2 )  && (j1 < j3)) {

      /* interchange 1 and 2 */

    temp = j1; j1 = j2; j2 = temp;
    temp = m1; m1 = m2; m2 = temp;
    change = -1;
  }
  else if(   (j3 < j1) && (j3 < j2))  {

      /* interchange 2 and 3 */

    temp = j2; j2 = j3; j3 = temp;
    m2 = -(m1 + m2);
    change = 1;
  }
  else change = 0;

       /* Test of angular momentum relations  */

  jz1 = (j1 + j2 - j3)/2;
  jz2 = (j1 + j3 - j2)/2;
  jz3 = (j2 + j3 - j1)/2;

  if((jz1 < 0) || (jz2 < 0) || (jz3 < 0))  return  value ;

       /* Test of angular projection relations */

  if((j1 < abs(m1)) || (j2 < abs(m2)) || (j3 < abs(m1 + m2)))  return  value;

       /* Definition of loop parameters */

  jt1 = (j1 - j3 + m2)/2;
  jt2 = (j2 - j3 - m1)/2;
  jt3 = (j1 - m1)/2;
  jt4 = (j2 + m2)/2;
  jt5 = (j2 - m2)/2;
  j4  = j1/2;
  j5  = j3/2;

        /* Loop limits */

  loop_min = MAX( MAX(jt1, jt2) , 0);
  loop_max = MIN( MIN(jt3, jt4) , jz1);

          /* Loop test */

  if( loop_min > loop_max) return  value;

           /*  test for maximum factorials   */

  if ((jz1 > FAC_LIM) || (jz3 > FAC_LIM) || (jt4 > FAC_LIM) || (jt5 > FAC_LIM)) {
    printf("\n\nTermination from VCC. - too large factorials\n");
  MPI_Abort(MPI_COMM_WORLD,Rank);
  }
  for (loop = loop_min,phase = PHASE(loop_min); loop <= loop_max;
                              phase *= -1, loop++) {
    value +=  (fac_ratio(jt3-loop,j4) / fac [jt4-loop] )
          * (fac_ratio(loop-jt2,j5) / fac[loop-jt1] )
          * ((phase / fac[jz1-loop]) / fac[loop] );
  }
  factor =  fac_ratio(j4,(j1 + m1)/2) * fac_ratio(j4,jt3)
           * fac_ratio((j1 + j2 + j3)/2+1,jz2)
           * fac_ratio(j5,(j3 + m1 + m2)/2)
           * fac_ratio(j5,(j3 - m1 - m2)/2)
           * fac[jt4] * fac[jt5] * fac[jz1] * fac[jz3]
           * (j3 + 1);
  value *= sqrt(factor);

         /* Add the angular momentum interchange factor */

  if(change == -1) value *= PHASE(jz1);
  else if (change == 1)
    value *= PHASE(jt3) * sqrt(  (j2+1.0)
                      / (j3+1.0));
  return (value);
} /* End: function clebsch_gordan() */

       /*
       ** The function               
       **      double fac_ratio()    
       ** calculates and returns the ratio (n! / m!).           
       */

static double fac_ratio(int m, int n)
{
  int i, high;
  double value;

  value = 1.0;      /* initialization */
  if( n == m)    return (value);          /* Return for n equal m */
  high = MAX(n,m); 
  for (i = MIN(n,m)+1; i <= high; i++) value*= i;
  if (n < m) return (1.0/value);       /* Return for n less  m  */
  return (value);                  /* Return for n greater than m */
} /* End: function fac_ratio() */

    /*
    ** The function                           
    **      TID run_clock(..)                 
    ** calculates the current wall clock time
    ** and return the time used from a reference point
    ** in format: sec, min, hour
    ** Input data:                            
    **    double   reference_time
    ** Return data 
    **    TID  run_time
    */

TID run_clock(double *reference_time)
{
  unsigned long long int 
                          num_sec;
  double
                           retval;
  TID
                         run_time;

  hpcwall(&retval);

  num_sec       = (unsigned long long int) (retval - (*reference_time));
  run_time.sec  =  num_sec % 60;
  run_time.min  = num_sec / 60;
  run_time.hour = run_time.min/ 60;
  run_time.min  = run_time.min % 60;

  return run_time;

} /* End: of function run_clock(..) */

     /*
     ** The function 
     **       hpcwall()
     ** returns the current time of day 
     */

void hpcwall(double *retval)
{
  static long
                   zsec = 0, zusec = 0;
  struct timeval
                   tp;
  struct timezone
                   tpz;


  gettimeofday(&tp, &tpz);
  if(zsec == 0)   zsec  = tp.tv_sec;
  if(zusec == 0)  zusec = tp.tv_usec;

  *retval = (double)(tp.tv_sec - zsec) + (tp.tv_usec - zusec) * 0.000001;

} /* End: function hpcwall() */

    /*
    ** The function                           
    **      TID run_time(..)                 
    ** calculates start and stop time and returns the difference.                
    ** Input data:                            
    **    int number      - = 1 for start     
    **                        2 for stop
    **    TID run  -  returns the time difference. 
    */

void run_time(int number, TID *run)
{
  clock_t          num, stop = 0;
  static clock_t   start;

  if(number == 1) {
    start = clock();
    num       = (UL)0;
    run->sec  = (UL)0;
    run->min  = (UL)0;
    run->hour = (UL)0;
    run->min  = (UL)0;
  }
  else if(number == 2)  {
    stop = clock();
    num       = ((stop - start) /CLOCK_UNIT);
    run->sec  = (UL) (num % 60);
    run->min  = (UL) (num / 60);
    run->hour = (UL) (run->min/ 60);
    run->min  = (UL) (run->min % 60);
  }
  return;

} /* End: of function run_time(..) */

    /*
    ** The function                           
    **      TID wallClock(..)                 
    ** calculates and sum up wall time used between e start
    ** and stop signal and returns the sum on request.                
    ** Input data:                            
    **    int clock_num    may include up to 100 different clocks
    **    int type  = 0    Total initialization     
    **              = 1    initialization for each start time     
    **              = 2    for stop  time from last start
    **    int       = 3    retuned total ex_time used 
    **    TID run        - returns ex_time - for type = 3 only
    */

TID wallClock(int clock_num, int type)
{
  unsigned long long int   num_sec;
  static long              zsec[100], zusec[100];
  static double            retval[100];
  TID                      ex_time;
  struct timeval           tp;
  struct timezone          tpz;

  if(type == 0) {       // total initialization 
    retval[clock_num] = 0.0;
  }
  else if(type == 1) { 
    gettimeofday(&tp, &tpz);  // start time

    zsec[clock_num]  = tp.tv_sec;
    zusec[clock_num] = tp.tv_usec;
  }
  else if(type == 2) {
    gettimeofday(&tp, &tpz);

    retval[clock_num] +=  (double)(  tp.tv_sec - zsec[clock_num]) 
                        + (tp.tv_usec - zusec[clock_num]) * 0.000001;
  }
  else if(type == 3) {

    num_sec = (unsigned long long int)retval[clock_num];
    ex_time.sec  = num_sec % 60;
    ex_time.min  = num_sec / 60;
    ex_time.hour = ex_time.min/ 60;
    ex_time.min  = ex_time.min % 60;
  }
  else {
    printf("\n\nError in function cum_time(): ");
    printf("\nInput data type = %d is wrong !!\n\n", type);
  MPI_Abort(MPI_COMM_WORLD,Rank);
  }
  return ex_time;

} // End: function wallClock()

    /*
    ** The function                           
    **      TID cpuClock(..)                 
    ** calculates and sum up cpu time used between a start and
    ** stop signal and returns the summed CPU time on request.                
    ** Input data:                            
    **    int clock_num    may include up to 100 different clocks
    **    int type  = 0    Total initialization     
    **              = 1    initialization for each start time     
    **              = 2    for stop  time from last start
    **    int       = 3    retuned total used time
    **    TID run        - returns ex_time. 
    */

TID cpuClock(int clockNum, int type)
{
  UL               num_sec;
  static double    retval[100];
  static clock_t   startClock[100], stopClock[100];
  TID              ex_time;

  if(type == 0) {       // total initialization 
    retval[clockNum] = 0.0;
  }
  else if(type == 1) { 
    startClock[clockNum] = clock();
  }
  else if(type == 2) {
    stopClock[clockNum] = clock();
    retval[clockNum]
       += ((double)( stopClock[clockNum]
                    -startClock[clockNum])/CLOCK_UNIT); 
  }
  else if(type == 3) {
    num_sec      = (UL)retval[clockNum]; 
    ex_time.sec  = num_sec % 60;
    ex_time.min  = num_sec / 60;
    ex_time.hour = ex_time.min/ 60;
    ex_time.min  = ex_time.min % 60;
  }
  else {
    printf("\n\nError in function cpuClock(): ");
    printf("\nInput data type = %d is wrong !!\n\n", type);
  MPI_Abort(MPI_COMM_WORLD,Rank);
  }
  return ex_time;

} // End: function cpuClock()

    /*
    ** The function                           
    **      TID time_step(..)                 
    ** calculates start and stop wall time and returns
    **  the difference
    ** Input data:                            
    **    int clock_num    may include up to 100 different clocks
    **    int type  = 1    for start - no return values     
    **              = 2    for stop  - return time from last start
    **    TID run        - returns the time difference. 
    */

TID time_step(int clock_num, int type)
{
  unsigned long long int num_sec;
  static long            zsec[100], zusec[100];
  double                 retval;
  TID                    ex_time;
  struct timeval         tp;
  struct timezone        tpz;

  if(type == 1) { 
    gettimeofday(&tp, &tpz);  // start time
    zsec[clock_num]  = tp.tv_sec;
    zusec[clock_num] = tp.tv_usec;
  }
  else if(type == 2) {
    gettimeofday(&tp, &tpz);

    retval = (double)(  tp.tv_sec - zsec[clock_num]) 
                      + (tp.tv_usec - zusec[clock_num]) * 0.000001;

    num_sec = (unsigned long long int)retval;
    ex_time.sec  = num_sec % 60;
    ex_time.min  = num_sec / 60;
    ex_time.hour = ex_time.min/ 60;
    ex_time.min  = ex_time.min % 60;
  }
  else {
    printf("\n\nError in function time_step(): ");
    printf("\nInput data type = %d is wrong !!\n\n", type);
  MPI_Abort(MPI_COMM_WORLD,Rank);
  }
  return ex_time;

} // End: function time_step()

       /*
       ** The function                                               
       **                    add_a_vector()                          
       ** multiply |init_vec> by a constant, add the result to 
       ** a vector |final_vec> and store the result in final_vec>.   
       **        final_vec[] <--- final_vec[] + factor * init_vec[]              
       */

void add_a_vector(double factor, int dim, double *init_vec, double *final_vec)
{
  do  {
    *final_vec += factor * (*init_vec);
    init_vec++;
    final_vec++;
  } while (--dim);

} /* End: function add_a_vector() */

       /*
       ** The function                                   
       **             scale_a_vector()                   
       ** multiply each component of a vector by  factor.
       **        ptr[] <--- factor * ptr[]                 
       */

void  scale_a_vector( int dim, double factor, double *ptr)
{
  do  {
    *ptr *= factor;
    ptr++;
  } while (--dim);
} /* End: function scale_a_vector */

     /*
     ** The function                                           
     **            append_basis_vector_to_file( )                    
     ** opens the file "file_name" and store 
     **     int      n             - local vector number on current file
     **     int      dim           - dim of vector[]
     **     double   vector[]      - list of vector elementsa
     ** If n > 0 the the wave functions are appended to the previous ones.
     ** If n = 0 a new file is created if "lanc-store.dat" does  not exists.
     ** If n = 0 and "lanc-store.dat" does exist the vector overwrites
     ** the previous ones.
     ** Note that the first vector on file is numbered n = ZERO and the
     ** vector is transferred to type float and stored on file
     */

void append_basis_vector_to_file(char *file_name, int n, int dim, double *vector) 
{
  char    *func = {" append_basis_vector_to_file(): "};
  int     k;
  float   *temp_mem, *final_ptr;
  double  *init_ptr;
  FILE    *file_ptr;

  if((file_ptr = fopen(file_name,(n == 0) ? "wb": "ab+")) == NULL) {
    printf("\nError in function append_basis_vector_to_file();");
    printf("\nWrong file = %s to store a vector no %d\n",file_name, n);
  MPI_Abort(MPI_COMM_WORLD,Rank);
  }

  // vector number

  if(fwrite((const void *)&n,(size_t) sizeof(int), 1, file_ptr) != 1) {
    printf("\nError in function append_basis_vector_to_file()");
    printf("\nIn writing vector number =  %d", n);
    printf("\nto file %s\n",file_name);
  MPI_Abort(MPI_COMM_WORLD,Rank);
  }   /* end of if-test  */

  // dimension 

  if(fwrite((const void *)&dim,(size_t) sizeof(int),1,file_ptr)!=1) { 
    printf("\nError in function append_basis_vector_to_file()");
    printf("\nIn writing numSD =  %d", dim);
    printf("\nto file %s\n",file_name);
  MPI_Abort(MPI_COMM_WORLD,Rank);
  }  // end of if-test

  temp_mem = MALLOC(dim, float, func, "temp_mem");
  init_ptr = vector;
  final_ptr = temp_mem;
  for(k = 0; k < dim; k++)  {
    *final_ptr = (float)(*init_ptr);
    init_ptr++;
    final_ptr++;
  }
  if(fwrite((const void *)temp_mem,(size_t) sizeof(float), 
                              (size_t) dim, file_ptr) != dim) { 
    printf("\n\nError in function id_write_diag_SD_matr_elem():");
    printf("\nIn writing %d diag matrix elements to file %s\n\n",
                                                  dim, file_name);
  MPI_Abort(MPI_COMM_WORLD,Rank);
  }
  free(temp_mem);
 
  fclose(file_ptr);

} // End: function append_basis_vector_to_file()

     /*
     ** The function                                           
     **          writeLanczosEigenVec( )                    
     ** opens the file "fileName" and store 
     **     EIGEN_DATA eigenData
     **     double     vector[] converted to float
     ** If n > 0 the vector[] is appended to the previous ones.
     ** If n = 0 a new file is created if "filenane" does  not exists.
     ** If n = 0 and "filename" does exist vector[] overwrites previous ones.
     ** Note that the first vector on file is numbered n = ZERO.
     */

void writeLanczosEigenVec(char *fileName,EIGEN_DATA eigenData,double *vector) 
{
  char    *func = {"writeLanczosEigenVec(): "};
  int      n, k;
  float   *tempMem, *finalPtr;
  double  *initPtr;
  FILE     *filePtr;

  n = eigenData.vecNum;

  if((filePtr = fopen(fileName,(n == 0) ? "wb": "ab+")) == NULL) {
    printf("\nError in function writeLanczosEigenVec();");
    printf("\nWrong file = %s to store eigenVec no %d\n",fileName, n);
  MPI_Abort(MPI_COMM_WORLD,Rank);
  }

  if(fwrite((const void *)&eigenData,(size_t)sizeof(EIGEN_DATA),1,filePtr)!= 1) { 
    printf("\nError in function  writeLanczosEigenVec()");
    printf("\nIn writing vector number %d to file %s\n", n, fileName);
  MPI_Abort(MPI_COMM_WORLD,Rank);
  }

  tempMem = MALLOC(eigenData.dim, float, func, "tempMem");
  initPtr = vector;
  finalPtr = tempMem;
  for(k = 0; k < eigenData.dim; k++)  {
    *finalPtr = (float)(*initPtr);
    initPtr++;
    finalPtr++;
  }
  if(fwrite((const void *)tempMem,(size_t)sizeof(float),(size_t)eigenData.dim, 
                                                      filePtr) != eigenData.dim) { 
    printf("\n\nError in function writeLanczosEigenVec():");
    printf("\nIn writing eigenvector elements for vec num = %d\n", n);
  MPI_Abort(MPI_COMM_WORLD,Rank);
  }
  free(tempMem);   // local memory
  fclose(filePtr); // local file

} // End: function writeLanczosEigenVec()

     /*
     ** The function                                           
     **             readLanczosEigenVecFromOpenFile()
     ** reads from open file Lanczos eigenvectors
     ** in the format EIGEN_DAtA
     **   eigenData.vecNum   - local vector number and
     **   eigenData.dim      - number of components
     **   eigenData.eigenVal - energy eigenvalue
     **   eigenData.angMom   _ angular momentum 
     **   double vec_ptr[]   - eigenvector amplitudes 
     ** The vector is stored on file as type float - convert to double
     */

void readLanczosEigenVecFromOpenFile(char *file_name,FILE *file_ptr,int vecNum,
                      EIGEN_DATA *eigenData,double *vec_ptr,float *temp_mem)
{
  int    k;

  if( fread((void *)eigenData,(size_t)sizeof(EIGEN_DATA),1,file_ptr) != 1)  {
    printf("\nError in function  readLanczosEigenVecFromOpenFile()");
    printf("\nIn reading struct eigenDatafrom file %s\n",file_name);
  MPI_Abort(MPI_COMM_WORLD,Rank);
  } // end of if-test

   if(eigenData->vecNum != vecNum) {
     printf("\nError in function readLanczosEigenVecFromOpenFile()");
     printf("\nread vector number = %d is not the wanted n = %d\n",
                                         eigenData->vecNum, vecNum);
   MPI_Abort(MPI_COMM_WORLD,Rank);
   }

   if(fread((void *)temp_mem,(size_t) sizeof(float),   
       (size_t) eigenData->dim, file_ptr) != (size_t) eigenData->dim)  {
      printf("\nError in function readLanczosEigenVecFromOpenFile()");
      printf("\nIn reading vector number %d from file %s\n",
                                          eigenData->vecNum,file_name);
    MPI_Abort(MPI_COMM_WORLD,Rank);
   }
   for(k = 0; k < eigenData->dim; k++)  {
     *vec_ptr = (double)(*temp_mem);
     temp_mem++;
     vec_ptr++;
   }

} // End: function  readLanczosEigenVecFromOpenFile()

     /*
     ** The function                                           
     **             readShellmodelEigenVec()
     ** reads shell model eigenvectors number vecNum
     ** stored on file in format EIGEN_DAtA
     **   eigenData.vecNum   - local vector number and
     **   eigenData.dim      - number of components
     **   eigenData.eigenVal - energy eigenvalue
     **   eigenData.angMom   _ angular momentum 
     **   double vec_ptr[]   - eigenvector amplitudes 
     ** The eigenvectors  are stored on file as type float - convert to double
     */

void readShellmodelEigenVec(char *filename, int vecNum, int dim, 
           EIGEN_DATA *eigenData, double *vec_ptr,float *temp_mem)
{
  int       k;
  long     position;
  FILE     *filePtr;

 if((filePtr= fopen(filename,"rb")) == NULL)   {
    printf("\n\nError in function readShellmodelEigenVec()");
    printf("\nWrong file = %s to read the eigenvectors\n",filename);
  MPI_Abort(MPI_COMM_WORLD,Rank);
  }
  position = (long)(vecNum*(sizeof(EIGEN_DATA) + dim*sizeof(float)));
  fseek(filePtr, position,SEEK_SET);

  if(fread((void *)eigenData,(size_t)sizeof(EIGEN_DATA),1,filePtr) != 1)  {
    printf("\nError in function readShellmodelEigenVec()");
    printf("\nIn reading struct eigenDatafrom file %s\n",filename);
  MPI_Abort(MPI_COMM_WORLD,Rank);
  } // end of if-test

   if(  (eigenData->vecNum != vecNum) 
      ||(eigenData->dim != dim)) {
     printf("\nError in function readShellmodelEigenVec()");
     printf("\nread vector number %d, wanted %d\n",eigenData->vecNum,vecNum);
    printf("\nread dimension= %d   wanted %d\n",eigenData->dim, dim);
   MPI_Abort(MPI_COMM_WORLD,Rank);
   }
   if(fread((void *)temp_mem,(size_t) sizeof(float),   
       (size_t) eigenData->dim, filePtr) != (size_t) eigenData->dim)  {
      printf("\nError in function  readShellmodelEigenVec()");
      printf("\nIn reading vector number %d from file %s\n",
                                          eigenData->vecNum,filename);
    MPI_Abort(MPI_COMM_WORLD,Rank);
   }
   for(k = 0; k < eigenData->dim; k++)  {
     *vec_ptr = (double)(*temp_mem);
     temp_mem++;
     vec_ptr++;
   }

   fclose(filePtr);

} // End: function  readShellmodelEigenVec()

     /*
     ** The function                                           
     **            read_basis_vector_from_open_file( )                    
     ** reads from open file seqiencially  the Lanczos vectors
     ** in the format
     **   int n            - as local vector number and
     **   int dim          - as number of components
     **   double vec_ptr[] - dim amplitudes 
     ** The vector is stored on file as type float - convert to double
     */

void read_basis_vector_from_open_file(char *file_name, FILE *file_ptr,
                      int n, int dim, double *vec_ptr, float *temp_mem)
{
  int    k, vec_no, num_elem;
  float  *init_ptr;
  double *final_ptr;

  // vector number 

  if( fread((void *)&vec_no,(size_t) sizeof(int), 1, file_ptr) != 1)  {
    printf("\nError in function read_basis_vector_from_open__file()");
    printf("\nIn reading the vector number from file %s",file_name);
    printf("\nThe wanted basis vector number = %d\n",n);
  MPI_Abort(MPI_COMM_WORLD,Rank);
  } // end of if-test

   if(n != vec_no) {
     printf("\nError in function read_basis_vector_from_open__file()");
     printf("\nread vector number = %d is not the wanted n = %d\n", vec_no, n);
   MPI_Abort(MPI_COMM_WORLD,Rank);
   }

   // read dim 

   if( fread((void *)&num_elem,(size_t) sizeof(int),1,file_ptr)!=1) {
      printf("\nError in function read_basis_vector_from_open__file()");
      printf("\nIn reading dimension from  file %s\n",file_name);
    MPI_Abort(MPI_COMM_WORLD,Rank);
   } // end of if-test

   if(dim  != num_elem) {
     printf("\nError in function read_basis_vector_from_open__file()");
     printf("\nread vector dimension = %d is not the wanted dimension = %d\n", num_elem, dim);
   MPI_Abort(MPI_COMM_WORLD,Rank);
   }
   init_ptr = temp_mem;
   final_ptr = vec_ptr;

   if(fread((void *)temp_mem,(size_t) sizeof(float),   
       (size_t) dim, file_ptr) != (size_t) dim)  {
      printf("\nError in function read_basis_vector_from_open__file()");
      printf("\nIn reading vector number %d from file %s\n",n, file_name);
    MPI_Abort(MPI_COMM_WORLD,Rank);
   }
   for(k = 0; k < dim; k++)  {
     *final_ptr = (double)(*init_ptr);
     init_ptr++;
     final_ptr++;
   }

} // End: function read_basis_vector_from_open_file()

      /* 
      ** The function                                    
      **      eigenvalues()                          
      ** diagonalyzes a general real symmetrix matrix of dimension dim 
      ** stored in h_matrix[] as a vector with element sequence
      **      h(0,0), h(0,1), h(1,1), h(2,0),......
      ** The function returns all eigenvalues in increasing order
      ** in eigenvalue[].
      */

void eigenvalues(int dim, double *h_matrix, double *eigenvalue)
{
  char          
                    *func = {" eigenvalues(): "};
  register int
                    loop, left, right;
  double
                    *d, *e, **a;


    /* temporary memory allocation */

  d = MALLOC(dim + 1, double, func, "d[]");
  e = MALLOC(dim + 1, double, func, "e[]");
  a = MALLOC(dim + 1, double *, func,"a[][]");
  for(loop = 0; loop < dim + 1; loop++) {
    a[loop] = MALLOC(dim + 1, double, func,"a[]");
  }

   /* transfer the Lanczo energy matrix to a[][] */

  for (loop = 0, right = 1; right <= dim; right++)  {
    for (left = 1; left <= right; left++)   {
       a[left][right] = h_matrix[loop];
       a[right][left] = h_matrix[loop];
       loop++;
    }
  }
    /* diagonalize the matrix a[][] */

  eigenvalue_tred2(a, dim, d, e);       /* transfer to tri-diag. form */
  
  for( loop = dim; loop >= 0; loop--) free(a[loop]);
  free(a);

  eigenvalue_tqli(d, e, dim);   /* diagonalize */

    /* return eigenvalues */

  for(left = 1; left <= dim; left++) eigenvalue[left - 1] = d[left];

    /*
     * Sort eigenvalues and the corresponding eigenvectors 
     * after increasing values of the eigenvalues.         
     */

  qsort((void *) eigenvalue, (size_t) dim , (size_t) sizeof(double),
      (int(*)(const void *, const void *))eigenvalue_comp);

    /* release temporary memory */

  free(e);
  free(d);

} /* End: function eigenvalues() */

    /*
    ** The function                           
    **            eigenvalue_comp()                
    ** uses by the library function qsort() to put
    ** eigenvalues in an increasing sequence.
    */

int eigenvalue_comp(const double *one, const double *two)
{
  if(*one < *two)       return -1;
  else if(*one > *two)  return +1;
  else                  return  0;
} /* End: function eigenvalue_comp() */

      /* 
      ** The function                                    
      **      eigenvectors()                          
      ** diagonalyzes a general real symmetrix matrix of dimension dim 
      ** stored in h_matrix[] as a vector with element sequence
      **      h(0,0), h(0,1), h(1,1), h(2,0),......
      ** The function returns all eigenvalues and corresponding
      ** eigenvectors  in increasing order in structure 
      ** EIGEN eigen[]
      */

void eigenvectors(int dim, double *h_matrix, EIGEN *eigen)
{
  char
                    *func = {"eigenvectors(): "};
  register int
                    loop, left, right;
  double
                   *d, *e, **a;


    /* temporary memory allocation */

  
  d = MALLOC(dim + 1, double, func, "d[]");
  e = MALLOC(dim + 1, double, func, "e[]");
  a = MALLOC(dim + 1, double *, func,"a[][]");
  for(loop = 0; loop < dim + 1; loop++) {
    a[loop] = MALLOC(dim + 1, double, func,"a[]");
  }

   /* transfer the Lanczo energy matrix to a[][] */

  for (loop = 0, right = 1; right <= dim; right++)  {
    for (left = 1; left <= right; left++)   {
       a[left][right] = h_matrix[loop];
       a[right][left] = h_matrix[loop];
       loop++;
    }
  }
    /* diagonalize the matrix a[][] */

  eigenvector_tred2(a, dim, d, e);       /* transfer to tri-diag. form */
  eigenvector_tqli(d, e, dim, a);        /* diagonalize */

                /* return eigenvalues */

  for(left = 1; left <= dim; left++) eigen[left - 1].value = d[left];

    /*
     * Sort eigenvalues and the corresponding eigenvectors 
     * after increasing values of the eigenvalues.         
     */

                /* return eigenvactors */

  for(left = 1; left <= dim; left++) {

                /* transfer eigenvector no left */

    for(right = 1; right <= dim; right++) {
      eigen[left - 1].vector[right - 1] = a[right][left];
    }
  }
  qsort((void *) eigen, (size_t) dim , (size_t) sizeof(EIGEN),
      (int(*)(const void *, const void *))eigenvector_comp);

    /* release temporary memory */

  for( loop = dim; loop >= 0; loop--) free(a[loop]);
  free(a);
  free(e);
  free(d);

} /* End: function eigenvectors() */

    /*
    ** The function                           
    **            eigenvector_comp()                
    ** uses by the library function qsort() to put eigenvalues and 
    ** corresponding eigenvectors in an increasing sequence.
    */

int eigenvector_comp(const EIGEN *one, const EIGEN *two)
{
  if(one->value < two->value)       return -1;
  else if(one->value > two->value)  return +1;
  else                              return  0;
} /* End: function eigenvector_comp() */


/* (C) Copr. 1986-92 Numerical Recipes Software )%. */

void eigenvalue_tred2(double **a, int n, double d[], double e[])
{
  register int
                  l,k,j,i;
  double
                  scale,hh,h,g,f;

  for (i=n;i>=2;i--) {
    l=i-1;
    h=scale=0.0;
    if (l > 1) {
      for (k=1;k<=l;k++)
        scale += fabs(a[i][k]);
      if (scale == 0.0)
        e[i]=a[i][l];
      else {
        for (k=1;k<=l;k++) {
          a[i][k] /= scale;
          h += a[i][k]*a[i][k];
        }
        f=a[i][l];
        g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
        e[i]=scale*g;
        h -= f*g;
        a[i][l]=f-g;
        f=0.0;
        for (j=1;j<=l;j++) {
          a[j][i]=a[i][j]/h;
          g=0.0;
          for (k=1;k<=j;k++)
            g += a[j][k]*a[i][k];
          for (k=j+1;k<=l;k++)
            g += a[k][j]*a[i][k];
          e[j]=g/h;
          f += e[j]*a[i][j];
        }
        hh=f/(h+h);
        for (j=1;j<=l;j++) {
          f=a[i][j];
          e[j]=g=e[j]-hh*f;
          for (k=1;k<=j;k++)
            a[j][k] -= (f*e[k]+g*a[i][k]);
        }
      }
    } else
      e[i]=a[i][l];
    d[i]=h;
  }
  d[1]=0.0;
  e[1]=0.0;

       /* Contents of this loop can be omitted if eigenvectors not
          wanted except for statement d[i]=a[i][i]; */

  for (i=1;i<=n;i++) {
    d[i]=a[i][i];
  }
} /* End: function eigenvalue_tred2() */

void eigenvalue_tqli(double d[], double e[], int n)
{
  register int
                 m,l,iter,i;
  double
                 s,r,p,g,f,dd,c,b;

  for (i=2;i<=n;i++) e[i-1]=e[i];
  e[n]=0.0;
  for (l=1;l<=n;l++) {
    iter=0;
    do {
      for (m=l;m<=n-1;m++) {
        dd=fabs(d[m])+fabs(d[m+1]);
        if ((double)(fabs(e[m])+dd) == dd) break;
      }
      if (m != l) {
        if (iter++ == 30) {
          printf("\n\nToo many iterations in tqli.\n");
        MPI_Abort(MPI_COMM_WORLD,Rank);
        }
        g=(d[l+1]-d[l])/(2.0*e[l]);
        r=pythag(g,1.0);
        g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
        s=c=1.0;
        p=0.0;
        for (i=m-1;i>=l;i--) {
          f=s*e[i];
          b=c*e[i];
          e[i+1]=(r=pythag(f,g));
          if (r == 0.0) {
            d[i+1] -= p;
            e[m]=0.0;
            break;
          }
          s=f/r;
          c=g/r;
          g=d[i+1]-p;
          r=(d[i]-g)*s+2.0*c*b;
          d[i+1]=g+(p=s*r);
          g=c*r-b;
        }
        if (r == 0.0 && i >= l) continue;
        d[l] -= p;
        e[l]=g;
        e[m]=0.0;
      }
    } while (m != l);
  }
} /* End: function eigenvalue_tqli() */

static double pythag(double a, double b)
{
  double
           absa,absb;

  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
} /* End: function pythag() */


/* (C) Copr. 1986-92 Numerical Recipes Software )%. */


void eigenvector_tred2(double **a, int n, double d[], double e[])
{
  register int
                  l,k,j,i;
  double
                 scale,hh,h,g,f;

  for (i=n;i>=2;i--) {
    l=i-1;
    h=scale=0.0;
    if (l > 1) {
      for (k=1;k<=l;k++)
        scale += fabs(a[i][k]);
      if (scale == 0.0)
        e[i]=a[i][l];
      else {
        for (k=1;k<=l;k++) {
          a[i][k] /= scale;
          h += a[i][k]*a[i][k];
        }
        f=a[i][l];
        g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
        e[i]=scale*g;
        h -= f*g;
        a[i][l]=f-g;
        f=0.0;
        for (j=1;j<=l;j++) {
          a[j][i]=a[i][j]/h;
          g=0.0;
          for (k=1;k<=j;k++)
            g += a[j][k]*a[i][k];
          for (k=j+1;k<=l;k++)
            g += a[k][j]*a[i][k];
          e[j]=g/h;
          f += e[j]*a[i][j];
        }
        hh=f/(h+h);
        for (j=1;j<=l;j++) {
          f=a[i][j];
          e[j]=g=e[j]-hh*f;
          for (k=1;k<=j;k++)
            a[j][k] -= (f*e[k]+g*a[i][k]);
        }
      }
    } else
      e[i]=a[i][l];
    d[i]=h;
  }
  d[1]=0.0;
  e[1]=0.0;

        /* 
	** Contents of this loop can be omitted if eigenvectors
	** not wanted except for statement d[i]=a[i][i];
        */

  for (i=1;i<=n;i++) {
    l=i-1;
    if (d[i]) {
      for (j=1;j<=l;j++) {
        g=0.0;
        for (k=1;k<=l;k++)
          g += a[i][k]*a[k][j];
        for (k=1;k<=l;k++)
          a[k][j] -= g*a[k][i];
      }
    }
    d[i]=a[i][i];
    a[i][i]=1.0;
    for (j=1;j<=l;j++) a[j][i]=a[i][j]=0.0;
  }

} /* End: function eigenvector_tred2() */

void eigenvector_tqli(double d[], double e[], int n, double **z)
{
  register int
                 m,l,iter,i,k;
  double
                 s,r,p,g,f,dd,c,b;

  for (i=2;i<=n;i++) e[i-1]=e[i];
  e[n]=0.0;
  for (l=1;l<=n;l++) {
    iter=0;
    do {
      for (m=l;m<=n-1;m++) {
        dd=fabs(d[m])+fabs(d[m+1]);
        if ((double)(fabs(e[m])+dd) == dd) break;
      }
      if (m != l) {
        if (iter++ == 30) {
          printf("\n\nToo many iterations in tqli.\n");
        MPI_Abort(MPI_COMM_WORLD,Rank);
        }
        g=(d[l+1]-d[l])/(2.0*e[l]);
        r=pythag(g,1.0);
        g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
        s=c=1.0;
        p=0.0;
        for (i=m-1;i>=l;i--) {
          f=s*e[i];
          b=c*e[i];
          e[i+1]=(r=pythag(f,g));
          if (r == 0.0) {
            d[i+1] -= p;
            e[m]=0.0;
            break;
          }
          s=f/r;
          c=g/r;
          g=d[i+1]-p;
          r=(d[i]-g)*s+2.0*c*b;
          d[i+1]=g+(p=s*r);
          g=c*r-b;

  	     /*
	     ** Next loop can be obmitted if
<	     ** eigenvectors are not wanted
	     */
 
          for (k=1;k<=n;k++) {
            f=z[k][i+1];
            z[k][i+1]=s*z[k][i]+c*f;
            z[k][i]=c*z[k][i]-s*f;
          }
        }
        if (r == 0.0 && i >= l) continue;
        d[l] -= p;
        e[l]=g;
        e[m]=0.0;
      }
    } while (m != l);
  }
} /* End: function eigenvector_tqli() */
