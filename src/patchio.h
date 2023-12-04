/* header file for patchio.a which handles file-io for patch clamp data       */

#ifndef patchio_h
#define patchio_h

#include <stdio.h>
#include <iostream>	//.h weg
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>	//.h weg
#include <errno.h>
//#include <d_nrutil.h>
//#include <d_nr.h>
#include "daten.h"


#define WRITE 1
#define READ 0
#define OPEN 1
#define DONT_OPEN 0

class Tchannel_type
{
   public:
   double **k, **k_T, **k_start, *level, *symbol, *symbol_q, *q, *initial_state_value;
   int no_states,
       no_symbols,
       no_fit,
       *symbol_index,
       **fit_index,
       **error_index;
   char *name;
   
   int nm;
   int nf;
   
   Tchannel_type(int No_states, int No_fit = 0, char *Name = (char *)"");
   ~Tchannel_type();
   void init();
   void init_steady_state();
   void init_symbols();
   Tchannel_type * operator + (Tchannel_type &channel2);
/* creates a macro channel out of the two single channel by applying the     */
/* Kronecker sum (see Albertsen and Hansen 1994)                             */
   Tchannel_type * operator * (int factor);
/* creates a macro channel by adding the same channel factor times           */

   double **kronecker_add_rate_constants(double **channel1, int N1,
					 double **channel2, int N2);

   Tchannel_type *reorder_states(Tchannel_type *channel);
   Tchannel_type *reduce_states(Tchannel_type *channel);
   Tchannel_type *merge_channels(Tchannel_type *channel1, Tchannel_type *channel2);
   
   friend Tchannel_type *copy(Tchannel_type *d,const Tchannel_type& r);
   
   friend std::ostream& operator << (std::ostream& str, const Tchannel_type &channel);
/* writes a channel in an ASCII stream                                       */
//   friend istream& operator >> (istream& str, Tchannel_type& channel);
/* reads a channel from an ASCII stream                                      */
};

unsigned short int * load_data(const char *datafile_name, 
                               unsigned long int &laenge);
/* The datafile_name is looked for in the actual directory or in
   the path, specified by the DATEN environment variable. If found, an array
   for the data is allocated, and the data is loaded and returned.           */

int save_data(const char *datafile_name,short int *time_series,
	      unsigned long int laenge);
/* The datafile_name is looked for in the actual directory or in
   the path, specified by the DATEN environment variable. If found, an array
   for the data is allocated, and the data is loaded and returned.           */

FILE* find_file(char *name, const char *envvar, int write, int open_file);
/* Used by "load_data" to look for files on disk                          */

inline  extern unsigned short int convert_short_int(short int value)
/* Used by" load_data" to sort the short int values of the data
   in the right order, independent of the platform (pc, workstation, etc)    */
{
  unsigned char *ByteA = (unsigned char*)&value;

  return ByteA[0] + (ByteA[1] << 8);
}

inline  extern float convert_float(float value)
{
  unsigned char *ByteA = (unsigned char*)&value;

  return ByteA[0] + (ByteA[1] << 8);
}

inline extern unsigned long int convert_long_int(long int value)
/* Used by "load_data" to sort the long int number of the length of the data
   in the right order, independent of the platform (pc, workstation, etc)    */
{
  unsigned short int *siA = (unsigned short int*)&value;

  return  convert_short_int(siA[0]) + (convert_short_int(siA[1]) << 16);
}



inline extern double convert_double(double value)
{
  float *siA = (float*)&value;
  double tmp;
  std::cout << "input " << std::hex << value << std::endl;
  tmp = (double) ((unsigned char) convert_float(siA[1]) << 16);
  return  convert_float(siA[0]) + tmp;
}



inline extern unsigned short int *inv_convert_si(short int value)
/* Used by" save_data" to sort the short int values of the data
   in the right order, independent of the platform (pc, workstation, etc)    */
{
  static unsigned short int return_value;
  static unsigned char *ByteA = (unsigned char*)&return_value;
  //cout << "short_int_value befor:" << hex << value << endl;
  ByteA[0] = value & 0x00FF;
  ByteA[1] = value >> 8;
  //cout << "short_int_value after inv_conversion:" << hex << return_value << endl;
  return &return_value;
}



inline extern unsigned long int *inv_convert_li(long int value)
/* Used by "save_data" to sort the long int number of the length of the data
   in the right order, independent of the platform (pc, workstation, etc)    */
{
  static unsigned long int return_value;
  //static unsigned short int *siA = (unsigned short int*)&return_value;
  static unsigned char *ByteA = (unsigned char*)&return_value;


  /*  siA[0] = value & 0x0000FFFF;
      siA[1] = value >> 16;*/

  //Durch die richtige PC-Kodierung ersetzt fw 5.11.99
  ByteA[0] = value & 0x000000FF;
  ByteA[1] = (value & 0x0000FF00)>>8;
  ByteA[2] = (value & 0x00FF0000)>>16; 
  ByteA[3] = (value & 0xFF000000)>>24;
#ifdef debug_patchio
  cout << "long_int_value befor: " << hex << value 
       << "\tand after inv_conversion: " << hex << return_value << endl;
#endif
  return &return_value;
}

inline extern float *inv_convert_float(float fvalue)
/* Used by "save_data" to sort the long int number of the length of the data
   in the right order, independent of the platform (pc, workstation, etc)    */
{

  //cout << "sizeof(float) : " << sizeof(float) << "sizeof(long int) : " << sizeof(long int) << "sizeof(double) : " << sizeof(double) <<endl;
#ifdef debug_patchio
  float *fdummy;  
  fdummy = (float *) inv_convert_li((unsigned long int) fvalue);
  cout << "inv_conv_float vorher : " << hex << (unsigned long int)fvalue << "\thinterher : "
       << hex << (unsigned long int )*fdummy << endl;
  cout << hex << convert_float(*fdummy) << endl;
#endif  
  return (float *) inv_convert_li((unsigned long int) fvalue);
}

//*****************************************************************************************************//
inline double *double_vector(int a, int b)
{
	double *vec;

	vec = (double *) malloc( (unsigned)(b-a+1) * sizeof(double) );

	return vec-a;
}

inline void free_double_vector(double *vec, int a, int b)
{
	//free( (char *)(vec + a) );
}

inline int *integer_vector(int a, int b)
{
	int *vec;
	vec = (int *)malloc((unsigned) (b-a+1)*sizeof(int));
	return vec-a;
}

inline void free_integer_vector(int *vec, int a, int b)
{
	//free( (char *)(vec + a) );
}

inline double **real_matrix(int a, int b, int c, int d)
{
	int i;
	double **matrix;

	matrix = (double **) malloc( (unsigned)(b - a + 1) * sizeof(double *) );
	matrix -= a;

	for (i = a; i <= b; i++)
	{
		matrix[i] = (double *) malloc( (unsigned)(d - c + 1)*sizeof(double) );

		matrix[i] -= c;
	}
	return matrix;
}

inline void free_real_matrix(double **matrix, int a, int b, int c, int d)
{
	int i;

	//for (i = b; i >= a; i--)
		//free( (char *)(matrix[i] + c) );
	//free( (char *)(matrix + a) );
}

inline int **integer_matrix(long a,long b,long c,long d)
{
	long i,row,col; 
	row=b-a+1;
	col=d-c+1;
	int **matrix;

	matrix=(int **) malloc((unsigned int)((row+1)*sizeof(int*)));
	matrix += 1;
	matrix -= a;

	matrix[a]=(int *) malloc((unsigned int)((row*col+1)*sizeof(int)));

	matrix[a] += 1;
	matrix[a] -= c;

	for(i=a+1;i<=b;i++) matrix[i]=matrix[i-1]+col;

	return matrix;
}

inline void free_integer_matrix(int **matrix,long a,long b,long c,long d)
{
	//free((char*) (matrix[a]+c-1));
	//free((char*) (matrix+a-1));
}

#endif /* patchio_h */



