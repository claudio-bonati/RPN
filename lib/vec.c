#ifndef VEC_C
#define VEC_C

#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#include"../include/endianness.h"
//#include"../include/flavour_matrix.h"
#include"../include/macro.h"
#include"../include/vec.h"


// A=1
void one_Vec(Vec * restrict A);


// A=0
void zero_Vec(Vec * restrict A);


// A=B
void equal_Vec(Vec * restrict A, Vec const * const restrict B);


// A+=B
void plus_equal_Vec(Vec * restrict A, Vec const * const restrict B);


// A-=B
void minus_equal_Vec(Vec * restrict A, Vec const * const restrict B);


// A-=(r*B)
void minus_equal_times_real_Vec(Vec * restrict A, Vec const * const restrict B, double r);


// A=b*B+c*C
void lin_comb_Vec(Vec * restrict A,
                  double b, Vec const * const restrict B,
                  double c, Vec const * const restrict C);


// A*=r
void times_equal_real_Vec(Vec * restrict A, double r);


// l2 norm
double norm_Vec(Vec const * const restrict A);


// unitarize
void unitarize_Vec(Vec * restrict A);


// random vector (normalized)
void rand_vec_Vec(Vec * restrict A);


// real part of the scalar product re(v_1^{\dag}v_2)
double scal_prod_Vec(Vec const * const restrict A, Vec const * const restrict B);


// random rotation close to identity
void rand_rot_Vec(Vec * restrict A, Vec const * const restrict B, double epsilon)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  int i, j, counter;
  double theta, tmp1, tmp2;

  equal_Vec(A, B);

  for(counter=0; counter<NFLAVOUR*NFLAVOUR/2; counter++)
     {
     i=(int)(casuale()*((double)NFLAVOUR - MIN_VALUE));
     j=(int)(casuale()*((double)NFLAVOUR - MIN_VALUE));
     if(i!=j)
       {
       tmp1=A->comp[i];
       tmp2=A->comp[j];

       theta=(2.0*casuale()-1.0)*epsilon*PI;

       A->comp[i]=  cos(theta)*tmp1 +sin(theta)*tmp2;
       A->comp[j]= -sin(theta)*tmp1 +cos(theta)*tmp2;
       }
     }
  }

// print on file
int print_on_file_Vec(FILE *fp, Vec const * const A)
  {
  int i, err;

  for(i=0; i<NFLAVOUR; i++)
     {
     err=fprintf(fp, "%.16f ", A->comp[i]);
     if(err<0)
       {
       fprintf(stderr, "Problem in writing on file a vector (%s, %d)\n", __FILE__, __LINE__);
       return 1;
       }
     }
  fprintf(fp, "\n");

  return 0;
  }


// print on binary file without changing endiannes
int print_on_binary_file_noswap_Vec(FILE *fp, Vec const * const A)
  {
  int i;
  size_t err;
  double aux;

  for(i=0; i<NFLAVOUR; i++)
     {
     aux=A->comp[i];
     err=fwrite(&aux, sizeof(double), 1, fp);
     if(err!=1)
       {
       fprintf(stderr, "Problem in binary writing on file a vector (%s, %d)\n", __FILE__, __LINE__);
       return 1;
       }
     }

  return 0;
  }


// print on binary file changing endiannes
int print_on_binary_file_swap_Vec(FILE *fp, Vec const * const A)
  {
  int i;
  size_t err;
  double aux;

  for(i=0; i<NFLAVOUR; i++)
     {
     aux=A->comp[i];

     SwapBytesDouble(&aux);

     err=fwrite(&aux, sizeof(double), 1, fp);
     if(err!=1)
       {
       fprintf(stderr, "Problem in binary writing on file a vector (%s, %d)\n", __FILE__, __LINE__);
       return 1;
       }
     }

  return 0;
  }


// print on binary file in bigendian
int print_on_binary_file_bigen_Vec(FILE *fp, Vec const * const A)
  {
  int err;

  if(endian()==0) // little endian machine
    {
    err=print_on_binary_file_swap_Vec(fp, A);
    }
  else
    {
    err=print_on_binary_file_noswap_Vec(fp, A);
    }

  return err;
  }


// read from file
int read_from_file_Vec(FILE *fp, Vec *A)
  {
  int i, err;
  double aux;

  for(i=0; i<NFLAVOUR; i++)
     {
     err=fscanf(fp, "%lg", &aux);
     if(err!=1)
       {
       fprintf(stderr, "Problems reading vector from file (%s, %d)\n", __FILE__, __LINE__);
       return 1;
       }
     A->comp[i]=aux;
     }

  return 0;
  }


// read from binary file without changing endiannes
int read_from_binary_file_noswap_Vec(FILE *fp, Vec *A)
  {
  size_t err;
  int i;
  double aux;

  for(i=0; i<NFLAVOUR; i++)
     {
     err=fread(&aux, sizeof(double), 1, fp);
     if(err!=1)
       {
       fprintf(stderr, "Problems reading ector from file (%s, %d)\n", __FILE__, __LINE__);
       return 1;
       }

     A->comp[i]=aux;
     }

  return 0;
  }


// read from binary file changing endianness
int read_from_binary_file_swap_Vec(FILE *fp, Vec *A)
  {
  int i;
  size_t err;
  double aux;

  for(i=0; i<NFLAVOUR; i++)
     {
     err=fread(&aux, sizeof(double), 1, fp);
     if(err!=1)
       {
       fprintf(stderr, "Problems reading vector from file (%s, %d)\n", __FILE__, __LINE__);
       return 1;
       }

     SwapBytesDouble(&aux);

     A->comp[i]=aux;
     }

  return 0;
  }


// read from binary file written in bigendian
int read_from_binary_file_bigen_Vec(FILE *fp, Vec *A)
  {
  int err;

  if(endian()==0) // little endian machine
    {
    err=read_from_binary_file_swap_Vec(fp, A);
    }
  else
    {
    err=read_from_binary_file_noswap_Vec(fp, A);
    }

  return err;
  }


#endif

