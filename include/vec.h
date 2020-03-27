#ifndef VEC_H
#define VEC_H

#include<math.h>
#include<stdio.h>

//#include"flavour_matrix.h"
#include"macro.h"
#include"random.h"

typedef struct Vec {
   double comp[NFLAVOUR] __attribute__((aligned(DOUBLE_ALIGN)));
} Vec;


// A=1
inline void one_Vec(Vec * restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;

  A->comp[0]=1.0;
  for(i=1; i<NFLAVOUR; i++)
     {
     A->comp[i]=0.0;
     }
  }


// A=0
inline void zero_Vec(Vec * restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NFLAVOUR; i++)
     {
     A->comp[i]=0.0;
     }
  }


// A=B
inline void equal_Vec(Vec * restrict A, Vec const * const restrict B)
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

  int i;

  for(i=0; i<NFLAVOUR; i++)
     {
     A->comp[i]=B->comp[i];
     }
  }


// A+=B
inline void plus_equal_Vec(Vec * restrict A, Vec const * const restrict B)
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

  int i;

  for(i=0; i<NFLAVOUR; i++)
     {
     A->comp[i]+=B->comp[i];
     }
  }


// A-=B
inline void minus_equal_Vec(Vec * restrict A, Vec const * const restrict B)
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

  int i;

  for(i=0; i<NFLAVOUR; i++)
     {
     A->comp[i]-=B->comp[i];
     }
  }


// A-=(r*B)
inline void minus_equal_times_real_Vec(Vec * restrict A, Vec const * const restrict B, double r)
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

  int i;

  for(i=0; i<NFLAVOUR; i++)
     {
     A->comp[i]-=(r*B->comp[i]);
     }
  }


// A=b*B+c*C
inline void lin_comb_Vec(Vec * restrict A,
                         double b, Vec const * const restrict B,
                         double c, Vec const * const restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  __assume_aligned(&(C->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NFLAVOUR; i++)
     {
     A->comp[i]=b*(B->comp[i])+c*(C->comp[i]);
     }
  }


// A*=r
inline void times_equal_real_Vec(Vec * restrict A, double r)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;

  for(i=0; i<NFLAVOUR; i++)
     {
     A->comp[i]*=r;
     }
  }


// l2 norm
inline double norm_Vec(Vec const * const restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;
  double ris;

  ris=0.0;
  for(i=0; i<NFLAVOUR; i++)
     {
     ris+=A->comp[i]*A->comp[i];
     }
  return sqrt(ris);
  }


// unitarize
inline void unitarize_Vec(Vec * restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  double norm;

  norm=norm_Vec(A);
  times_equal_real_Vec(A, 1./norm);
  }


// random vector (normalized)
inline void rand_vec_Vec(Vec * restrict A)
  {
  int i;

  for(i=0; i<NFLAVOUR; i++)
     {
     A->comp[i]=2.0*casuale()-1.0;
     }

  unitarize_Vec(A);
  }


// real part of the scalar product re(v_1^{\dag}v_2)
inline double scal_prod_Vec(Vec const * const restrict A, Vec const * const restrict B)
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

  int i;
  double ris=0.0;

  for(i=0; i<NFLAVOUR; i++)
     {
     ris+=(A->comp[i]) * B->comp[i];
     }

  return ris;
  }


// random rotation close to identity
void rand_rot_Vec(Vec * restrict A, Vec const * const restrict B, double epsilon);


// print on screen
void print_on_screen_Vec(Vec const * const A);


// print on file
int print_on_file_Vec(FILE *fp, Vec const * const A);


// print on binary file without changing endiannes
int print_on_binary_file_noswap_Vec(FILE *fp, Vec const * const A);


// print on binary file changing endiannes
int print_on_binary_file_swap_Vec(FILE *fp, Vec const * const A);


// print on binary file in bigendian
int print_on_binary_file_bigen_Vec(FILE *fp, Vec const * const A);


// read from file
int read_from_file_Vec(FILE *fp, Vec *A);


// read from binary file without changing endiannes
int read_from_binary_file_noswap_Vec(FILE *fp, Vec *A);


// read from binary file changing endianness
int read_from_binary_file_swap_Vec(FILE *fp, Vec *A);


// read from binary file written in bigendian
int read_from_binary_file_bigen_Vec(FILE *fp, Vec *A);



#endif // SON_H
