#ifndef VEC_H
#define VEC_H

#include<math.h>
#include<stdio.h>

//#include"flavour_matrix.h"
#include"macro.h"
#include"random.h"

typedef struct Vec {
   double comp[NCOLOR] __attribute__((aligned(DOUBLE_ALIGN)));
} Vec;


// A=1
inline void one_Vec(Vec * restrict A)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  int i;

  A->comp[0]=1.0;
  for(i=1; i<NCOLOR; i++)
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

  for(i=0; i<NCOLOR; i++)
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

  for(i=0; i<NCOLOR; i++)
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

  for(i=0; i<NCOLOR; i++)
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

  for(i=0; i<NCOLOR; i++)
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

  for(i=0; i<NCOLOR; i++)
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

  for(i=0; i<NCOLOR; i++)
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

  for(i=0; i<NCOLOR; i++)
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
  double aux, ris;

  ris=0.0;
  for(i=0; i<NCOLOR; i++)
     {
     aux=fabs(A->comp[i]);
     ris+=aux*aux;
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

  for(i=0; i<NCOLOR; i++)
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

  for(i=0; i<NCOLOR; i++)
     {
     ris+=(A->comp[i]) * B->comp[i];
     }

  return ris;
  }


// random rotation close to identity
inline void rand_rot_Vec(Vec * restrict A, Vec const * const restrict B, double epsilon)
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

  int i, j;
  double theta, tmp1, tmp2;

  equal_Vec(A, B);

  for(i=0; i<NCOLOR-1; i++)
     {
     for(j=i+1; j<NCOLOR; j++)
        {
        tmp1=A->comp[i];
        tmp2=A->comp[j];

        theta=(2.0*casuale()-1.0)*epsilon;

        A->comp[i]=  cos(theta)*tmp1 +sin(theta)*tmp2;
        A->comp[j]= -sin(theta)*tmp1 +cos(theta)*tmp2;
        }
     }
  }

/*
// initialize the flavour matrix with a vector
// FM[mf(i,j)]=\sum_{on_gauge}conj(v1[i])v1[j] - delta^{ij}/N
// i, j are the flavour indices
inline void init_FMatrix_SoNVecs(FMatrix * restrict fmatrix, SoNVecs const * const restrict v1)
  {
  #ifdef __INTEL_COMPILER
  __assume_aligned(&(fmatrix->comp), DOUBLE_ALIGN);
  __assume_aligned(&(v1->comp), DOUBLE_ALIGN);
  #endif

  int i, j, k;

  zero_FMatrix(fmatrix);

  for(i=0; i<NHIGGS; i++)
     {
     for(j=0; j<NHIGGS; j++)
        {
        for(k=0; k<NCOLOR; k++)
           {
           fmatrix->comp[mf(i,j)]+=( (v1->comp[NCOLOR*i+k])*(v1->comp[NCOLOR*j+k]) + 0.0*I );
           }
        }
     }

  for(i=0; i<NHIGGS; i++)
     {
     fmatrix->comp[mf(i,i)]-=( 1.0/(double)NHIGGS + 0.0*I);
     }
  }
*/


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
