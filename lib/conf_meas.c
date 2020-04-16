#ifndef CONF_MEAS_C
#define CONF_MEAS_C

#include"../include/macro.h"

#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#include"../include/flavour_matrix.h"
#include"../include/gparam.h"
#include"../include/geometry.h"
#include"../include/conf.h"

// computation of the plaquette in position r and positive directions i,j
double plaquettep(Conf const * const GC,
                  Geometry const * const geo,
                  long r,
                  int i,
                  int j)
   {

//
//       ^ i
//       |   (2)
//       +---<---+
//       |       |
//   (3) V       ^ (1)
//       |       |
//       +--->---+---> j
//       r   (4)
//

   double ris;

   ris = GC->link[nnp(geo, r, j)][i];
   ris *= GC->link[nnp(geo, r, i)][j];
   ris *= GC->link[r][i];
   ris *= GC->link[r][j];

   return ris;
   }


// compute the average plaquettes
double plaquette(Conf const * const GC,
                 Geometry const * const geo,
                 GParam const * const param)
   {
   long r;
   double ris=0.0;

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : ris)
   #endif
   for(r=0; r<(param->d_volume); r++)
      {
      int i, j;
      i=0;
     
      for(i=0; i<STDIM; i++)
         {
         for(j=i+1; j<STDIM; j++)
            {
            ris+=plaquettep(GC, geo, r, i, j);
            }
         }
      }

   ris*=param->d_inv_vol;
   ris/=((double) (STDIM)*(STDIM-1)/2);

   return ris;
   }


// compute the average polyakov loop
double polyakov(Conf const * const GC,
                Geometry const * const geo,
                GParam const * const param)
   {
   long r;
   double ris=0.0;

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : ris)
   #endif
   for(r=0; r<(param->d_volume); r++)
      {
      int i;
      long raux=r;
      double poly=1.0;

      for(i=0; i<param->d_size[0]; i++)
         {
         poly*=GC->link[raux][0];
         raux=nnp(geo, raux, 0);
         }

      ris+=poly;
      }

   ris*=param->d_inv_vol;

   return ris;
   }


// compute the average value of phi_x U_{x,mu} phi_{x+mu}
double higgs_interaction(Conf const * const GC,
                         Geometry const * const geo,
                         GParam const * const param)
  {
  long r;
  double ris=0.0;

  #ifdef OPENMP_MODE
  #pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : ris)
  #endif
  for(r=0; r<(param->d_volume); r++)
     {
     int i;
     double aux=0.0;
     Vec v1;

     for(i=0; i<STDIM; i++)
        {
        equal_Vec(&v1, &(GC->phi[r]));

        aux+= GC->link[r][i] * scal_prod_Vec(&v1, &(GC->phi[nnp(geo, r, i)]));
        }

     ris+=aux;
     }

  ris/=(double) STDIM;
  ris*=param->d_inv_vol;

  return ris;
  }


// compute flavour related observables
//
// flavour matrices Qh HAVE TO BE INITIALIZED before calling this function
//
// tildeG0=Tr[(\sum_x Q_x)(\sum_y Q_y)]/volume
// tildeGminp=ReTr[(\sum_x Q_xe^{ipx})(\sum_y Q_ye^{-ipy)]/volume
//
// tildeG0 is the susceptibility, tildeGminp is used to compute the 2nd momentum correlation function
//
void compute_flavour_observables(Conf const * const GC,
                                 GParam const * const param,
                                 double *tildeG0,
                                 double *tildeGminp)
  {
  int coord[STDIM];
  long r;
  const double p = 2.0*PI/(double)param->d_size[1];
  FMatrix Q, Qp, Qmp, tmp1, tmp2;

  // Q =sum_x Q_x
  // Qp=sum_x e^{ipx}Q_x
  // Qmp=sum_x e^{-ipx}Q_x
  //

  zero_FMatrix(&Q);
  zero_FMatrix(&Qp);
  zero_FMatrix(&Qmp);
  for(r=0; r<(param->d_volume); r++)
     {
     equal_FMatrix(&tmp1, &(GC->Qh[r]));
     equal_FMatrix(&tmp2, &tmp1);

     plus_equal_FMatrix(&Q, &tmp1);

     si_to_cart(coord, r, param);

     times_equal_complex_FMatrix(&tmp1, cexp(I*((double)coord[1])*p));
     plus_equal_FMatrix(&Qp, &tmp1);

     times_equal_complex_FMatrix(&tmp2, cexp(-I*((double)coord[1])*p));
     plus_equal_FMatrix(&Qmp, &tmp2);
     }

  equal_FMatrix(&tmp1, &Q);
  times_equal_FMatrix(&tmp1, &Q);

  *tildeG0=retr_FMatrix(&tmp1)*param->d_inv_vol;

  equal_FMatrix(&tmp1, &Qp);
  times_equal_FMatrix(&tmp1, &Qmp);
  *tildeGminp=retr_FMatrix(&tmp1)*param->d_inv_vol;
  }


void perform_measures(Conf *GC,
                      GParam const * const param,
                      Geometry const * const geo,
                      FILE *datafilep)
   {
   long r;

   double tildeG0, tildeGminp;

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r)
   #endif
   for(r=0; r<(param->d_volume); r++)
      {
      init_FMatrix(&(GC->Qh[r]), &(GC->phi[r]));
      }

   if(param->d_beta<0)  // antiferromagnetic case
   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r)
   #endif
   for(r=0; r<(param->d_volume); r++)
      {
      times_equal_real_FMatrix(&(GC->Qh[r]), geo->d_parity[r]);
      }


   compute_flavour_observables(GC,
                               param,
                               &tildeG0,
                               &tildeGminp);

   fprintf(datafilep, "%.12g %.12g ", tildeG0, tildeGminp);

   fprintf(datafilep, "\n");

   fflush(datafilep);
   }


void perform_measures_z2(Conf *GC,
                         GParam const * const param,
                         Geometry const * const geo,
                         FILE *datafilep)
   {
   long r;

   double tildeG0, tildeGminp, plaq;

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r)
   #endif
   for(r=0; r<(param->d_volume); r++)
      {
      init_FMatrix(&(GC->Qh[r]), &(GC->phi[r]));
      }

   if(param->d_beta<0)  // antiferromagnetic case
   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r)
   #endif
   for(r=0; r<(param->d_volume); r++)
      {
      times_equal_real_FMatrix(&(GC->Qh[r]), geo->d_parity[r]);
      }


   compute_flavour_observables(GC,
                               param,
                               &tildeG0,
                               &tildeGminp);

   plaq=plaquette(GC, geo, param);

   fprintf(datafilep, "%.12g %.12g %.12g ", tildeG0, tildeGminp, plaq);

   fprintf(datafilep, "\n");

   fflush(datafilep);
   }


#endif












