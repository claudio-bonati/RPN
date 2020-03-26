#ifndef CONF_MEAS_C
#define CONF_MEAS_C

#include"../include/macro.h"

#include<math.h>
#include<stdio.h>
#include<stdlib.h>

//#include"../include/flavour_matrix.h"
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


// compute the average Polyakov loop
double polyakov(Conf const * const GC,
                Geometry const * const geo,
                GParam const * const param)
   {
   long r0;
   double ris=0.0;

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r0) reduction(+ : ris)
   #endif
   for(r0=0; r0<param->d_volume; r0++)
      {
      long r=r0;
      int i;
      double aux=1.0;

      for(i=0; i<param->d_size[0]; i++)
         {
         aux *= GC->link[r][0];
         r=nnp(geo, r, 0);
         }

      ris+=aux;
      }

   ris *= param->d_inv_vol;

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

/*

// compute flavour related observables
//
// flavour matrices Qh and Dh HAVE TO BE INITIALIZED before calling this function
//
// tildeG0=ReTr[(\sum_x Q_x)(\sum_y Q_y)]/volume/NHIGGS
// tildeGminp=ReTr[(\sum_x Q_xe^{ipx})(\sum_y Q_ye^{-ipy)]/volume/NHIGGS
//
// tildeG0 is NHIGGS*susceptibility, tildeGminp is used to compute the 2nd momentum correlation function
//
// tildeD0=conj(\sum_x D_x) (\sum_y D_y) / volume
// tildeDminp=(\sum_x D_x e^{ipx}) conj(\sum_y D_y e^{ipy}) /volume
//
// tildeD0 is a U1 susceptibility, tildeDminp is used to compute the 2nd momentum correlation function
void compute_flavour_observables(Gauge_Conf const * const GC,
                                 GParam const * const param,
                                 double *tildeG0,
                                 double *tildeGminp,
                                 double *tildeD0,
                                 double *tildeDminp)
  {
  int coord[STDIM];
  long r;
  const double p = 2.0*PI/(double)param->d_size[1];
  double complex D, Dp;
  FMatrix Q, Qp, Qmp, tmp1, tmp2;

  // Q =sum_x Q_x
  // Qp=sum_x e^{ipx}Q_x
  // Qmp=sum_x e^{-ipx}Q_x
  //
  // D, Dp and are the analogous of Q and Qp for D

  D=0.0+0.0*I;
  Dp=0.0+0.0*I;

  zero_FMatrix(&Q);
  zero_FMatrix(&Qp);
  zero_FMatrix(&Qmp);
  for(r=0; r<(param->d_volume); r++)
     {
     equal_FMatrix(&tmp1, &(GC->Qh[r]));
     equal_FMatrix(&tmp2, &tmp1);

     plus_equal_FMatrix(&Q, &tmp1);
     D+=(GC->Dh[r]);

     si_to_cart(coord, r, param);

     times_equal_complex_FMatrix(&tmp1, cexp(I*((double)coord[1])*p));
     plus_equal_FMatrix(&Qp, &tmp1);
     Dp+=((GC->Dh[r]) * cexp(I*((double)coord[1])*p) );

     times_equal_complex_FMatrix(&tmp2, cexp(-I*((double)coord[1])*p));
     plus_equal_FMatrix(&Qmp, &tmp2);
     }

  equal_FMatrix(&tmp1, &Q);
  times_equal_FMatrix(&tmp1, &Q);

  *tildeG0=retr_FMatrix(&tmp1)*param->d_inv_vol;
  *tildeD0=creal(conj(D)*D)*param->d_inv_vol;

  equal_FMatrix(&tmp1, &Qp);
  times_equal_FMatrix(&tmp1, &Qmp);
  *tildeGminp=retr_FMatrix(&tmp1)*param->d_inv_vol;
  *tildeDminp=creal(Dp*conj(Dp))*param->d_inv_vol;
  }


// compute correlators of flavour observables
//
// flavour matrices Qh and Dh HAVE TO BE INITIALIZED before calling this function
//
// corrQQ is the correlato ReTr[Q_x Q_{x+d}]/N_higgs
// corr0string0 is the correlator \sum_f Re[hf^{dag} U_{x,1}U_{x+1,1}....Q_{x+d-1,1} hf], where hf is the f-th flavour
// corr0string1 is the correlator Re[h0^{dag} U_{x,1}U_{x+1,1}....U_{x+d-1,1} h1], where h1 is the second flavour
void compute_flavour_observables_corr(Gauge_Conf const * const GC,
                                      Geometry const * const geo,
                                      GParam const * const param,
                                      double *corrQQ,
                                      double *corr0string0,
                                      double *corr0string1)
  {
  int dist;
  long r;
  double accumulator1, accumulator2;

  for(dist=0; dist<param->d_size[1]; dist++)
     {
     accumulator1=0.0;

     #ifdef OPENMP_MODE
     #pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : accumulator1)
     #endif
     for(r=0; r<(param->d_volume); r++)
        {
        int i;
        long r1;
        FMatrix tmp1;

        equal_FMatrix(&tmp1, &(GC->Qh[r]));
        r1=r;
        for(i=0; i<dist; i++)
           {
           r1=nnp(geo, r1, 1);
           }
        times_equal_FMatrix(&tmp1, &(GC->Qh[r1]));
        accumulator1+=retr_FMatrix(&tmp1);
        }
     accumulator1*=param->d_inv_vol;
     corrQQ[dist]=accumulator1;

     accumulator1=0.0;
     accumulator2=0.0;

     #ifdef OPENMP_MODE
     #pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : accumulator1) reduction(+ : accumulator2)
     #endif
     for(r=0; r<(param->d_volume); r++)
        {
        int i;
        long r1;
        GAUGE_VECS phi1, phi2;
        GAUGE_GROUP U;

        equal_vecs(&phi1, &(GC->higgs[r]));
        r1=r;
        one(&U);
        for(i=0; i<dist; i++)
           {
           times_equal(&U, &(GC->lattice[r1][1]));
           r1=nnp(geo, r1, 1);
           }
        matrix_times_vector_all_vecs(&phi2, &U, &(GC->higgs[r1]));
        accumulator1+=re_scal_prod_vecs(&phi1, &phi2);
        #if NHIGGS >1
         accumulator2+=re_scal_prod_single_vecs(&phi1, &phi2, 0, 1);
        #else
         accumulator2+=0.0;
        #endif
        }
     accumulator1*=param->d_inv_vol;
     accumulator2*=param->d_inv_vol;

     corr0string0[dist]=accumulator1;
     corr0string1[dist]=accumulator2;
     }
  }
*/

void perform_measures(Conf *GC,
                      Geometry const * const geo,
                      GParam const * const param,
                      FILE *datafilep)
   {
//   double plaqs, plaqt, polyre, polyim, he, tildeG0, tildeGminp, tildeD0, tildeDminp;
//   long r;

   double plaq, poly, he;

   plaq = plaquette(GC, geo, param);
   poly = polyakov(GC, geo, param);
   he = higgs_interaction(GC, geo, param);

   fprintf(datafilep, "%.12g %.12g %.12g", plaq, poly, he);

/*
   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r)
   #endif
   for(r=0; r<(param->d_volume); r++)
      {
      init_FMatrix_vecs(&(GC->Qh[r]), &(GC->higgs[r]));
      GC->Dh[r] = HiggsU1Obs_vecs(&(GC->higgs[r]));
      }

   compute_flavour_observables(GC,
                               param,
                               &tildeG0,
                               &tildeGminp,
                               &tildeD0,
                               &tildeDminp);
*/

//   fprintf(datafilep, "%.12g %.12g ", tildeG0, tildeGminp);
//   fprintf(datafilep, "%.12g %.12g ", tildeD0, tildeDminp);

   fprintf(datafilep, "\n");

   fflush(datafilep);
   }



#endif












