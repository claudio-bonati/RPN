#ifndef GAUGE_CONF_UPD_C
#define GAUGE_CONF_UPD_C

#include"../include/macro.h"

#include<math.h>
#ifdef OPENMP_MODE
  #include<omp.h>
#endif
#include<stdlib.h>

#include"../include/conf.h"
#include"../include/gparam.h"
#include"../include/random.h"


// perform an update with heatbath of the link variables
void heatbath_for_link(Conf *GC,
                      Geometry const * const geo,
                      GParam const * const param,
                      long r,
                      int i)
  {
  double v1v2=scal_prod_Vec(&(GC->phi[r]), &(GC->phi[nnp(geo, r, i)]));

  // energy = link^2 - 2 link * v1v2
  // probability = exp (-beta*energy) \propto \exp(-beta*(link-v1v2)^2 )

  GC->link[r][i] = gauss1()/sqrt(2*param->d_beta) + v1v2;
  }


// perform an update with overrelaxation of the link variables
void overrelaxation_for_link(Conf *GC,
                             Geometry const * const geo,
                             long r,
                             int i)
  {
  double v1v2=scal_prod_Vec(&(GC->phi[r]), &(GC->phi[nnp(geo, r, i)]));

  // energy = link^2 - 2 link * v1v2
  // probability = exp (-beta*energy) \propto \exp(-beta*(link-v1v2)^2 )

  GC->link[r][i] = 2*v1v2-GC->link[r][i];
  }


// compute the staple for the phi field
void calcstaples_for_phi(Conf *GC,
                         Geometry const * const geo,
                         long r,
                         Vec *staple)
  {
  int i;
  Vec v1, v2, tmp;

  zero_Vec(staple);

  for(i=0; i<STDIM; i++)
     {
     // forward
     equal_Vec(&v1, &(GC->phi[nnp(geo, r, i)]) );

     // backward
     equal_Vec(&v2, &(GC->phi[nnm(geo, r, i)]) );

     lin_comb_Vec(&tmp, GC->link[r][i], &v1, GC->link[nnm(geo, r, i)][i], &v2);

     plus_equal_Vec(staple, &tmp);
     }
  }


// perform an update of the higgs field with overrelaxation
void overrelaxation_for_phi(Conf *GC,
                            Geometry const * const geo,
                            long r)
  {
  double norm;
  Vec staple;

  calcstaples_for_phi(GC, geo, r, &staple);
  norm=norm_Vec(&staple);

  if(norm>MIN_VALUE)
    {
    double aux1, aux2;
    Vec v1;

    equal_Vec(&v1, &GC->phi[r]);

    aux1=scal_prod_Vec(&v1, &staple);
    aux2=aux1/(norm*norm);

    times_equal_real_Vec(&v1, 2.0);
    minus_equal_times_real_Vec(&v1, &staple, aux2);

    equal_Vec(&GC->phi[r], &v1);
    }
  else
    {
    rand_vec_Vec(&(GC->phi[r]) );
    }
  }


// perform an update of the phi field with metropolis
// retrn 0 if the trial state is rejected and 1 otherwise
int metropolis_for_phi(Conf *GC,
                       Geometry const * const geo,
                       GParam const * const param,
                       long r)
  {
  int acc=0;
  double old_energy, new_energy;
  Vec staple, new_vector;

  calcstaples_for_phi(GC, geo, r, &staple);

  old_energy=2.0*param->d_beta * scal_prod_Vec(&(GC->phi[r]), &staple);

  rand_rot_Vec(&new_vector, &(GC->phi[r]), param->d_epsilon_metro);

  new_energy=2.0*param->d_beta * scal_prod_Vec(&new_vector, &staple);

  if(casuale()< exp(old_energy-new_energy))
    {
    equal_Vec(&(GC->phi[r]), &new_vector);
    acc+=1;
    }

  return acc;
  }


// perform a complete update
void update(Conf * GC,
            Geometry const * const geo,
            GParam const * const param,
            double *acc)
   {
   int err, *a;
   long r, asum;
   int j, dir;

   err=posix_memalign((void**)&a, (size_t)INT_ALIGN, (size_t) param->d_volume * sizeof(int));
   if(err!=0)
     {
     fprintf(stderr, "Problems in allocating a vector! (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }

   for(r=0; r<param->d_volume; r++)
      {
      a[r]=0;
      }

   // heatbath on links
   for(dir=0; dir<STDIM; dir++)
      {
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<(param->d_volume)/2; r++)
         {
         heatbath_for_link(GC, geo, param, r, dir);
         }

      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=(param->d_volume)/2; r<(param->d_volume); r++)
         {
         heatbath_for_link(GC, geo, param, r, dir);
         }
      }

   // overrelax links and higgs
   for(j=0; j<param->d_overrelax; j++)
      {
      for(dir=0; dir<STDIM; dir++)
         {
         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(r)
         #endif
         for(r=0; r<(param->d_volume)/2; r++)
            {
            overrelaxation_for_link(GC, geo, r, dir);
            }

         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(r)
         #endif
         for(r=(param->d_volume)/2; r<(param->d_volume); r++)
            {
            overrelaxation_for_link(GC, geo, r, dir);
            }
         }

      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<(param->d_volume)/2; r++)
         {
         overrelaxation_for_phi(GC, geo, r);
         }

      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=(param->d_volume)/2; r<(param->d_volume); r++)
         {
         overrelaxation_for_phi(GC, geo, r);
         }
      }

   // metropolis on higgs
   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r)
   #endif
   for(r=0; r<(param->d_volume)/2; r++)
      {
      a[r]+=metropolis_for_phi(GC, geo, param, r);
      }

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r)
   #endif
   for(r=(param->d_volume)/2; r<(param->d_volume); r++)
      {
      a[r]+=metropolis_for_phi(GC, geo, param, r);
      }

   // acceptance computation
   asum=0;
   #ifdef OPENMP_MODE
   #pragma omp parallel for reduction(+:asum) private(r)
   #endif
   for(r=0; r<param->d_volume; r++)
      {
      asum+=(long)a[r];
      }

   *acc=((double)asum)*param->d_inv_vol;

   // final unitarization
   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r)
   #endif
   for(r=0; r<(param->d_volume); r++)
      {
      unitarize_Vec(&(GC->phi[r]));
      }

   free(a);

   GC->update_index++;
   }



#endif
