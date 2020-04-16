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

// perform an update with metropolis of the link variables
void metropolis_for_link(Conf *GC,
                         Geometry const * const geo,
                         GParam const * const param,
                         long r,
                         int i)
  {
  double old_energy, new_energy;
  double old_link, new_link;

  const double v1v2=scal_prod_Vec(&(GC->phi[r]), &(GC->phi[nnp(geo, r, i)]));

  old_link = GC->link[r][i];
  old_energy = old_link*old_link/param->d_beta - 2.0 * old_link * v1v2;

  new_link = old_link + 0.5*(2.0*casuale()-1);
  new_energy = new_link*new_link/param->d_beta - 2.0 * new_link * v1v2;

  if(casuale()< exp(old_energy-new_energy))
    {
    GC->link[r][i] = new_link;
    }
  }


// perform an update with heatbath of the link variables
void heatbath_for_link(Conf *GC,
                       Geometry const * const geo,
                       GParam const * const param,
                       long r,
                       int i)
  {
  const double v1v2=scal_prod_Vec(&(GC->phi[r]), &(GC->phi[nnp(geo, r, i)]));

  // probability = exp (-link^2/beta+2*link*v1v2) \propto \exp(-(link-beta*v1v2)^2/beta )
  // gaussian with sigma=sqrt(\beta/2)
  // generic gaussian = sigma*normal + mu

  GC->link[r][i] = gauss1()*sqrt(param->d_beta/2.0) + v1v2*param->d_beta;
  }


// perform an update with overrelaxation of the link variables
void overrelaxation_for_link(Conf *GC,
                             Geometry const * const geo,
                             GParam const * const param,
                             long r,
                             int i)
  {
  const double v1v2=scal_prod_Vec(&(GC->phi[r]), &(GC->phi[nnp(geo, r, i)]));

  // probability = exp (-link^2/beta+2*link*v1v2) \propto \exp(-(link-beta*v1v2)^2/beta )

  GC->link[r][i] = 2.0*v1v2*param->d_beta - GC->link[r][i];
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


// perform an update of the higgs field with overrelaxation (ok also in the Z_2 link case)
void overrelaxation_for_phi(Conf *GC,
                            Geometry const * const geo,
                            long r)
  {
  double norm, aux1;
  Vec staple, newlink;

  calcstaples_for_phi(GC, geo, r, &staple);
  norm=norm_Vec(&staple);

  #ifdef DEBUG
  double prod_before=scal_prod_Vec(&(GC->phi[r]), &staple);
  #endif

  if(norm>MIN_VALUE)
    {
    equal_Vec(&newlink, &staple);
    times_equal_real_Vec(&newlink, 1./norm);
    aux1=scal_prod_Vec(&(GC->phi[r]), &newlink);
    times_equal_real_Vec(&newlink, 2.0*aux1);

    minus_equal_Vec(&newlink, &(GC->phi[r]));

    equal_Vec(&(GC->phi[r]), &newlink);
    }
  else
    {
    rand_vec_Vec(&(GC->phi[r]) );
    }

  #ifdef DEBUG
  double prod_after=scal_prod_Vec(&(GC->phi[r]), &staple);

  if(fabs(prod_before-prod_after)>MIN_VALUE)
    {
    fprintf(stderr, "Problem in overrelaxation (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  if(fabs(norm_Vec(&(GC->phi[r]))-1)>MIN_VALUE)
    {
    fprintf(stderr, "Problem in overrelaxation (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif
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

  old_energy=-2.0 * scal_prod_Vec(&(GC->phi[r]), &staple);

  rand_rot_Vec(&new_vector, &(GC->phi[r]), param->d_epsilon_metro);

  new_energy=-2.0 * scal_prod_Vec(&new_vector, &staple);

  if(casuale()< exp(old_energy-new_energy))
    {
    equal_Vec(&(GC->phi[r]), &new_vector);
    acc+=1;
    }

  #ifdef DEBUG
  if(fabs(norm_Vec(&(GC->phi[r]))-1)>MIN_VALUE)
    {
    fprintf(stderr, "Problem in metropolis (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

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
      for(r=0; r<(param->d_volume); r++)
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
         for(r=0; r<(param->d_volume); r++)
            {
            overrelaxation_for_link(GC, geo, param, r, dir);
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

   // metropolis on phi
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


// perform an update of the phi field with metropolis
// retrn 0 if the trial state is rejected and 1 otherwise
int metropolis_for_phi_without_links(Conf *GC,
                                     Geometry const * const geo,
                                     GParam const * const param,
                                     long r)
  {
  int i, acc=0;
  double old_energy, new_energy, tmp;
  Vec old_vector, new_vector;

  equal_Vec(&old_vector, &(GC->phi[r]));

  old_energy=0.0;
  for(i=0; i<STDIM; i++)
     {
     tmp = scal_prod_Vec(&old_vector, &(GC->phi[nnp(geo, r, i)]) );
     old_energy -= param->d_beta * tmp * tmp;

     tmp = scal_prod_Vec(&old_vector, &(GC->phi[nnm(geo, r, i)]) );
     old_energy -= param->d_beta * tmp * tmp;
     }

  rand_rot_Vec(&new_vector, &old_vector, param->d_epsilon_metro);

  if(casuale()<0.5)
    {
    times_equal_real_Vec(&new_vector, -1.0);
    }

  new_energy=0.0;
  for(i=0; i<STDIM; i++)
     {
     tmp = scal_prod_Vec(&new_vector, &(GC->phi[nnp(geo, r, i)]) );
     new_energy -= param->d_beta * tmp * tmp;

     tmp = scal_prod_Vec(&new_vector, &(GC->phi[nnm(geo, r, i)]) );
     new_energy -= param->d_beta * tmp * tmp;
     }

  if(casuale()< exp(old_energy-new_energy))
    {
    equal_Vec(&(GC->phi[r]), &new_vector);
    acc+=1;
    }

  return acc;
  }


// perform a complete update without using link variables
void update_without_links(Conf * GC,
                          Geometry const * const geo,
                          GParam const * const param,
                          double *acc)
   {
   int err, *a;
   long r, asum;

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

   // metropolis on phi
   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r)
   #endif
   for(r=0; r<(param->d_volume)/2; r++)
      {
      a[r]+=metropolis_for_phi_without_links(GC, geo, param, r);
      }

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r)
   #endif
   for(r=(param->d_volume)/2; r<(param->d_volume); r++)
      {
      a[r]+=metropolis_for_phi_without_links(GC, geo, param, r);
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


// perform an update with metropolis of the Z_2 link variables
void metropolis_for_link_z2(Conf *GC,
                            Geometry const * const geo,
                            GParam const * const param,
                            long r,
                            int i)
  {
  double old_energy;
  double old_link, new_link;

  const double v1v2=scal_prod_Vec(&(GC->phi[r]), &(GC->phi[nnp(geo, r, i)]));

  old_link = GC->link[r][i];
  old_energy = - param->d_beta * old_link * v1v2;

  new_link = -old_link;
  // new_energy = - old_energy;

  if(casuale()< exp(2*old_energy))
    {
    GC->link[r][i] = new_link;
    }
  }


// perform an update of the phi field with metropolis using Z_2 links
// retrn 0 if the trial state is rejected and 1 otherwise
int metropolis_for_phi_z2(Conf *GC,
                          Geometry const * const geo,
                          GParam const * const param,
                          long r)
  {
  int acc=0;
  double old_energy, new_energy;
  Vec staple, new_vector;

  calcstaples_for_phi(GC, geo, r, &staple);

  old_energy=-param->d_beta * scal_prod_Vec(&(GC->phi[r]), &staple);

  rand_rot_Vec(&new_vector, &(GC->phi[r]), param->d_epsilon_metro);

  new_energy=-param->d_beta * scal_prod_Vec(&new_vector, &staple);

  if(casuale()< exp(old_energy-new_energy))
    {
    equal_Vec(&(GC->phi[r]), &new_vector);
    acc+=1;
    }

  #ifdef DEBUG
  if(fabs(norm_Vec(&(GC->phi[r]))-1)>MIN_VALUE)
    {
    fprintf(stderr, "Problem in metropolis (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  return acc;
  }


// perform a complete update with Z_2 links
void update_z2(Conf * GC,
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

   // metropolis on Z_2 links
   for(dir=0; dir<STDIM; dir++)
      {
      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<(param->d_volume); r++)
         {
          metropolis_for_link_z2(GC, geo, param, r, dir);
         }
      }

   // overrelax higgs
   for(j=0; j<param->d_overrelax; j++)
      {
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

   // metropolis on phi
   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r)
   #endif
   for(r=0; r<(param->d_volume)/2; r++)
      {
      a[r]+=metropolis_for_phi_z2(GC, geo, param, r);
      }

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r)
   #endif
   for(r=(param->d_volume)/2; r<(param->d_volume); r++)
      {
      a[r]+=metropolis_for_phi_z2(GC, geo, param, r);
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
