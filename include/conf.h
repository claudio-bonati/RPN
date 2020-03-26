#ifndef CONF_H
#define CONF_H

#include"macro.h"

#include<openssl/md5.h>
#include<stdio.h>

//#include"flavour_matrix.h"
#include"gparam.h"
#include"geometry.h"
#include"vec.h"


typedef struct Conf {
  long update_index;

  double **link;       // [volume] [STDIM]
  Vec *phi;         // [volume]

//  // for higgs field & co
//  GAUGE_VECS *higgs;    // [volume]
//  FMatrix *Qh;          // [volume]
//  double *Dh;   // [volume]
  } Conf;


// in conf_def.c
void init_conf(Conf *GC,
               GParam const * const param);
void read_conf(Conf *GC,
               GParam const * const param);
void free_conf(Conf *GC,
               GParam const * const param);
void write_conf_on_file_with_name(Conf const * const GC,
                                  GParam const * const param,
                                  char const * const namefile);
void write_conf_on_file(Conf const * const GC,
                        GParam const * const param);
void write_conf_on_file_back(Conf const * const GC,
                             GParam const * const param);
void compute_md5sum_conf(char *res,        // the lenght is 2*MD5_DIGEST_LENGTH
                         Conf const * const GC,
                         GParam const * const param);


// in gauge_conf_meas.c
double plaquettep(Conf const * const GC,
                  Geometry const * const geo,
                  GParam const * const param,
                  long r,
                  int i,
                  int j);
double plaquette(Conf const * const GC,
                 Geometry const * const geo,
                 GParam const * const param);
double polyakov(Conf const * const GC,
                Geometry const * const geo,
                GParam const * const param);

/*
void higgs_interaction(Gauge_Conf const * const GC,
                       Geometry const * const geo,
                       GParam const * const param,
                       double *he);
void compute_flavour_observables(Gauge_Conf const * const GC,
                                 GParam const * const param,
                                 double *tildeG0,
                                 double *tildeGminp,
                                 double *tildeD0,
                                 double *tildeDminp);
void compute_flavour_observables_corr(Gauge_Conf const * const GC,
                                      Geometry const * const geo,
                                      GParam const * const param,
                                      double *corrQQ,
                                      double *corr0string0,
                                      double *corr0string1);
void perform_measures_higgs(Gauge_Conf * GC,
                            Geometry const * const geo,
                            GParam const * const param,
                            FILE *datafilep);
void perform_measures_higgs_for_testing(Gauge_Conf *GC,
                                        Geometry const * const geo,
                                        GParam const * const param,
                                        FILE *datafilep);
*/

// in gauge_conf_upd.c
int metropolis_for_link(Conf *GC,
                        Geometry const * const geo,
                        GParam const * const param,
                        long r,
                        int i);
void calcstaples_for_phi(Conf *GC,
                         Geometry const * const geo,
                         GParam const * const param,
                         long r,
                         Vec *staple);
void overrelaxation_for_phi(Conf *GC,
                            Geometry const * const geo,
                            GParam const * const param,
                            long r);
int metropolis_for_phi(Conf *GC,
                       Geometry const * const geo,
                       GParam const * const param,
                       long r);
void update(Conf *GC,
            Geometry const * const geo,
            GParam const * const param);


#endif
