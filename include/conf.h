#ifndef CONF_H
#define CONF_H

#include"macro.h"

#include<openssl/md5.h>
#include<stdio.h>

#include"flavour_matrix.h"
#include"gparam.h"
#include"geometry.h"
#include"vec.h"


typedef struct Conf {
  long update_index;

  double **link;    // [volume] [STDIM]
  Vec *phi;         // [volume]

  FMatrix *Qh;      // [volume]
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


// in conf_meas.c
double plaquettep(Conf const * const GC,
                  Geometry const * const geo,
                  long r,
                  int i,
                  int j);
double plaquette(Conf const * const GC,
                 Geometry const * const geo,
                 GParam const * const param);
double polyakov(Conf const * const GC,
                Geometry const * const geo,
                GParam const * const param);

double higgs_interaction(Conf const * const GC,
                         Geometry const * const geo,
                         GParam const * const param);
void compute_flavour_observables(Conf const * const GC,
                                 GParam const * const param,
                                 double *tildeG0,
                                 double *tildeGminp);
void perform_measures(Conf * GC,
                      Geometry const * const geo,
                      GParam const * const param,
                      FILE *datafilep);


// in gauge_conf_upd.c
void metropolis_for_link(Conf *GC,
                         Geometry const * const geo,
                         GParam const * const param,
                         long r,
                         int i);
void heatbath_for_link(Conf *GC,
                       Geometry const * const geo,
                       GParam const * const param,
                       long r,
                       int i);
void overrelaxation_for_link(Conf *GC,
                        Geometry const * const geo,
                        GParam const * const param,
                        long r,
                        int i);
void calcstaples_for_phi(Conf *GC,
                         Geometry const * const geo,
                         long r,
                         Vec *staple);

void overrelaxation_for_phi(Conf *GC,
                            Geometry const * const geo,
                            long r);
int metropolis_for_phi(Conf *GC,
                       Geometry const * const geo,
                       GParam const * const param,
                       long r);
void update(Conf * GC,
            Geometry const * const geo,
            GParam const * const param,
            double *acc);
int metropolis_for_phi_without_links(Conf *GC,
                                     Geometry const * const geo,
                                     GParam const * const param,
                                     long r);
void update_without_links(Conf * GC,
                          Geometry const * const geo,
                          GParam const * const param,
                          double *acc);


#endif
