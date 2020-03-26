bin_PROGRAMS = rpn

dSFMT = ./include/dSFMT-common.h \
        ./include/dSFMT.h \
        ./include/dSFMT-params11213.h \
        ./include/dSFMT-params1279.h \
        ./include/dSFMT-params132049.h \
        ./include/dSFMT-params19937.h \
        ./include/dSFMT-params216091.h \
        ./include/dSFMT-params2203.h \
        ./include/dSFMT-params4253.h \
        ./include/dSFMT-params44497.h \
        ./include/dSFMT-params521.h \
        ./include/dSFMT-params86243.h \
        ./include/dSFMT-params.h \
        ./lib/dSFMT.c

endian = ./include/endianness.h  \
         ./lib/endianness.c 

#flavourmatrix = ./include/flavour_matrix.h \
#                ./lib/flavour_matrix.c

conf = ./include/conf.h \
       ./lib/conf_def.c \
       ./lib/conf_meas.c 
#             ./lib/gauge_conf_meas_multilevel.c \
#             ./lib/gauge_conf_multilevel.c \
#             ./lib/gauge_conf_upd.c 

geometry = ./include/geometry.h \
           ./lib/geometry.c

gparam = ./include/gparam.h \
         ./lib/gparam.c

#son = ./include/son.h \
#      ./include/son_monopoles.h \
#      ./include/son_upd.h \
#      ./lib/son.c \
#      ./lib/son_monopoles.c \
#      ./lib/son_upd.c

random = ./include/myrngparam.h \
         ./include/random.h \
         ./lib/random.c 

vec = ./include/vec.h \
      ./lib/vec.c 

rpn_SOURCES=./src/rpn.c \
            $(conf) \
            $(dSFMT) \
            $(endian) \
            $(geometry) \
            $(gparam) \
            ./include/macro.h \
            $(random) \
            $(vec)

AM_CFLAGS = -std=c99 -O2 -Wall -Wextra -Werror -pedantic -Wconversion @OPENMP_CFLAGS@

CLEANFILES = dati.dat log.dat dati.log conf.dat* template_input.in 