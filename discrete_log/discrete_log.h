#ifndef DISCRETELOG_H
#define DISCRETELOG_H

#include <gmp.h>

void babystep_giantstep(mpz_t cyclic_group_order, mpz_t generator, mpz_t h,
                        mpz_t log);

#endif