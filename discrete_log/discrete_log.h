#ifndef DISCRETELOG_H
#define DISCRETELOG_H

#include <gmp.h>

void babystep_giantstep(mpz_t cyclic_group_order, mpz_t generator, mpz_t pow,
                        mpz_t log);

void pollard_rho(mpz_t cyclic_group_order, mpz_t generator, mpz_t pow,
                 mpz_t log);

#endif